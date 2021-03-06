#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <stdexcept>

#include "pileup.hpp"

const std::string MALFORMED = "Malformed pileup line";
const std::string MALFORMED_OR_MISSING = "Malformed pileup line or missing mapping qualities";
const char* FIELD_SEPARATORS = " \t";

PileupLine parsePileupLine(char* line, bool parse_base_qualities, bool parse_mapping_qualities) {
    char* saveptr = nullptr;
    PileupLine result;

    char* chromosome_name = strtok_r(line, FIELD_SEPARATORS, &saveptr);
    result.chromosome_name = chromosome_name;

    char* position = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);
    if (position == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    result.position = atoi(position);

    char* reference = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);
    if (reference == nullptr || strlen(reference) != 1) {
        throw std::invalid_argument {MALFORMED};
    }
    result.reference_base = reference[0];

    char* coverage_str = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);
    if (coverage_str == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    int coverage = atoi(coverage_str);

    char* read_bases = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);
    if (read_bases == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    ReadStack stack = parseReadBases(read_bases, result.reference_base, coverage);
    result.bases = std::move(stack.bases);
    result.strands = std::move(stack.strands);
    result.base_counts = std::move(stack.counts);

    // ensure pointer is advanced if only mapping qualities are to be parsed
    char* base_qualities = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);

    // parse base qualities if necessary
    if (parse_base_qualities) {
        if (read_bases == nullptr) {
            throw std::invalid_argument {MALFORMED};
        }
        result.base_qualities = parseQualities(base_qualities, coverage);
    }

    // parse mapping qualities if necessary
    if (parse_mapping_qualities) {
        char* mapping_qualities = strtok_r(nullptr, FIELD_SEPARATORS, &saveptr);
        if (mapping_qualities == nullptr) {
            throw std::invalid_argument {MALFORMED_OR_MISSING};
        }
        result.mapping_qualities = parseQualities(mapping_qualities, coverage);
    }

    return result;
}

ReadStack parseReadBases(const char* read_bases, char reference, int coverage) {
    ReadStack result;
    result.bases.reserve(coverage);
    result.strands.reserve(coverage);
    result.counts = {0,0,0,0};

    for(size_t i = 0; i < strlen(read_bases); ++i) {
        char base = read_bases[i];
        if (base == '.') {
            base = char(std::toupper(reference));
        }
        else if (base == ',') {
            base = char(std::tolower(reference));
        }
        switch (base) {
            case 'a':
                result.bases.push_back('A');
                result.strands.push_back(0);
                ++result.counts[0];
                break;
            case 'A':
                result.bases.push_back('A');
                result.strands.push_back(1);
                ++result.counts[0];
                break;
            case 'c':
                result.bases.push_back('C');
                result.strands.push_back(0);
                ++result.counts[1];
                break;
            case 'C':
                result.bases.push_back('C');
                result.strands.push_back(1);
                ++result.counts[1];
                break;
            case 'g':
                result.bases.push_back('G');
                result.strands.push_back(0);
                ++result.counts[2];
                break;
            case 'G':
                result.bases.push_back('G');
                result.strands.push_back(1);
                ++result.counts[2];
                break;
            case 't':
                result.bases.push_back('T');
                result.strands.push_back(0);
                ++result.counts[3];
                break;
            case 'T':
                result.bases.push_back('T');
                result.strands.push_back(1);
                ++result.counts[3];
                break;
            case '^':
                // skip next char
                ++i; break;
            case '+':
            case '-': {
                // parse following number, which indicates a range of insert/del bases
                if (!isdigit(read_bases[i + 1])) {
                    break;
                }
                char* first_after_number;
                // number is always positive since '-' is skipped
                unsigned long length = (unsigned long)strtol(read_bases + i + 1, &first_after_number, 10);

                // overflow handled manually, surrounding for loop will then terminate
                if (std::numeric_limits<decltype(i)>::max() - length < i) {
                    i = std::numeric_limits<decltype(i)>::max();
                } else {
                    // skip parsed number + that number of bases after the number
                    // (-1 because i is incremented in the surrounding loop)
                    i = (first_after_number - read_bases) + length - 1;
                }
                break;
            }
            default:
                break;
        }
    }
    return result;
}

std::vector<uint8_t> parseQualities(const char* base_qualities, int coverage) {
    std::vector<uint8_t> result;
    result.reserve(coverage);
    for (const char* q = base_qualities; *q != '\0' && *q != '\t' && *q != '\n'; ++q) {
        uint8_t quality = uint8_t(*q - 33);
        // FIXME lower threshold for qualities necessary? -> adapt processing code if 0 is allowed
        if (quality < 1) {
            quality = 1;
        }
        result.push_back(quality);
    }
    return result;
}

std::vector<UniqueProfile> countUniqueProfiles(const std::vector<PileupLine>& pileup) {
    if (pileup.size() == 0) {
        return {};
    }
    std::vector<const profile_t*> profiles;
    profiles.reserve(pileup.size());
    std::transform(pileup.begin(), pileup.end(), std::back_inserter(profiles),
        [](const PileupLine& line) {
            return &line.base_counts;
        });
    std::sort(profiles.begin(), profiles.end(),
        [](const profile_t* a, const profile_t* b) {
            return *a < *b;
        });

    std::vector<UniqueProfile> unique_profiles;
    // initialize with first profile, will be incremented in first loop pass
    unique_profiles.emplace_back(*profiles[0], 0);
    for (const auto p : profiles) {
        if (*p == unique_profiles.back().profile) {
            unique_profiles.back().count += 1;
        } else {
            unique_profiles.emplace_back(*p, 1);
        }
    }

    return unique_profiles;
}

std::array<double, 4> computeNucleotideDistribution(const std::vector<UniqueProfile>& profiles) {
    std::array<uint64_t, 4> acc {0,0,0,0};
    uint64_t total = 0;
    for (UniqueProfile p : profiles) {
        total += p.count * p.coverage;
        for (int i = 0; i < 4; ++i) {
            acc[i] += p.count * p.profile[i];
        }
    }

    if (total != 0) {
        return {
            double(acc[0]) / double(total),
            double(acc[1]) / double(total),
            double(acc[2]) / double(total),
            double(acc[3]) / double(total)};
    } else {
        return {0.25,0.25,0.25,0.25};
    }
}
