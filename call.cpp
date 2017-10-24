#include <algorithm>
#include <cstring>
#include <limits>

#include "call.hpp"

std::vector<PileupLine> readFile(std::istream& in, const bool parse_base_qualities, const bool parse_mapping_qualities) {
    std::vector<PileupLine> result;

    // read whole input file into buffer
    const auto start_pos = in.tellg();
    in.ignore(std::numeric_limits<std::streamsize>::max());
    const auto char_count = in.gcount();
    in.seekg(start_pos);
    char* s = new char[char_count+1];
    in.read(&s[0], char_count);
    s[char_count] = '\0';

    std::vector<char*> lines;
    char* saveptr;
    char* line = strtok_r(s, "\n", &saveptr);
    lines.push_back(line);
    while((line = strtok_r(nullptr, "\n", &saveptr)) != nullptr) {
        lines.push_back(line);
    }
    result.reserve(lines.size());

    #pragma omp parallel for default(none) shared(s, lines, result)
    for (size_t i = 0; i < lines.size(); ++i) {
        PileupLine plp = parsePileupLine(lines[i], parse_base_qualities, parse_mapping_qualities);
        result.push_back(std::move(plp));
    }
    delete[] s;
    return result;
}

// std::vector<PileupLine> parseLines(const std::vector<char*>& lines, const bool parse_base_qualities, const bool parse_mapping_qualities) {
//     std::vector<PileupLine> result;
//     result.reserve(lines.size());

//     #pragma omp parallel for default(none) shared(lines, lines, result)
//     for (size_t i = 0; i < lines.size(); ++i) {
//         PileupLine plp = parsePileupLine(lines[i], parse_base_qualities, parse_mapping_qualities);
//         result.push_back(std::move(plp));
//     }
//     delete[] s;
//     return result;
// }

std::vector<UniqueProfile> countUniqueProfiles(const std::vector<PileupLine>& pileup) {
    if (pileup.size() == 0) {
        return {};
    }
    using profile = std::array<int, 4>;
    std::vector<const profile*> profiles;
    profiles.reserve(pileup.size());
    std::transform(pileup.begin(), pileup.end(), std::back_inserter(profiles),
        [](const PileupLine& line) {
            return &line.base_counts;
        });
    std::sort(profiles.begin(), profiles.end(),
        [](const profile* a, const profile* b) {
            return *a < *b;
        });

    std::vector<UniqueProfile> unique_profiles;
    unique_profiles.reserve(pileup.size());
    // initialize with first profile, will be incremented in first loop pass
    unique_profiles.push_back({*profiles[0], 0});
    for (const profile* p : profiles) {
        if (*p == unique_profiles.back().profile) {
            unique_profiles.back().count += 1;
        } else {
            unique_profiles.push_back({*p, 1});
        }
    }

    return unique_profiles;
}

// std::vector<VCFRecord> callLikelihoodRatio(std::istream& in) {

// }

std::vector<OutputRecord> callBayes(std::istream& in) {
    auto inputRecords = readFile(in);
    auto unique_profiles = countUniqueProfiles(inputRecords);

    // estimate genotype likelihoods

    // build reverse index

    // calculate genotype likelihoods

    // build output lines
    std::vector<OutputRecord> outputRecords;
    outputRecords.reserve(inputRecords.size());

    return outputRecords;
}

// std::vector<VCFRecord> callSiteMLError(std::istream& in, bool use_prior) {

// }

// std::vector<VCFRecord> callQualityBasedSimple(std::istream& in) {

// }

// std::vector<VCFRecord> callQualityBasedSamtools(std::istream& in) {

// }
