#include <cctype>
#include <stdexcept>

PileupLine parsePileupLine(char* line, bool parse_base_qualities, bool parse_mapping_qualities) {
    const string MALFORMED = "Malformed pileup line";
    const string MALFORMED_OR_MISSING = "Malformed pileup line or missing mapping qualities";
    const char* DELIM = " \t";
    char* saveptr;

    PileupLine result;

    // chromosome name
    char* chrom = strtok_r(line, DELIM, &saveptr);
    result.chrom = chrom;

    // read position
    char* pos = strtok_r(nullptr, DELIM, &saveptr);
    if (pos == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    result.pos = atoi(pos);

    // read reference base
    char* ref = strtok_r(nullptr, DELIM, &saveptr);
    if (ref == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    result.reference_base = ref[0];

    // read coverage
    char* cov = strtok_r(nullptr, DELIM, &saveptr);
    if (cov == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    int coverage = atoi(cov);

    // parse bases
    char* read_bases = strtok_r(nullptr, DELIM, &saveptr);
    if (read_bases == nullptr) {
        throw std::invalid_argument {MALFORMED};
    }
    ReadStack stack = parseReadBases(read_bases, result.reference_base, coverage);
    result.bases = std::move(stack.bases);
    result.strands = std::move(stack.strands);
    result.base_counts = std::move(stack.counts);

    // ensure pointer is advanced if only mapping qualities are to be parsed
    char* base_qualities = strtok_r(nullptr, DELIM, &saveptr);

    // parse base qualities if necessary
    if (parse_base_qualities) {
        if (read_bases == nullptr) {
            throw std::invalid_argument {MALFORMED};
        }
        result.base_qualities = parseQualities(base_qualities, coverage);
    }

    // parse mapping qualities if necessary
    if (parse_mapping_qualities) {
        char* mapping_qualities = strtok_r(nullptr, DELIM, &saveptr);
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
        if(std::isupper(base)) {
            result.strands.push_back(1);
        } else if (std::islower(base)) {
            result.strands.push_back(0);
        }
        switch (base) {
            case 'a':
            case 'A':
                result.bases.push_back('A');
                ++result.counts[0];
                break;
            case 'c':
            case 'C':
                result.bases.push_back('C');
                ++result.counts[1];
                break;
            case 'g':
            case 'G':
                result.bases.push_back('G');
                ++result.counts[2];
                break;
            case 't':
            case 'T':
                result.bases.push_back('T');
                ++result.counts[3];
                break;
            case ',':
                result.strands.push_back(0);
                result.bases.push_back(toupper(reference)); break;
            case '.':
                result.strands.push_back(1);
                result.bases.push_back(toupper(reference)); break;
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
    for (char* q = base_qualities; *q != '\0' && *q != '\t' && *q != '\n'; ++q) {
        uint8_t quality = uint8_t(*q) - 33;
        // FIXME lower threshold for qualities necessary? -> adapt processing code if 0 is allowed
        if (quality < 1) {
            quality = 1;
        }
        result.push_back(quality);
    }
    return result;
}
