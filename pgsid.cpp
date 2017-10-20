#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <limits>
#include <utility>
#include <string>
#include <vector>

#include "stats.hpp"

std::pair<std::vector<char>, std::array<int, 4>> parseReadBases(const char* read, char reference, int coverage) {
    std::vector<char> result;
    std::array<int, 4> counts {0,0,0,0};
    result.reserve(coverage);

    for(size_t i = 0; i < strlen(read); ++i) {
        char base = read[i];
        switch (base) {
            case 'a':
            case 'A':
                result.push_back('A'); ++counts[0]; break;
            case 'c':
            case 'C':
                result.push_back('C'); ++counts[1]; break;
            case 'g':
            case 'G':
                result.push_back('G'); ++counts[2]; break;
            case 't':
            case 'T':
                result.push_back('T'); ++counts[3]; break;
            case '^':
                // skip next char
                ++i; break;
            case '+':
            case '-': {
                // parse following number, which indicates a range of insert/del bases
                if (!isdigit(read[i + 1])) {
                    break;
                }
                char* first_after_number;
                // number is always positive since '-' is skipped
                unsigned long length = (unsigned long)strtol(read + i + 1, &first_after_number, 10);

                // overflow handled manually, surrounding for loop will then terminate
                if (std::numeric_limits<decltype(i)>::max() - length < i) {
                    i = std::numeric_limits<decltype(i)>::max();
                } else {
                    // skip parsed number + that number of bases after the number
                    // (-1 because i is incremented in the surrounding loop)
                    i = (first_after_number - read) + length - 1;
                }
                break;
            }
            case '.':
            case ',':
                result.push_back(toupper(reference)); break;
            default:
                break;
        }
    }
    return {result, counts};
}

std::vector<uint8_t> parseQualities(char* base_qualities, int coverage) {
    std::vector<uint8_t> result;
    result.reserve(coverage);
    for (char* q = base_qualities; *q != '\0' && *q != '\t' && *q != '\n'; ++q) {
        uint8_t quality = uint8_t(*q) - 33;
        if (quality < 1) {
            quality = 1;
        }
        result.push_back(quality);
    }
    return result;
}

typedef struct {
    int pos {-1};
    char reference_base {'N'};
    std::array<int, 4> base_counts;
    std::vector<char> bases;
    std::vector<uint8_t> base_qualities;
    std::vector<uint8_t> mapping_qualities;
} PileupLine;

PileupLine parsePileupLine(char* line) {
    PileupLine result;
    const char* DELIM = " \t";
    char* saveptr;

    // drop name
    strtok_r(line, DELIM, &saveptr);

    // read position
    char* pos = strtok_r(nullptr, DELIM, &saveptr);
    if (pos == nullptr) {
        return result;
    }
    result.pos = atoi(pos);

    // read reference base
    char* ref = strtok_r(nullptr, DELIM, &saveptr);
    if (ref == nullptr) {
        result.pos = -1;
        return result;
    }
    result.reference_base = ref[0];

    // read coverage
    char* cov = strtok_r(nullptr, DELIM, &saveptr);
    if (cov == nullptr) {
        result.pos = -1;
        return result;
    }
    int coverage = atoi(cov);

    // parse bases
    char* read_bases = strtok_r(nullptr, DELIM, &saveptr);
    if (read_bases == nullptr) {
        result.pos = -1;
        return result;
    }
    auto parsed = parseReadBases(read_bases, result.reference_base, coverage);
    result.bases = parsed.first;
    result.base_counts = parsed.second;

    // parse base qualities
    char* base_qualities = strtok_r(nullptr, DELIM, &saveptr);
    if (read_bases == nullptr) {
        result.pos = -1;
        return result;
    }
    result.base_qualities = parseQualities(base_qualities, coverage);

    char* mapping_qualities = strtok_r(nullptr, DELIM, &saveptr);
    if (mapping_qualities == nullptr) {
        std::cerr << "Missing mapping qualities!" << std::endl;
        result.pos = -1;
        return result;
    }
    result.mapping_qualities = parseQualities(mapping_qualities, coverage);
    return result;
}

/* void processFile(FILE* input, double het_prior) { */
void processFile(std::istream& in, const double het_prior) {

    auto const start_pos = in.tellg();
    in.ignore(std::numeric_limits<std::streamsize>::max());
    auto const char_count = in.gcount();
    in.seekg(start_pos);
    /* auto s = std::string(char_count, char{}); */
    char* s = new char[char_count+1];
    /* in.read(&s[0], s.size()); */
    in.read(&s[0], char_count);
    s[char_count] = '\0';

    std::vector<char*> lines;
    char* saveptr;
    char* line = strtok_r(s, "\n", &saveptr);
    lines.push_back(line);
    while((line = strtok_r(nullptr, "\n", &saveptr)) != nullptr) {
        lines.push_back(line);
    }

    /* std::vector<char*> lines {nullptr}; */
    /* size_t line_length = 0; */

    /* while(getline(&lines.back(), &line_length, input) != -1) { */
    /*     if (line_length > 0) { */
    /*         lines.push_back(NULL); */
    /*         line_length = 0; */
    /*     } */
    /* } */
    /* std::vector<std::string> outputs (lines.size()); */
    std::string* outputs = new std::string[lines.size()];
    /* unsigned int n_lines = 0; */
#pragma omp parallel for default(none) shared(s, lines, outputs)
    for (int i = 0; i < lines.size(); ++i) {
        const std::array<char,4> base_from_index {'A', 'C', 'G', 'T'};
        char* line = lines[i];
    /* for (char* line : lines) { */
        PileupLine plp = parsePileupLine(line);
        /* ::free(line); */
        if (plp.pos < 0) {
            outputs[i] = "";
        } else {
            std::array<int, 4> indices {0,1,2,3};
            std::sort(indices.begin(), indices.end(),
                    [&counts = plp.base_counts](int a, int b) {
                        return counts[a] > counts[b];
                    });
            auto& bases = plp.bases;
            auto& base_qualities = plp.base_qualities;
            auto& mapping_qualities = plp.mapping_qualities;

            std::vector<long double> errors;
            errors.reserve(bases.size());
            for (size_t i = 0; i < bases.size(); ++i) {
                long double base_e = pow(10L, -(long double)(base_qualities[i]) / 10L);
                long double mapping_e = pow(10L, -(long double)(mapping_qualities[i]) / 10L);
                long double e = std::max(base_e, mapping_e);
                errors.push_back(e);
            }
            std::vector<size_t> base_ranks(bases.size());
            int j = 0;
            std::generate(base_ranks.begin(), base_ranks.end(),
                    [&j] () { return j++; });
            std::sort(base_ranks.begin(), base_ranks.end(),
                    [&errors](auto a, auto b) {
                        return errors[a] < errors[b];
                    });
            std::vector<long double> dependent_errors(errors.size());
            std::transform(errors.begin(), errors.end(),
                    base_ranks.begin(), dependent_errors.begin(),
                    [](long double e, size_t rank) -> long double {
                        return powl(e, pow(0.5, rank));
                    });

            long double logPhom = 0;
            long double logPhet = 0;
            for (size_t i = 0; i < bases.size(); ++i) {
                long double e = errors[i];
                long double de = dependent_errors[i];
                if (base_from_index[indices[0]] == bases[i]) {
                    logPhom += log(1L - e);
                } else {
                    /* logPhom += log(de); */
                    logPhom += log(e);
                }
                if (base_from_index[indices[0]] == bases[i]
                        || base_from_index[indices[1]] == bases[i]) {
                    /* logPhet += std::max(log(0.5L - e/3.L), 0.0); */
                    logPhet += log(1L - 2.L * e/3.L);
                } else {
                    /* logPhet += log(2.L/3.L * de); */
                    logPhet += log(2.L/3.L * e);
                }
            }
            /* int n1 = plp.base_counts[indices[1]] + plp.base_counts[indices[2]] + plp.base_counts[indices[3]]; */
            /* logPhom += lgamma(n1 + 1); */
            /* logPhom -= lgamma(plp.base_counts[indices[1]] + 1); */
            /* logPhom -= lgamma(plp.base_counts[indices[2]] + 1); */
            /* logPhom -= lgamma(plp.base_counts[indices[3]] + 1); */
            /* logPhom -= n1 * log(3); */

            auto logbinom = [](int n, int k) -> double {
               return lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1);
            };
            // distribution of "correct" bases
            int n = plp.base_counts[indices[0]] + plp.base_counts[indices[1]];
            int k = plp.base_counts[indices[1]];
            logPhet += logbinom(n, k) - n * M_LN2l;

            // distribution of "error" bases
            /* int n2 = plp.base_counts[indices[2]] + plp.base_counts[indices[3]]; */
            /* int k2 = plp.base_counts[indices[2]]; */
            /* logPhet += logbinom(n2, k2) - n2 * M_LN2l; */

            long double phom = exp(logPhom) * (1 - (long double)het_prior);
            long double phet = exp(logPhet) * (long double)het_prior;
            /* phom = phom / (phom + phet); */
            /* phet = phet / (phom + phet); */

            double pval_hom = likelihoodRatioTest(phet, phom);
            double pval_het = likelihoodRatioTest(phom, phet);

            char buffer [256];
            sprintf(buffer, "%d,%d:%d:%d:%d,%s,%e,%e\n",
                    plp.pos,
                    plp.base_counts[0],
                    plp.base_counts[1],
                    plp.base_counts[2],
                    plp.base_counts[3],
                    pval_het < 0.00001 ? "het" : "hom",
                    pval_hom,
                    pval_het);
            /* sprintf(buffer, "%d,%d:%d:%d:%d,%s,%Le,%Le\n", */
            /*         plp.pos, */
            /*         plp.base_counts[0], */
            /*         plp.base_counts[1], */
            /*         plp.base_counts[2], */
            /*         plp.base_counts[3], */
            /*         phet > phom ? "het" : "hom", */
            /*         phom, */
            /*         phet); */
            outputs[i] = std::string(buffer);
            /* ++n_lines; */
        }
    }
    /* for (auto o : outputs) { */
    for (size_t i = 0; i < lines.size(); ++i) {
        std::cout << outputs[i];
    }
    delete[] outputs;
    delete[] s;
    /* std::cout << "Processed " << n_lines << " lines\n"; */
}

int main(int argc, char** argv) {
    double het_prior {1};
    if (argc > 1) {
        het_prior = atof(argv[1]);
        if (het_prior == 0.0) {
            std::cerr << "# Invalid heterozygous prior: " << het_prior << " (parsed from " << argv[1] << ')' << std::endl;
            ::exit(EXIT_FAILURE);
        }
    }
    if (argc > 2) {
        /* FILE* f = ::fopen(argv[2], "r"); */
        std::ifstream f {argv[2]};
        /* if (f != nullptr) { */
        if (f) {
            processFile(f, het_prior);
        } else {
            std::cerr << "# Error opening file: " << argv[2] << std::endl;
        }
    } else {
        /* processFile(stdin, het_prior); */
        processFile(std::cin, het_prior);
    }
}
