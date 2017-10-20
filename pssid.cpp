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

typedef struct {
    std::vector<char> bases;
    std::vector<bool> strands;
    std::array<int, 4> counts;
} ReadStack;

ReadStack parseReadBases(const char* read, char reference, int coverage) {
    ReadStack result;
    result.bases.reserve(coverage);
    result.strands.reserve(coverage);
    result.counts = {0,0,0,0};

    for(size_t i = 0; i < strlen(read); ++i) {
        char base = read[i];
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
            case ',':
                result.strands.push_back(0);
                result.bases.push_back(toupper(reference)); break;
            case '.':
                result.strands.push_back(1);
                result.bases.push_back(toupper(reference)); break;
            default:
                break;
        }
    }
    return result;
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
    std::vector<bool> strands;
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
    auto stack = parseReadBases(read_bases, result.reference_base, coverage);
    result.bases = std::move(stack.bases);
    result.strands = std::move(stack.strands);
    result.base_counts = std::move(stack.counts);

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

long double inline logbinom(int n, int k) {
    return lgammal(n+1) - lgammal(k+1) - lgammal(n-k+1);
}

long double logalpha(int n, int k, long double e) {
    return logbinom(n, k) + k * logl(e) + (n-k) * logl(1L-e);
}

long double logbeta(int n, int k, long double e) {
    long double alpha_sum = 0;
    for(int i = k+1; i <= n; ++i) {
        alpha_sum += expl(logalpha(n, i, e));
    }
    return logl(alpha_sum / (alpha_sum + expl(logalpha(n, k, e))));
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

    std::string* outputs = new std::string[lines.size()];
#pragma omp parallel for default(none) shared(s, lines, outputs) schedule(dynamic, 100000)
    for (int i = 0; i < lines.size(); ++i) {
        const std::array<char,4> base_from_index {'A', 'C', 'G', 'T'};
        char* line = lines[i];
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
            auto& strands = plp.strands;

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
            int num_bases = bases.size();
            int num_hom_errors = 0;
            int num_het_errors = 0;
            int het_rank_forward = 0;
            int het_rank_backward = 0;
            int hom_rank_forward = 0;
            int hom_rank_backward = 0;
            for (size_t ii = 0; ii < num_bases; ++ii) {
                int i = base_ranks[ii];
                long double e = errors[i];
                long double de = dependent_errors[i];
                if (base_from_index[indices[0]] != bases[i]) {
                    if (strands[i]) {
                        logPhom += pow(0.85, hom_rank_forward) * logbeta(num_bases, num_hom_errors, e);
                        ++hom_rank_forward;
                    } else {
                        logPhom += pow(0.85, hom_rank_backward) * logbeta(num_bases, num_hom_errors, e);
                        ++hom_rank_backward;
                    }
                    ++num_hom_errors;
                }
                if (base_from_index[indices[0]] == bases[i]
                        || base_from_index[indices[1]] == bases[i]) {
                    // pass
                } else {
                    if (strands[i]) {
                        logPhet += pow(0.85, het_rank_forward) * logbeta(num_bases, num_het_errors, e);
                        ++het_rank_forward;
                    } else {
                        logPhet += pow(0.85, het_rank_backward) * logbeta(num_bases, num_het_errors, e);
                        ++het_rank_backward;
                    }
                    ++num_het_errors;
                }
            }

            auto logbinom = [](int n, int k) -> double {
               return lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1);
            };
            // distribution of "correct" bases
            int n = plp.base_counts[indices[0]] + plp.base_counts[indices[1]];
            int k = plp.base_counts[indices[1]];
            logPhet += logbinom(n, k) - n * M_LN2l;

            long double phom = exp(logPhom) * (1 - (long double)het_prior);
            long double phet = exp(logPhet) * (long double)het_prior;

            char buffer [256];
            sprintf(buffer, "%d,%d:%d:%d:%d,%s,%e,%e\n",
                    plp.pos,
                    plp.base_counts[0],
                    plp.base_counts[1],
                    plp.base_counts[2],
                    plp.base_counts[3],
                    phet > phom ? "het" : "hom",
                    double(phom),
                    double(phet));
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
