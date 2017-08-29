#include <numeric>
#include <limits>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "profiles.hpp"

using namespace std;

bool inline tryIncrementBaseCount(char base, Profile& p) {
    switch (base) {
        case 'a':
        case 'A':
            ++p[A]; break;
        case 'c':
        case 'C':
            ++p[C]; break;
        case 'g':
        case 'G':
            ++p[G]; break;
        case 't':
        case 'T':
            ++p[T]; break;
        default:
            return false;
    }
    return true;
}

Profile parseReadBases(const char* read, char reference) {
    Profile p {0, 0, 0, 0, 0};

    for(size_t i = 0; i < strlen(read); ++i) {
        char base = read[i];
        if (!tryIncrementBaseCount(base, p)) {
            switch (base) {
                case '^':
                    // skip next char
                    ++i;
                    break;
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
                    if (numeric_limits<decltype(i)>::max() - length < i) {
                        i = numeric_limits<decltype(i)>::max();
                    } else {
                        // skip parsed number + that number of bases after the number
                        // (-1 because i is incremented in the surrounding loop)
                        i = (first_after_number - read) + length - 1;
                    }
                    break;
                }
                case '.':
                case ',': {
                    tryIncrementBaseCount(reference, p);
                    break;
                }
                default:
                    break;
            }
        }
    }
    auto sum = p[A] + p[C] + p[G] + p[T];
    if (sum > numeric_limits<decltype(p)::value_type>::max()) {
        cerr << "Overflow error during base parsing!";
        exit(EXIT_FAILURE);
    } else {
        p[COV] = (decltype(p)::value_type)sum;
        return p;
    }
}

PileupLine parsePileupLine(char* line) {
    PileupLine result;
    const char* delim = " \t";

    // drop name
    strtok(line, delim);

    // read position
    char* pos = strtok(NULL, delim);
    if (pos == NULL) {
        return result;
    }
    result.pos = atoi(pos);

    // read reference base
    char* ref = strtok(NULL, delim);
    if (ref == NULL) {
        result.pos = -1;
        return result;
    }
    result.reference_base = ref[0];

    // drop coverage
    strtok(NULL, delim);

    // parse nucleotide profile
    char* read_bases = strtok(NULL, delim);
    if (read_bases == NULL) {
        result.pos = -1;
        return result;
    }
    result.profile = parseReadBases(read_bases, result.reference_base);
    return result;
}

PileupData readPileup(FILE* input) {
    PileupData result;

    size_t line_length = 10000;
    char* line = (char*)malloc(line_length*sizeof(char));

    int lines = 0;
    while(getline(&line, &line_length, input) != -1) {
        PileupLine plpline = parsePileupLine(line);
        if (plpline.pos < 0) {
            cerr << "# Read error: input must be valid 'samtools mpileup' output." << endl;
            exit(EXIT_FAILURE);
        }

        if (plpline.profile[COV] >= 4) {
            result.positions.push_back(plpline.pos);
            result.profiles.push_back(plpline.profile);
            result.num_sites += 1;
        }
        ++lines;
    }
    free(line);

    if(result.num_sites == 0) {
        cerr << "# No profiles found with required minimum coverage of 4!" << endl;
        exit(EXIT_FAILURE);
    }
    cerr << "# " << result.num_sites << " of " << lines << " sites with required coverage > 4" << endl;
    return result;
}

ostream& std::operator<<(ostream& os, const Profile& p) {
    return os << p[0] << "." << p[1] << "." << p[2] << "." << p[3];
}
