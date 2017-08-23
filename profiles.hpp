#ifndef PROFILES_H
#define PROFILES_H

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

/** Profile of a pileup site.
 *
 * The values in the array represent the number of occurences of each
 * nucleotide and their sum (the site coverage):
 * [#A, #C, #G, #T, #(A+C+G+T)]
 */
using Profile = std::array<uint16_t, 5>;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;
const int COV = 4;

namespace std {
    std::ostream& operator<<(std::ostream& os, const Profile& p);
}

/** Compute site nucleotide profile from read bases.
 *
 * Given the read base column of a line of `samtools mpileup` output,
 * compute the site nucleotide profile, while ignoring inserts and
 * deletions.  See the samtools manpage (subcommand mpileup) for the
 * format.
 */
Profile parseReadBases(const char* read_bases, char reference);

/** Alignment position, reference base and samtools mpileup-encoded read bases.
 */
typedef struct {
    int pos {-1};
    char reference_base {'N'};
    Profile profile;
} PileupLine;

/** Extract the relevant information (for us) from a line of samtools
 * mpileup output.
 */
PileupLine parsePileupLine(char* line);

/** The data read from `samtools mpileup` output.
 *
 * Contains the number of sites read and two vectors with that number
 * of elements, one containing the profiles and one their alignment
 * positions. Only sites with coverage >=4 are selected.
 */
typedef struct {
    int num_sites {0};
    std::vector<int> positions {};
    std::vector<Profile> profiles {};
} PileupData;

/** Read samtools mpileup output from given file.
 *  The file handle must not be NULL.
 */
PileupData readPileup(FILE* input);

#endif
