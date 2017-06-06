#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "profiles.hpp"
#include "likelihoods.hpp"

using namespace std;

const int NAME_BUFFER_SIZE = 256;
const int READ_BUFFER_SIZE = 8000;

int main() {
    char name [NAME_BUFFER_SIZE];
    int pos;
    char reference_base;
    int coverage;
    char read [READ_BUFFER_SIZE];

    map<int, Profile> positions;
    vector<Profile> profiles;
    map<Profile, int> profile_counts;

    while(!feof(stdin) && !ferror(stdin)) {
        int parsed_fields = fscanf(stdin, "%s %d %c %d %s %*s\n", &name[0], &pos, &reference_base, &coverage, &read[0]);
        if (ferror(stdin) || parsed_fields < 5) {
            cerr << "Read error: input must be valid 'samtools mpileup' output." << endl;
            exit(EXIT_FAILURE);
        }

        Profile p = parseRead((read));
        if (p[COV] >= 4) {
            positions.insert(make_pair(pos, p));
            profiles.push_back(p);

            bool new_element;
            map<Profile, int>::iterator position;
            tie(position, new_element) = profile_counts.emplace(p, 1);
            if (!new_element) {
                position->second += 1;
            }
        }
    }

    cerr << "Parsed " << positions.size() << " reads" << endl;
    cerr << "Found " << profile_counts.size() << " distict profiles" << endl;

    computeLikelihoods(profile_counts);
}
