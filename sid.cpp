#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <getopt.h>

#include "profiles.hpp"
#include "likelihoods.hpp"
#include "stats.hpp"

using namespace std;

const int NAME_BUFFER_SIZE = 256;
const int READ_BUFFER_SIZE = 8000;

void processFile(FILE* input) {
    char name [NAME_BUFFER_SIZE];
    int pos;
    char reference_base;
    int coverage;
    char read [READ_BUFFER_SIZE];

    map<int, Profile> positions;
    vector<Profile> profiles;
    map<Profile, int> profile_counts;

    while(!feof(input) && !ferror(input)) {
        int parsed_fields = fscanf(input, "%s %d %c %d %s %*s\n", &name[0], &pos, &reference_base, &coverage, &read[0]);
        if (ferror(input) || parsed_fields < 5) {
            cerr << "Read error: input must be valid 'samtools mpileup' output." << endl;
            exit(EXIT_FAILURE);
        }

        Profile p = parseRead(read, reference_base);
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


    vector<pair<double, double>> profile_likelihoods = computeLikelihoods(profile_counts);
    vector<double> p_values = likelihoodRatioTest(profile_likelihoods);
    vector<double> p_values_bonf = adjustBonferroni(p_values);
    vector<double> p_values_bh = adjustBenjaminiHochberg(p_values);

    vector<size_t> sorted = descending_sorted_indices(p_values);
    vector<pair<Profile, int>> profiles2 (profile_counts.begin(), profile_counts.end());
    cout << "# Profile\t\tcount\tl_ho\tl_het\tp\tp_bonf\tp_bh" << endl;
    cout.precision(2);
    for (int i = 0; i < profiles2.size(); ++i) {
        if (p_values[sorted[i]] <= 0.05) {
            //cout << profile << '\t' << count;
            cout << profiles2[sorted[i]].first << "\t\t" << profiles2[sorted[i]].second;
            cout << '\t' << profile_likelihoods[sorted[i]].first << '\t' << profile_likelihoods[sorted[i]].second;
            cout << '\t' << p_values[sorted[i]];
            cout << '\t' << p_values_bonf[sorted[i]];
            cout << '\t' << p_values_bh[sorted[i]];
            cout << endl;
        }
    }
}

void printHelp(char* program_name) {
    cout << "Usage: " << program_name << " [options] [input files]" << endl;
    cout << "Options:" << endl;
    vector<pair<string, string>> options = {
        {{"-h"}, {"print this help message and exit"}},
        {{"-o FILE"}, {"print output to FILE"}}
    };
    for (const auto& [option, description] : options) {
        cout << '\t' << option << "\t\t\t" << description << endl;
    }
}

int main(int argc, char** argv) {
    string optstring = "ho:";
    for (char opt = getopt(argc, argv, optstring.c_str()); opt != -1; opt = getopt(argc, argv, optstring.c_str())) {
        switch (opt) {
        case 'h':
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
        }
    }
    if (optind < argc) {
        FILE* f = fopen(argv[optind], "r");
        if (f != nullptr) {
            processFile(f);
        } else {
            cerr << "Error opening file: " << argv[optind] << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        processFile(stdin);
    }
}
