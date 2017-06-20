#include <algorithm>
#include <array>
#include <cmath>
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

struct arguments {
    double p_value_threshold {0.05};
} args {};

void processFile(FILE* input) {
    char name [NAME_BUFFER_SIZE];
    int pos;
    char reference_base;
    int coverage;
    char read [READ_BUFFER_SIZE];

    vector<pair<int, Profile>> positions;
    vector<Profile> profiles;
    int count = 0;
    while(!feof(input) && !ferror(input)) {
        int parsed_fields = fscanf(input, "%s %d %c %d %s %*s\n", &name[0], &pos, &reference_base, &coverage, &read[0]);
        if (ferror(input) || parsed_fields < 5) {
            cerr << "Read error: input must be valid 'samtools mpileup' output." << endl;
            exit(EXIT_FAILURE);
        }

        Profile p = parseRead(read, reference_base);
        if (p[COV] >= 4) {
            // mapping position -> profile
            positions.emplace_back(pos, p);
            profiles.push_back(p);
        }
        ++count;
    }

    if(positions.size() == 0) {
        cerr << "No profiles found with required minimum coverage of 4!" << endl;
        exit(EXIT_FAILURE);
    }

    sort(profiles.begin(), profiles.end());
    vector<int> counts {1};
    // count number of occurences of each profile
    for(auto it = profiles.begin() + 1; it != profiles.end(); ++it) {
        if (*it != *(it - 1)) {
            counts.push_back(1);
        } else {
            ++counts.back();
        }
    }
    // remove duplicate profiles
    auto new_end = unique(profiles.begin(), profiles.end());
    profiles.resize(distance(profiles.begin(), new_end));

    cerr << "# " << positions.size() << " of " << count << " reads with required coverage > 4, ";
    cerr << counts.size() << " distinct profiles" << endl;

    map<Profile, int> index_of {};
    for (int i = 0; i < profiles.size(); ++i) {
        index_of.emplace(profiles[i], i);
    }

    vector<pair<double, double>> profile_likelihoods = computeLikelihoods(profiles, counts);
    vector<double> p_values = likelihoodRatioTest(profile_likelihoods);
    vector<double> p_values_bonf = adjustBonferroni(p_values, positions.size());
    vector<double> p_values_bh = adjustBenjaminiHochberg(p_values);

    vector<size_t> sorted = descending_sorted_indices(p_values);

    cout << "# pos\tprofile\t\tl_ho\t\tl_het\t\tp\t\tp_bonf\t\tp_bh" << endl;
    cout << scientific;
    cout.precision(2);
    for (auto& [pos, profile] : positions) {
        if (p_values_bh[index_of[profile]] <= args.p_value_threshold || isnan(p_values[index_of[profile]])) {
            cout << pos << '\t' << profile;
            cout << '\t' << profile_likelihoods[index_of[profile]].first << '\t' << profile_likelihoods[index_of[profile]].second;
            cout << '\t' << p_values[index_of[profile]];
            cout << '\t' << p_values_bonf[index_of[profile]];
            cout << '\t' << p_values_bh[index_of[profile]];
            cout << endl;
        }
    }
}

void printHelp(char* program_name) {
    cout << "Usage: " << program_name << " [options] [input files]" << endl;
    cout << "Options:" << endl;
    vector<pair<string, string>> options = {
        {"-h", "print this help message and exit"},
        {"-p NUM", "p value threshold"},
        // {"-o FILE", {"print output to FILE"}
    };
    for (const auto& [option, description] : options) {
        cout << '\t' << option << "\t\t\t" << description << endl;
    }
}

int main(int argc, char** argv) {
    string optstring = "hp:";
    for (char opt = getopt(argc, argv, optstring.c_str()); opt != -1; opt = getopt(argc, argv, optstring.c_str())) {
        switch (opt) {
        case 'h':
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        case 'p':
            double p_value = atof(optarg);
            if (p_value <= 0.0) {
                cerr << "Invalid p value threshold or parse error: " << optarg << endl;
                exit(EXIT_FAILURE);
            } else {
                args.p_value_threshold = p_value;
            }
            break;
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
