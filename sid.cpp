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

enum class Correction {Bonferroni, BH, None};

struct arguments {
    Correction correction {Correction::Bonferroni};
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

    vector<double> p_het (profile_likelihoods.size());
    transform(profile_likelihoods.begin(), profile_likelihoods.end(), p_het.begin(),
              [](pair<double, double> ls) {
                  return likelihoodRatioTest(ls.first, ls.second);
              });

    vector<double> p_hom (profile_likelihoods.size());
    transform(profile_likelihoods.begin(), profile_likelihoods.end(), p_hom.begin(),
              [](pair<double, double> ls) {
                  return likelihoodRatioTest(ls.second, ls.first);
              });
    vector<double> p_het_adj;
    vector<double> p_hom_adj;
    if (args.correction == Correction::Bonferroni) {
        p_het_adj = adjustBonferroni(p_het, positions.size());
        p_hom_adj = adjustBonferroni(p_hom, positions.size());
    } else if (args.correction == Correction::BH) {
        p_het_adj = adjustBenjaminiHochberg(p_het);
        p_hom_adj = adjustBenjaminiHochberg(p_hom);
    } else {
        p_het_adj = p_het;
        p_hom_adj = p_hom;
    }

    vector<pair<double, double>> relative_likelihoods = relativeLikelihoods(profile_likelihoods);

    cerr << "# pos\tprofile\t\tclass\tp_ho\t\tp_het\t\tl_ho\t\tl_het\t\trl_ho\t\trl_het" << endl;
    cout << scientific;
    cout.precision(2);
    for (const pair<int, Profile>&  pos_profile : positions) {
        const int& pos = pos_profile.first;
        const Profile& profile = pos_profile.second;

        int i = index_of[profile];
        cout << pos << '\t' << profile;
        if (relative_likelihoods[i].first < 1.0) {
            cout << '\t' << "het";
        } else if(relative_likelihoods[i].second < 1.0) {
            cout << '\t' << "hom";
        } else {
            cout << '\t' << "inc";
        }
//        if ((p_het_adj[i] <= args.p_value_threshold) && !(p_hom_adj[i] <= args.p_value_threshold)) {
//            cout << '\t' << "het";
//        } else if(!(p_het_adj[i] <= args.p_value_threshold) && (p_hom_adj[i] <= args.p_value_threshold)) {
//            cout << '\t' << "hom";
//        } else {
//            cout << '\t' << "inc";
//        }
        cout << '\t' << p_hom[i];
        cout << '\t' << p_het[i];
        cout << '\t' << profile_likelihoods[i].first << '\t' << profile_likelihoods[i].second;
        cout << '\t' << relative_likelihoods[i].first << '\t' << relative_likelihoods[i].second;
        cout << endl;
    }
}

void printHelp(char* program_name) {
    cout << "Usage: " << program_name << " [options] [input files]" << endl;
    cout << "Options:" << endl;
    vector<pair<string, string>> options = {
        {"-h", "print this help message and exit"},
        {"-p NUM", "p value threshold"},
        {"-c CORRECTION", "correction for multiple testing, one of 'bonf' (Bonferroni), 'bh' (Benjamini-Hochberg), 'none'"}
        // {"-o FILE", {"print output to FILE"}
    };
    for (const auto& option_desc : options) {
        cout << '\t' << option_desc.first << "\t\t\t" << option_desc.second << endl;
    }
}

int main(int argc, char** argv) {
    string optstring = "hp:c:";
    for (char opt = getopt(argc, argv, optstring.c_str()); opt != -1; opt = getopt(argc, argv, optstring.c_str())) {
        switch (opt) {
        case 'h':
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        case 'c':
            if (string(optarg) ==  "bonf") {
                args.correction = Correction::Bonferroni;
            } else if (string(optarg) == "bh") {
                args.correction = Correction::BH;
            } else if (string(optarg) == "none") {
                args.correction = Correction::None;
            } else{
                cerr << "Warning: unknown correction " << optarg << ". Using default." << endl;
            }
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
