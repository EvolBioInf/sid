#include <algorithm>
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
    string selection = "rel";
} args {};

double inline aic(double likelihood, double num_params) {
    return 2 * num_params - 2 * log(likelihood);
}

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

    GenomeParameters gp = estimateGenomeParameters(profiles, counts);

    string sep = ",";
    cout << scientific;
    cout.precision(4);
    if (args.selection == "ratio") {
        // compute p-values
        vector<double> p_hom;
        vector<double> p_het;
        for (int i = 0; i < profiles.size(); ++i) {
            // p value for heterozygous, H_0: homozygous more likely
            p_het.push_back(likelihoodRatioTest(gp.hom_likelihoods[i], gp.het_likelihoods[i]));
            // p value for homozygous, H_0: heterozygous more likely
            p_hom.push_back(likelihoodRatioTest(gp.het_likelihoods[i], gp.hom_likelihoods[i]));
        }

        // adjust p-values
        vector<double> p_het_adj;
        vector<double> p_hom_adj;
        if (args.correction == Correction::Bonferroni) {
            p_het_adj = adjustBonferroni(p_het, profiles.size());
            p_hom_adj = adjustBonferroni(p_hom, profiles.size());
        } else if (args.correction == Correction::BH) {
            p_het_adj = adjustBenjaminiHochberg(p_het);
            p_hom_adj = adjustBenjaminiHochberg(p_hom);
        } else {
            p_het_adj = p_het;
            p_hom_adj = p_hom;
        }

        cout << "#pos" + sep + "profile" + sep + "class" + sep + "p_hom" + sep + "p_het" << endl;
        for (const pair<int, Profile>&  pos_profile : positions) {
            const int& pos = pos_profile.first;
            const Profile& profile = pos_profile.second;

            int i = index_of[profile];
            cout << pos << sep << profile << sep;

            if ((p_het_adj[i] <= args.p_value_threshold) && !(p_hom_adj[i] <= args.p_value_threshold)) {
                cout << "het";
            } else if(!(p_het_adj[i] <= args.p_value_threshold) && (p_hom_adj[i] <= args.p_value_threshold)) {
                cout << "hom";
            } else {
                cout << "inc";
            }
            cout << sep << p_hom_adj[i] << sep << p_het_adj[i] << endl;
        }
    } else if (args.selection == "rel") {
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "reL_hom" + sep + "reL_het" << endl;
        for (const pair<int, Profile>&  pos_profile : positions) {
            const int& pos = pos_profile.first;
            const Profile& profile = pos_profile.second;

            int i = index_of[profile];

            double het_aic = aic(gp.het_likelihoods[i], 1);
            double hom_aic = aic(gp.hom_likelihoods[i], 1);
            double het_reL = 1;
            double hom_reL = 1;
            if (het_aic < hom_aic) {
                hom_reL = exp((het_aic - hom_aic) / 2.0);
            } else {
                het_reL = exp((hom_aic - het_aic) / 2.0);
            }

            cout << pos << sep << profile << sep;
            if (hom_reL < 1.0) {
                cout << "het";
            } else if(het_reL < 1.0) {
                cout << "hom";
            } else {
                cout << "inc";
            }
            cout << sep << hom_reL << sep << het_reL << endl;
        }
    } else if (args.selection == "map") {
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "ap_hom" + sep + "ap_het" << endl;
        for (const pair<int, Profile>&  pos_profile : positions) {
            const int& pos = pos_profile.first;
            const Profile& profile = pos_profile.second;

            int i = index_of[profile];

            double hom_ap = gp.hom_likelihoods[i] * (1 - gp.heterozygosity);
            double het_ap = gp.het_likelihoods[i] * gp.heterozygosity;

            cout << pos << sep << profile << sep;
            if (het_ap > hom_ap) {
                cout << "het";
            } else {
                cout << "hom";
            }
            cout << sep << hom_ap << sep << het_ap << endl;
        }
    } else if (args.selection == "bayes") {
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "p_hom" + sep + "p_het" << endl;
        for (const pair<int, Profile>&  pos_profile : positions) {
            const int& pos = pos_profile.first;
            const Profile& profile = pos_profile.second;

            int i = index_of[profile];

            double hom_ap = gp.hom_likelihoods[i] * (1 - gp.heterozygosity);
            double het_ap = gp.het_likelihoods[i] * gp.heterozygosity;

            double p_hom = hom_ap / (hom_ap + het_ap);
            double p_het = het_ap / (hom_ap + het_ap);

            cout << pos << sep << profile << sep;
            if (p_het > p_hom) {
                cout << "het";
            } else {
                cout << "hom";
            }
            cout << sep << p_hom << sep << p_het << endl;
        }
    } else if (args.selection == "local") {
        vector<double> errors;
        vector<string> genotypes;

        for (const Profile& p : profiles) {
            double max_likelihood = -1.0;
            string ml_genotype = "";
            double ml_error = -1.0;
            // homozygous
            // ml_error = n2+n3+n4 / n1+n2+n3+n4
            for (int i = 0; i < 4; ++i) {
                double error = (double)(p[COV] - p[i]) / p[COV];
                double l = profileLikelihoodHomozygous(p, error, i);
                if (l > max_likelihood) {
                    max_likelihood = l;
                    ml_genotype = to_string(i) + to_string(i);
                    ml_error = error;
                }
            }
            // heterozygous
            // ml_error = 1.5 * (n3+n4)/(n1+n2+n3+n4)
            for (int i = 0; i < 4; ++i) {
                for (int j = i+1; j < 4; ++j) {
                    double error = 1.5 * (double)(p[COV] - p[i] - p[j])/p[COV];
                    // error per base cannot exceed 1.0
                    error = min(error, 1.0);
                    double l = profileLikelihoodHeterozygous(p, error, i, j);
                    if (l > max_likelihood) {
                        max_likelihood = l;
                        ml_genotype = to_string(i) + to_string(j);
                        ml_error = error;
                    }
                }
            }
            errors.push_back(ml_error);
            genotypes.push_back(ml_genotype);
        }
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "gen" + sep + "err" << endl;
        for (const pair<int, Profile>& pos_profile : positions) {
            const int& pos = pos_profile.first;
            const Profile& profile = pos_profile.second;

            int i = index_of[profile];

            cout << pos << sep << profile << sep;
            if (genotypes[i][0] == genotypes[i][1]) {
                cout << "hom";
            } else {
                cout << "het";
            }
            cout << sep << genotypes[i] << sep << errors[i];
            cout << endl;
        }
    }
}

void printHelp(const char* program_name, const vector<struct option>& options, const vector<string>& descriptions) {
    cout << "Usage: " << program_name << " [options] [input files]" << endl;
    cout << "Options:" << endl;
    for (int i = 0; i < options.size() && i < descriptions.size(); ++i) {
        cout << '\t';
        cout << '-' << (char) options[i].val;
        cout << ", --" << options[i].name;
        cout << "\t" << descriptions[i] << endl;
    }
}

int main(int argc, char** argv) {
    vector<struct option> options;
    vector<string> descriptions;

    options.push_back({"help", no_argument, nullptr, 'h'});
    descriptions.push_back("Print this help message and exit");

    options.push_back({"correction", required_argument, nullptr, 'c'});
    descriptions.push_back("Correction for multiple testing, one of 'bonf' (Bonferroni), 'bh' (Benjamini-Hochberg), or 'none'");

    options.push_back({"significance", required_argument, nullptr, 'p'});
    descriptions.push_back("Significance level to use for likelihood ratio testing");

    options.push_back({"selection", required_argument, nullptr, 's'});
    descriptions.push_back("Model selection procedure, one of 'rel' (relative likelihood), 'ratio' (likelihood ratio test)");

    // end marker for getopt_long
    options.push_back({0,0,0,0});

    string optstring;
    for (struct option opt : options) {
        optstring += opt.val;
        // use the fact that no_argument = 0, required_argument = 1, optional_argument = 2
        for (int i = 0; i < opt.has_arg; ++i) {
            optstring += ':';
        }
    }

    int optindex = -1;
    char opt = 0;

    double p_value = -1.0;
    while ((opt = getopt_long(argc, argv, optstring.c_str(), options.data(), &optindex)) != -1) {
        string value = optarg == nullptr ? "" : string(optarg);
        switch (opt) {
        case 'h':
        printHelp(argv[0], options, descriptions);
            exit(EXIT_SUCCESS);
            break;
        case 'c':
            if (value ==  "bonf") {
                args.correction = Correction::Bonferroni;
            } else if (value == "bh") {
                args.correction = Correction::BH;
            } else if (value == "none") {
                args.correction = Correction::None;
            } else{
                cerr << "Warning: unknown correction " << optarg << ". Using default." << endl;
            }
            break;
        case 'p':
            p_value = atof(optarg);
            if (p_value <= 0.0) {
                cerr << "Invalid p value threshold or parse error: " << optarg << endl;
                exit(EXIT_FAILURE);
            } else {
                args.p_value_threshold = p_value;
            }
            break;
        case 's':
            if (value != "rel" && value != "ratio" && value != "bayes" && value != "map" && value != "local") {
                cerr << "Unknown model selection procedure: " << value << endl;
                exit(EXIT_FAILURE);
            }
            args.selection = value;
            break;
        default:
            cerr << "Unknown option character: " << opt << " (" << (int)opt << ")" << endl;
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
