#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "profiles.hpp"
#include "likelihoods.hpp"
#include "stats.hpp"

using namespace std;

enum class Correction {Bonferroni, BH, None};

struct arguments {
    Correction correction {Correction::Bonferroni};
    double p_value_threshold {0.05};
    string selection = "rel";
    double error_threshold = 0.1;
} args {};

void callVariants(const PileupData& pileup, const vector<Profile>& profiles, const GenomeParameters& gp) {
    auto t1 = chrono::high_resolution_clock::now();

    map<Profile, unsigned int> index_of {};
    for (unsigned int i = 0; i < profiles.size(); ++i) {
        index_of.emplace(profiles[i], i);
    }
    auto t2 = chrono::high_resolution_clock::now();
    auto tind = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cerr << "# Building of reverse index took " << tind.count() << " seconds" << endl;

    //overwritten later, only here for scope and auto
    auto t3 = chrono::high_resolution_clock::now();

    string sep = ",";
    cout << scientific;
    cout.precision(4);
    if (args.selection == "ratio") {
        // compute p-values
        vector<double> p_hom;
        p_hom.reserve(profiles.size());
        vector<double> p_het;
        p_het.reserve(profiles.size());
        for (unsigned int i = 0; i < profiles.size(); ++i) {
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

        vector<string> output;
        output.reserve(profiles.size());
        for (size_t i = 0; i < profiles.size(); ++i) {
            const Profile& p = profiles[i];
            string label ("inc");
            if ((p_het_adj[i] <= args.p_value_threshold) && !(p_hom_adj[i] <= args.p_value_threshold)) {
                label = "het";
            } else if(!(p_het_adj[i] <= args.p_value_threshold) && (p_hom_adj[i] <= args.p_value_threshold)) {
                label = "hom";
            }
            char buffer[256];
            sprintf(buffer, "%d.%d.%d.%d,%s,%e,%e", p[A], p[C], p[G], p[T], label.c_str(), p_hom_adj[i], p_het_adj[i]);
            output.emplace_back(buffer);
        }

        cout << "#pos" + sep + "profile" + sep + "class" + sep + "p_hom" + sep + "p_het" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

            int i = index_of[profile];
            printf("%d,%s\n", pos, output[i].c_str());
        }
    } else if (args.selection == "rel") {
        vector<string> output;
        output.reserve(profiles.size());
        for (size_t i = 0; i < profiles.size(); ++i) {
            const Profile& p = profiles[i];
            // compute Akaike Information Criterion with 1 degree of freedom
            auto aic = [] (double likelihood) {
                return 2 * 1 - 2 * log(likelihood);
            };
            double het_aic = aic((double)gp.het_likelihoods[i]);
            double hom_aic = aic((double)gp.hom_likelihoods[i]);

            // compute relative likelihoods
            double het_reL = 1;
            double hom_reL = 1;
            if (het_aic < hom_aic) {
                hom_reL = exp((het_aic - hom_aic) / 2.0);
            } else {
                het_reL = exp((hom_aic - het_aic) / 2.0);
            }
            string label ("inc");
            if (hom_reL < 1.0) {
                label = "het";
            } else if(het_reL < 1.0) {
                label = "hom";
            }
            char buffer [256];
            sprintf(buffer, "%d.%d.%d.%d,%s,%e,%e", p[A], p[C], p[G], p[T], label.c_str(), hom_reL, het_reL);
            output.emplace_back(buffer);

        }
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "reL_hom" + sep + "reL_het" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& p = pileup.profiles[ii];

            int i = index_of[p];

            printf("%d,%s\n", pos, output[i].c_str());
        }
    } else if (args.selection == "map") {
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "ap_hom" + sep + "ap_het" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

            int i = index_of[profile];

            auto hom_ap = gp.hom_likelihoods[i] * (1 - gp.heterozygosity);
            auto het_ap = gp.het_likelihoods[i] * gp.heterozygosity;

            cout << pos << sep << profile << sep;
            if (het_ap > hom_ap) {
                cout << "het";
            } else {
                cout << "hom";
            }
            cout << sep << hom_ap << sep << het_ap << endl;
        }
    } else if (args.selection == "bayes") {
        vector<string> output;
        output.reserve(profiles.size());
        for (int i = 0; i < profiles.size(); ++i) {
            auto hom_ap = gp.hom_likelihoods[i] * (1 - gp.heterozygosity);
            auto het_ap = gp.het_likelihoods[i] * gp.heterozygosity;
            auto p_hom = hom_ap / (hom_ap + het_ap);
            auto p_het = het_ap / (hom_ap + het_ap);

            string label = "hom";
            if (p_het > p_hom) {
                label = "het";
            }
            char buffer[256];
            sprintf(buffer, ",%d.%d.%d.%d,%s,%Le,%Le", profiles[i][0], profiles[i][1], profiles[i][2],profiles[i][3], label.c_str(), p_hom, p_het);
            output.emplace_back(buffer);
        }

        t3 = chrono::high_resolution_clock::now();
        auto tcall = chrono::duration_cast<chrono::duration<double>>(t3 - t2);
        cerr << "# SNP calling took " << tcall.count() << " seconds" << endl;

        cout << "#pos" + sep + "profile" + sep + "class" + sep + "p_hom" + sep + "p_het" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];
            int i = index_of[profile];

            printf("%d%s\n", pos, output[i].c_str());
        }
    } else if (args.selection == "local") {
        vector<string> output;
        output.reserve(profiles.size());
        for (const Profile& p : profiles) {
            int largest = -1;
            int snd_largest = -1;
            int largest_i = -1;
            int snd_largest_i = -1;
            for (int i = 0; i < 4; ++i) {
                if (p[i] > largest) {
                    snd_largest = largest;
                    snd_largest_i = largest_i;
                    largest = p[i];
                    largest_i = i;
                } else if (p[i] > snd_largest) {
                    snd_largest = p[i];
                    snd_largest_i = i;
                }
            }

            // homozygous
            // ml_error = n2+n3+n4 / n1+n2+n3+n4
            double error1 = (double)(p[COV] - p[largest_i]) / p[COV];
            error1 = min(args.error_threshold, error1);
            long double l1 = profileLikelihoodHomozygous(p, error1, largest_i);

            // heterozygous
            // ml_error = 1.5 * (n3+n4)/(n1+n2+n3+n4)
            double error2 = 1.5 * (double)(p[COV] - p[largest_i] - p[snd_largest_i])/p[COV];
            error2 = min(args.error_threshold, error2);
            long double l2 = profileLikelihoodHeterozygous(p, error2, largest_i, snd_largest_i);

            double p1 = likelihoodRatioTest(l2, l1);
            double p2 = likelihoodRatioTest(l1, l2);

            double error = -1.0;
            double pvalue = -1.0;
            string genotype;
            if (p1 < args.p_value_threshold && p2 > args.p_value_threshold) {
                error = error1;
                genotype = to_string(largest_i) + to_string(largest_i);
                pvalue = p1;
            } else if (p1 > args.p_value_threshold && p2 < args.p_value_threshold) {
                error = error2;
                genotype = to_string(largest_i) + to_string(snd_largest_i);
                pvalue = p2;
            } else {
                if (l1 > l2) {
                    error = error1;
                    genotype = to_string(largest_i) + to_string(largest_i);
                    pvalue = p1;
                } else {
                    error = error2;
                    genotype = to_string(largest_i) + to_string(snd_largest_i);
                    pvalue = p2;
                }
            }

            string label {"hom"};
            if (pvalue > args.p_value_threshold) {
                label = "inc";
            } else if (genotype[0] != genotype[1]) {
                label = "het";
            }
            char buffer[256];
            sprintf(buffer, "%d.%d.%d.%d,%s,%s,%e,%e", p[A], p[C], p[G], p[T], label.c_str(), genotype.c_str(), error, pvalue);
            output.emplace_back(buffer);
        }
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "gen" + sep + "err" + sep + "p" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

            int i = index_of[profile];
            printf("%d,%s\n", pos, output[i].c_str());
        }
    } else {
        cerr << "# Unknown model selection procedure: " << args.selection << endl;
        exit(EXIT_FAILURE);
    }
    auto t4 = chrono::high_resolution_clock::now();
    auto tout = chrono::duration_cast<chrono::duration<double>>(t4 - t3);
    cerr << "# Output took " << tout.count() << " seconds" << endl;
}


void processFile(FILE* input) {
    auto t1 = chrono::high_resolution_clock::now();
    PileupData pileup = readPileup(input);

    auto t2 = chrono::high_resolution_clock::now();
    vector<Profile> unique_profiles (pileup.profiles);

    sort(unique_profiles.begin(), unique_profiles.end());
    vector<int> counts {1};
    // count number of occurences of each profile
    for(auto it = unique_profiles.begin() + 1; it != unique_profiles.end(); ++it) {
        if (*it != *(it - 1)) {
            counts.push_back(1);
        } else {
            ++counts.back();
        }
    }
    // remove duplicate profiles
    auto new_end = unique(unique_profiles.begin(), unique_profiles.end());
    unique_profiles.resize(distance(unique_profiles.begin(), new_end));
    unique_profiles.shrink_to_fit();
    auto t3 = chrono::high_resolution_clock::now();

    cerr << "# " << counts.size() << " distinct profiles" << endl;

    GenomeParameters gp = estimateGenomeParameters(unique_profiles, counts);

    auto t4 = chrono::high_resolution_clock::now();

    auto tread = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    auto tsort = chrono::duration_cast<chrono::duration<double>>(t3 - t2);
    auto testimate = chrono::duration_cast<chrono::duration<double>>(t4 - t3);
    cerr << fixed;
    cerr << "# Reading data took " << tread.count() << " seconds" << endl;
    cerr << "# Sorting data took " << tsort.count() << " seconds" << endl;
    cerr << "# Esimating parameters took " << testimate.count() << " seconds" << endl;

    callVariants(pileup, unique_profiles, gp);
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

    options.push_back({"error_threshold", required_argument, nullptr, 'e'});
    descriptions.push_back("Largest allowed error per site for local selection procedure");

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
    double error_threshold = 0.1;
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
                cerr << "# Warning: unknown correction " << optarg << ". Using default." << endl;
            }
            break;
        case 'p':
            p_value = atof(optarg);
            if (p_value <= 0.0) {
                cerr << "# Invalid p value threshold or parse error: " << optarg << endl;
                exit(EXIT_FAILURE);
            } else {
                args.p_value_threshold = p_value;
            }
            break;
        case 's':
            /* if (value != "rel" && value != "ratio" && value != "bayes" && value != "map" && value != "local") { */
            /*     cerr << "Unknown model selection procedure: " << value << endl; */
            /*     exit(EXIT_FAILURE); */
            /* } */
            args.selection = value;
            break;
        case 'e':
            error_threshold = atof(optarg);
            args.error_threshold = error_threshold;
            break;
        default:
            cerr << "# Unknown option character: " << opt << " (" << (int)opt << ")" << endl;
        }
    }
    if (optind < argc) {
        FILE* f = fopen(argv[optind], "r");
        if (f != nullptr) {
            processFile(f);
        } else {
            cerr << "# Error opening file: " << argv[optind] << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        processFile(stdin);
    }
}
