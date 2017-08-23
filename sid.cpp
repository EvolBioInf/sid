#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_set>
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
} args {};

typedef struct {
    int pos {-1};
    char reference_base {'N'};
    char* read;
} PileupLine;

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

    result.read = strtok(NULL, delim);
    if (result.read == NULL) {
        result.pos = -1;
    }
    return result;
}

typedef struct {
    int num_sites {0};
    vector<int> positions {};
    vector<Profile> profiles {};
} PileupData;

PileupData readlinePileup(FILE* input) {
    PileupData result;

    int pos;
    char reference_base;
    int count = 0;

    size_t read_buffer_size = 10000;
    char* read = (char*)malloc(read_buffer_size*sizeof(char));

    char* line = (char*)malloc(10000*sizeof(char));
    size_t line_length = 0;
    ssize_t num_read = -1;
    while((num_read = getline(&line, &line_length, input)) != -1) {
        if (line_length > read_buffer_size) {
            read = (char*)realloc((void*)read, line_length*sizeof(char));
            read_buffer_size = line_length;
        }

        PileupLine plpline = parsePileupLine(line);
        if (plpline.pos < 0) {
            cerr << "# Read error: input must be valid 'samtools mpileup' output." << endl;
            exit(EXIT_FAILURE);
        }
        Profile p = parseRead(plpline.read, plpline.reference_base);

        if (p[COV] >= 4) {
            result.positions.push_back(plpline.pos);
            result.profiles.push_back(p);
            result.num_sites += 1;
        }
        ++count;
    }
    free(line);
    free(read);

    if(result.num_sites == 0) {
        cerr << "# No profiles found with required minimum coverage of 4!" << endl;
        exit(EXIT_FAILURE);
    }
    cerr << "# " << result.num_sites << " of " << count << " sites with required coverage > 4" << endl;
    return result;
}

void callVariants(const PileupData& pileup, const vector<Profile>& profiles, const GenomeParameters& gp) {
    auto t1 = chrono::high_resolution_clock::now();

    map<Profile, int> index_of {};
    for (int i = 0; i < profiles.size(); ++i) {
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
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

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
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

            int i = index_of[profile];

            // compute Akaike Information Criterion with 1 degree of freedom
            auto aic = [] (double likelihood) {
                return 2 * 1 - 2 * log(likelihood);
            };
            double het_aic = aic(gp.het_likelihoods[i]);
            double hom_aic = aic(gp.hom_likelihoods[i]);

            // compute relative likelihoods
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
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

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
        vector<double> p_het (profiles.size());
        vector<double> p_hom (profiles.size());
        vector<string> output (profiles.size());
        for (int i = 0; i < profiles.size(); ++i) {
            double hom_ap = gp.hom_likelihoods[i] * (1 - gp.heterozygosity);
            double het_ap = gp.het_likelihoods[i] * gp.heterozygosity;
            p_hom[i] = hom_ap / (hom_ap + het_ap);
            p_het[i] = het_ap / (hom_ap + het_ap);

            string label = "hom";
            if (p_het[i] > p_hom[i]) {
                label = "het";
            }
            char buffer[256];
            sprintf(buffer, ",%d.%d.%d.%d,%s,%e,%e", profiles[i][0], profiles[i][1], profiles[i][2],profiles[i][3], label.c_str(), p_hom[i], p_het[i]);
            output[i] = string(buffer);
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
    } else if (args.selection == "localmax") {
        const double ERROR_THRESHOLD = 0.1;
        vector<double> errors (profiles.size());
        vector<string> genotypes (profiles.size());

        int i = 0;
        for (const Profile& p : profiles) {
            double max_likelihood = -1.0;
            string ml_genotype = "";
            double ml_error = -1.0;
            // homozygous
            // ml_error = n2+n3+n4 / n1+n2+n3+n4
            for (int i = 0; i < 4; ++i) {
                double error = (double)(p[COV] - p[i]) / p[COV];
                error = min(ERROR_THRESHOLD, error);
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
                    error = min(ERROR_THRESHOLD, error);
                    double l = profileLikelihoodHeterozygous(p, error, i, j);
                    if (l > max_likelihood) {
                        max_likelihood = l;
                        ml_genotype = to_string(i) + to_string(j);
                        ml_error = error;
                    }
                }
            }
            errors[i] = ml_error;
            genotypes[i] = ml_genotype;
            ++i;
        }
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "gen" + sep + "err" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

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
    } else if (args.selection == "local") {
        const double ERROR_THRESHOLD = 0.1;
        vector<double> errors (profiles.size());
        vector<string> genotypes (profiles.size());
        vector<double> pvalues (profiles.size());

        int ii = 0;
        for (const Profile& p : profiles) {
            array<int, 4> sorted_bases;
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
            error1 = min(ERROR_THRESHOLD, error1);
            double l1 = profileLikelihoodHomozygous(p, error1, largest_i);

            // heterozygous
            // ml_error = 1.5 * (n3+n4)/(n1+n2+n3+n4)
            double error2 = 1.5 * (double)(p[COV] - p[largest_i] - p[snd_largest_i])/p[COV];
            error2 = min(ERROR_THRESHOLD, error2);
            double l2 = profileLikelihoodHeterozygous(p, error2, largest_i, snd_largest_i);

            double p1 = likelihoodRatioTest(l2, l1);
            double p2 = likelihoodRatioTest(l1, l2);

            if (p1 < args.p_value_threshold && p2 > args.p_value_threshold) {
                errors[ii] = error1;
                genotypes[ii] = to_string(largest_i) + to_string(largest_i);
                pvalues[ii] = likelihoodRatioTest(l2, l1);
            } else if (p1 > args.p_value_threshold && p2 < args.p_value_threshold) {
                errors[ii] = error2;
                genotypes[ii] = to_string(largest_i) + to_string(snd_largest_i);
                pvalues[ii] = likelihoodRatioTest(l1, l2);
            } else {
                if (l1 > l2) {
                    errors[ii] = error1;
                    genotypes[ii] = to_string(largest_i) + to_string(largest_i);
                    pvalues[ii] = p1;
                } else {
                    errors[ii] = error2;
                    genotypes[ii] = to_string(largest_i) + to_string(snd_largest_i);
                    pvalues[ii] = p2;
                }
            }
            ++ii;
        }
        cout << "#pos" + sep + "profile" + sep + "class" + sep + "gen" + sep + "err" + sep + "p" << endl;
        for (int ii = 0; ii < pileup.num_sites; ++ii) {
            const int& pos = pileup.positions[ii];
            const Profile& profile = pileup.profiles[ii];

            int i = index_of[profile];

            cout << pos << sep << profile << sep;
            if (pvalues[i] > args.p_value_threshold) {
                cout << "inc";
            } else if (genotypes[i][0] == genotypes[i][1]) {
                cout << "hom";
            } else {
                cout << "het";
            }
            cout << sep << genotypes[i] << sep << errors[i] << sep << pvalues[i];
            cout << endl;
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
    PileupData pileup = readlinePileup(input);

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
