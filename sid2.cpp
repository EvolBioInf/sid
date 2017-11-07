#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>

#include <getopt.h>

#include "call.hpp"

typedef struct {
    std::string method {"local"};
    bool estimate_prior = false;
    double snp_prior = -1;
    double significance_level = 0.05;
    double site_error_threshold = 0.1;
} GlobalOptions;

typedef struct {
    const char* name;
    const int has_arg;
    const char* description;
    const std::function<void(GlobalOptions&, const char*)> update;
} OptionDescription;

const static std::map<char, OptionDescription> OPTIONS {
    {'h', {"help", no_argument, "Print this help message", 
        [](GlobalOptions& o __attribute__((unused)), const char* optarg __attribute__((unused))) {
            std::cout << "sid [flags] input_file" << '\n';
            for (auto o : OPTIONS) {
                std::cout << '\t' << '-' << o.first;
                if (o.second.has_arg > 0) {
                    std::cout << ' ' << o.second.name;
                }
                std::cout << "\t" << o.second.description << '\n';
            } 
        }}},
    {'m', {"METHOD", required_argument, "Select the method to use for SNP calling: 'likelihood_ratio' , 'bayes', 'local' or 'quality', default: local",
        [](GlobalOptions& o, const char* optarg) { 
            o.method = optarg;
        }}},
    {'r', {"PRIOR", required_argument, "Use the given prior for SNPs, applicable for methods 'local', 'quality'. Conflicts -R. Default: no prior",
        [](GlobalOptions& o, const char* optarg) {
            o.snp_prior = atof(optarg);
        }}},
    {'R', {"", no_argument, "Estimate SNP prior from data, applicable for methods 'likelihood_ratio', 'local', 'quality'. Conflicts -r.",
        [](GlobalOptions& o, const char* optarg __attribute__((unused))) {
            o.estimate_prior = true;
        }}},
    {'p', {"LEVEL", required_argument, "Significance level for statistical tests, only applicable for methods 'likelihood_ratio', 'local'. Default: 0.05",
        [](GlobalOptions& o, const char* optarg) {
            o.significance_level = atof(optarg);
        }}},
    {'E', {"ERROR", required_argument, "Maximum allowed site error rate for 'local' method. Default: 0.1",
        [](GlobalOptions& o, const char* optarg) {
            o.site_error_threshold = atof(optarg);
        }}}
};

std::string buildOptString(const decltype(OPTIONS)& options) {
    std::string optstring;
    for (auto o : options) {
        optstring += o.first;
        for (int i = 0; i < o.second.has_arg; ++i) {
            optstring += ':';
        }
    }
    return optstring;
}

int main(int argc, char** argv) {
    GlobalOptions global_options;
    auto optstring = buildOptString(OPTIONS);
    int flag = 0;
    while((flag = getopt(argc, argv, optstring.c_str())) != -1) {
        auto option = OPTIONS.find(static_cast<char>(flag));
        if (option != OPTIONS.end()) {
            option->second.update(global_options, optarg);
        } else {
            exit(EXIT_FAILURE);
        }
    }
    if (optind < argc) {
        const char* input_file_path = argv[optind];
        auto in = std::ifstream(input_file_path);
        if (!in) {
            std::cerr << "Could not open file: " << input_file_path << std::endl;
            exit(EXIT_FAILURE);
        }

        std::vector<OutputRecord> output_records;
        if (global_options.method == "local") {
            output_records = callSiteMLError(in, global_options.estimate_prior, global_options.snp_prior, global_options.site_error_threshold);
        } else if (global_options.method == "bayes") {
            output_records = callBayes(in);
        } else if (global_options.method == "likelihood_ratio") {
            output_records = callLikelihoodRatio(in, global_options.estimate_prior);
        } else if (global_options.method == "quality") {
            output_records = callQualityBasedSimple(in, global_options.estimate_prior, global_options.snp_prior);
        }

        std::cout << "chrom,pos,label,gt,hom_conf,het_conf,conf_type" << std::endl;
        for (const auto& record : output_records) {
            std::cout << record << '\n';
        }
    } else {
        std::cerr << "No file name given!" << std::endl;
        exit(EXIT_FAILURE);
    }
}
