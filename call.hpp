#pragma once
// #ifndef CALL_H
// #define CALL_H

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "pileup.hpp"

std::vector<PileupLine> readFile(std::istream& in, bool parse_base_qualities, bool parse_mapping_qualities);

typedef struct {
	std::string label {"none"};
	std::string genotype {"NN"};
	double confidence_homozygous;
	double confidence_heterozygous;
	std::string confidence_type {"unspecified"};
	std::vector<std::string> additional_data;
} Classification;

typedef struct {
	std::string chromosome_name;
	int position;
	Classification classification;
} OutputRecord;

inline std::ostream& operator<<(std::ostream& os, const OutputRecord& r) {
	os << r.chromosome_name;
	os << ',' << r.position;
	os << ',' << r.classification.label;
	os << ',' << r.classification.genotype;
	os << ',' << r.classification.confidence_homozygous;
	os << ',' << r.classification.confidence_heterozygous;
	os << ',' << r.classification.confidence_type;
    return os;
}

std::vector<OutputRecord> callLikelihoodRatio(std::istream& in);
std::vector<OutputRecord> callBayes(std::istream& in);
std::vector<OutputRecord> callSiteMLError(std::istream& in, const bool estimate_prior, double prior, double error_threshold);
std::vector<OutputRecord> callQualityBasedSimple(std::istream& in, const bool estimate_prior, double prior);


// #endif
