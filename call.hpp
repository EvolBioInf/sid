#pragma once
// #ifndef CALL_H
// #define CALL_H

#include <array>
#include <iostream>
#include <vector>

#include "pileup_parser.hpp"

std::vector<PileupLine> readFile(std::istream& in, bool parse_base_qualities, bool parse_mapping_qualities);

typedef struct {
	std::array<uint16_t, 4> profile;
	uint32_t coverage;
	uint32_t count;
} UniqueProfile;

std::vector<UniqueProfile> countUniqueProfiles(const std::vector<PileupLine>& pileup);

typedef struct {
	std::string label {"none"};
} OutputRecord;

std::vector<OutputRecord> callBayes(std::istream& in);

// #endif
