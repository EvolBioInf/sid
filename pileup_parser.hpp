#pragma once
// #ifndef PLPPARSER_H
// #define PLPPARSER_H

#include <array>
#include <vector>

typedef struct {
	std::string chromosome_name;
	int position {-1};
	char reference_base {'N'};
	std::array<int, 4> base_counts;
	std::vector<char> bases;
	std::vector<bool> strands;
	std::vector<uint8_t> base_qualities;
	std::vector<uint8_t> mapping_qualities;
} PileupLine;

PileupLine parsePileupLine(char* line, bool parse_base_qualities, bool parse_mapping_qualities);

typedef struct {
	std::vector<char> bases;
	std::vector<bool> strands; // 1 ~ forward strand, 0 ~ reverse strand
	std::array<int, 4> counts;
} ReadStack;

ReadStack parseReadBases(const char* read_bases, char reference, int coverage);

std::vector<uint8_t> parseQualities(const char* base_qualities, int coverage);

// #endif
