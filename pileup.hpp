#pragma once

#include <array>
#include <numeric>
#include <vector>

using profile_t = std::array<uint16_t, 4>;

typedef struct {
	std::string chromosome_name;
	int position {-1};
	char reference_base {'N'};
	profile_t base_counts;
	std::vector<char> bases;
	std::vector<bool> strands;
	std::vector<uint8_t> base_qualities;
	std::vector<uint8_t> mapping_qualities;
} PileupLine;

PileupLine parsePileupLine(char* line, bool parse_base_qualities, bool parse_mapping_qualities);

typedef struct {
	std::vector<char> bases;
	std::vector<bool> strands; // 1 ~ forward strand, 0 ~ reverse strand
	profile_t counts;
} ReadStack;

ReadStack parseReadBases(const char* read_bases, char reference, int coverage);

std::vector<uint8_t> parseQualities(const char* base_qualities, int coverage);

typedef struct UniqueProfile {
    profile_t profile;
	uint32_t count;
    uint32_t coverage;
	UniqueProfile () : profile{0,0,0,0}, count{0}, coverage{0} {}
	UniqueProfile (profile_t p, uint32_t count) : profile{p}, count{count} {
		coverage = std::accumulate(profile.begin(), profile.end(), 0);
	}
} UniqueProfile;

std::vector<UniqueProfile> countUniqueProfiles(const std::vector<PileupLine>&);

std::array<double, 4> computeNucleotideDistribution(const std::vector<UniqueProfile>&);