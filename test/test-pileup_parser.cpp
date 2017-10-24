#include "catch.hpp"
#include "pileup_parser.hpp"

#include <cstring>
#include <string>
#include <vector>

TEST_CASE("Quality parsing", "[parseQualities]") {
	SECTION("general functionality") {
		const char* qualities = "+5D";
		std::vector<uint8_t> result {10, 20, 35};
		int coverage = strlen(qualities);
		REQUIRE(parseQualities(qualities, coverage) == result);
	}
	SECTION("empty string") {
		const char* qualities = "";
		std::vector<uint8_t> result {};
		int coverage = 0;
		REQUIRE(parseQualities(qualities, coverage) == result);
	}
}

TEST_CASE("Base parsing", "[parseReadBases]") {
	SECTION("general functionality") {
		const char* bases = "AgACgt";
		int coverage = strlen(bases);
		auto result_bases = std::vector<char> {'A', 'G', 'A', 'C', 'G', 'T'};
		auto result_strands = std::vector<bool> {1, 0, 1, 1, 0, 0};
		auto result_counts = std::array<int, 4> {2, 1, 2, 1};
		auto result = parseReadBases(bases, 'N', coverage);
		REQUIRE(result.bases == result_bases);
		REQUIRE(result.strands == result_strands);
		REQUIRE(result.counts == result_counts);
	}
}

TEST_CASE("Line parsing", "[parsePileupLine]") {
	SECTION("general functionality") {
		char line[] = "chr19\t1337\tA\t6\tAgACgt\t++5D5\tDD55D";

		std::array<int, 4> result_counts {2, 1, 2, 1};
		std::vector<char> result_bases {'A', 'G', 'A', 'C', 'G', 'T'};
		std::vector<bool> result_strands  {1, 0, 1, 1, 0, 0};
		std::vector<uint8_t> result_base_qualities {10, 10, 20, 35, 20};
		std::vector<uint8_t> result_mapping_qualities {35, 35, 20, 20, 35};

		auto result = parsePileupLine(line, true, true);
		REQUIRE(result.chromosome_name == "chr19");
		REQUIRE(result.position == int(1337));
		REQUIRE(result.base_counts == result_counts);
		REQUIRE(result.bases == result_bases);
		REQUIRE(result.strands == result_strands);
		REQUIRE(result.base_qualities == result_base_qualities);
		REQUIRE(result.mapping_qualities == result_mapping_qualities);
	}
}