#include "catch.hpp"
#include "call.hpp"
#include "pileup_parser.hpp"

bool operator==(const UniqueProfile& a, const UniqueProfile& b) {
    return a.count == b.count && a.profile == b.profile;
}

std::ostream& operator<<(std::ostream& os, const UniqueProfile& p) {
    for (auto x : p.profile) {
        os << x << ",";
    }
    return os << p.count;
}

TEST_CASE("Identification of  unique profiles", "[countUniqueProfiles]") {
	SECTION("general functionality") {
        std::vector<PileupLine> lines {
            {"name", 0, 'N', {1,1,1,1}, {}, {}, {}, {}},
            {"name", 0, 'N', {2,2,2,2}, {}, {}, {}, {}},
            {"name", 0, 'N', {1,1,1,1}, {}, {}, {}, {}}
        };
        std::vector<UniqueProfile> unique_profiles {
            {{1,1,1,1}, 2},
            {{2,2,2,2}, 1}
        };
		REQUIRE(countUniqueProfiles(lines) == unique_profiles);
	}
	SECTION("empty list") {
        std::vector<PileupLine> nolines {};
        std::vector<UniqueProfile> empty {};

		REQUIRE(countUniqueProfiles(nolines) == empty);
	}
}