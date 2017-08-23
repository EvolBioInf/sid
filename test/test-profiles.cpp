#include <cmath>
#include <limits>
#include <vector>

#include "profiles.hpp"

#include "catch.hpp"

using namespace std;

Profile makeProfile(uint16_t a, uint16_t c, uint16_t g, uint16_t t) {
    return Profile {a,c,g,t,(uint16_t)(a+c+g+t)};
}

TEST_CASE("read base column 'samtools mpileup' format is parsed correctly", "[parseReadBases]") {
    SECTION ("simple reads") {
        REQUIRE(parseReadBases("aA", 'n') == makeProfile(2,0,0,0));

        REQUIRE(parseReadBases("cC", 'n') == makeProfile(0,2,0,0));

        REQUIRE(parseReadBases("gG", 'n') == makeProfile(0,0,2,0));

        REQUIRE(parseReadBases("tT", 'n') == makeProfile(0,0,0,2));
    }
    SECTION ("empty read") {
        REQUIRE(parseReadBases("", 'n') == makeProfile(0,0,0,0));
    }
    SECTION ("ignore read end") {
        REQUIRE(parseReadBases("a$", 'n') == makeProfile(1,0,0,0));
    }
    SECTION("skip quality markers") {
        Profile p {1,0,0,0,1};
        REQUIRE(parseReadBases("a^a", 'n') == p);
        REQUIRE(parseReadBases("^aa", 'n') == p);
    }
    SECTION("skip indels") {
        Profile p {1,0,0,0,1};
        REQUIRE(parseReadBases("a+3act", 'n') == p);
        REQUIRE(parseReadBases("+3acta", 'n') == p);

        REQUIRE(parseReadBases("a-3act", 'n') == p);
        REQUIRE(parseReadBases("-3acta", 'n') == p);
    }
    SECTION("correctly handel reference bases") {
        Profile p {1,0,1,0,2};
        CHECK(parseReadBases("a.", 'g') == p);
        CHECK(parseReadBases(",g", 'a') == p);
        CHECK(parseReadBases("ag", 't') == p);
        CHECK(parseReadBases("ag", 'n') == p);
        CHECK(parseReadBases("ag", 'n') == p);
    }
}
