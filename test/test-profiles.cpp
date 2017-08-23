#include <cmath>
#include <limits>
#include <vector>

#include "profiles.hpp"

#include "catch.hpp"

using namespace std;

TEST_CASE("sequence reads in 'samtools mpileup' format are parsed correctly", "[parseRead]") {
    SECTION ("simple reads") {
        Profile aa {2,0,0,0,2};
        REQUIRE(parseRead("aA") == aa);

        Profile cc {0,2,0,0,2};
        REQUIRE(parseRead("cC") == cc);

        Profile gg {0,0,2,0,2};
        REQUIRE(parseRead("gG") == gg);

        Profile tt {0,0,0,2,2};
        REQUIRE(parseRead("tT") == tt);
    }
    SECTION ("empty read") {
        Profile p {0,0,0,0,0};
        REQUIRE(parseRead("") == p);
    }
    SECTION ("ignore read end") {
        Profile p {1,0,0,0,1};
        REQUIRE(parseRead("a$") == p);
    }
    SECTION("skip quality markers") {
        Profile p {1,0,0,0,1};
        REQUIRE(parseRead("a^a") == p);
        REQUIRE(parseRead("^aa") == p);
    }
    SECTION("skip indels") {
        Profile p {1,0,0,0,1};
        REQUIRE(parseRead("a+3act") == p);
        REQUIRE(parseRead("+3acta") == p);

        REQUIRE(parseRead("a-3act") == p);
        REQUIRE(parseRead("-3acta") == p);
    }
    SECTION("correctly handel reference bases") {
        Profile p {1,0,1,0,2};
        CHECK(parseRead("a.", 'g') == p);
        CHECK(parseRead(",g", 'a') == p);
        CHECK(parseRead("ag", 't') == p);
        CHECK(parseRead("ag", 'n') == p);
        CHECK(parseRead("ag") == p);
    }
}
