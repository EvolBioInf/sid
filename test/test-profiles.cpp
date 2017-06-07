#include <cmath>
#include <limits>
#include <vector>

#include "profiles.hpp"

#include "catch.hpp"

using namespace std;

TEST_CASE("arithmetic operatons on profiles work correctly", "[operators]") {
    SECTION("addition") {
        Profile a {0,1,0,1,2};
        Profile a_ {a};
        Profile b {1,1,0,0,2};
        Profile b_ {b};
        Profile result {1,2,0,1,4};

        REQUIRE(a + b == result);
        REQUIRE(a == a_);
        REQUIRE(b == b_);
    }
    SECTION("multiplication") {
        Profile a {0,1,0,1,2};
        Profile two_a {0,2,0,2,4};

        REQUIRE(a * 2 == two_a);
    }
}

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
    //SECTION("correctly handel reference bases") {
    //    Profile p {1,1,0,0,2};
    //    REQUIRE(parseRead("a.", 'g') == p);
    //    REQUIRE(parseRead(",g", 'g') == p);
    //    REQUIRE(parseRead("ag", 't') == p);
    //    REQUIRE(parseRead("ag", 'n') == p);
    //    REQUIRE(parseRead("ag") == p);
    //}
}
