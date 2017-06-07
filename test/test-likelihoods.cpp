#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include<sstream>

#include "likelihoods.hpp"

#include "catch.hpp"

using namespace std;

// enable test suite to print double arrays
namespace std {
    std::ostream& operator<<(std::ostream& os, const array<double, 4>& a) {
        return os << "[" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << "]";
    }
}

TEST_CASE("binomial probabilities are computed correctly", "[binom_probability_gamma]") {
    // precomputed with python scipy.stats.binom.pmf(k, n, p)
    SECTION("small numbers") {
        vector<int> ns {19, 26, 22, 16, 11, 29, 19, 16, 8, 20, 18, 5, 27, 19, 7, 5, 7, 17, 14, 20};
        vector<int> ks {3, 19, 14, 3, 1, 9, 15, 5, 2, 5, 11, 2, 23, 17, 5, 2, 2, 15, 7, 18};
        vector<double> ps {0.465, 0.79, 0.408, 0.246, 0.198, 0.78, 0.141, 0.458, 0.498, 0.726,
            0.371, 0.677, 0.944, 0.856, 0.391, 0.415, 0.53, 0.531, 0.2, 0.781};
        vector<double> results {0.0043887732991, 0.134444693935, 0.0170865213655, 0.212236876077, 
            0.23977370731, 7.55016097332e-08, 3.65294483313e-10, 0.104376365853, 0.11113197172, 
            1.15182183063e-05, 0.0227207304627, 0.154448930158, 0.0458548310537, 0.252217428847,
            0.0711766702131, 0.344797248656, 0.135288326179, 0.00225062122373, 0.00921270484992, 0.106502656473};
        for(unsigned int i = 0; i < ns.size(); ++i) {
            CHECK(binom_probability_gamma(ns[i], ks[i], ps[i]) == Approx(results[i]));
        }
    }

    SECTION("large numbers") {
        vector<int> ns {822, 2743, 2839, 1717, 817, 967, 3453, 508, 847, 3366};
        vector<int> ks {536, 2115, 2368, 761, 65, 776, 924, 51, 654, 895};
        vector<double> ps {0.319, 0.092, 0.301, 0.088, 0.233, 0.122, 0.267, 0.055, 0.133, 0.237};
        vector<double> results {2.81146007965e-85, 0.0, 0.0, 0.0, 3.06463896982e-31, 0.0,
            0.0152864467193, 1.62629093503e-05, 0.0, 8.05206271196e-06};

        for(unsigned int i = 0; i < ns.size(); ++i) {
            CHECK(binom_probability_gamma(ns[i], ks[i], ps[i]) == Approx(results[i]));
        }
    }
}

TEST_CASE("nucleotide distributions are computed correctly", "[computeNucleotideDistribution]") {
    SECTION("zero profiles") {
        map<Profile, int> profiles;
        array<double, 4> dist {0.25,0.25,0.25,0.25};

        array<double, 4> result = computeNucleotideDistribution(profiles);
        for(int i = 0; i < 4; ++i) {
            CHECK(result[i] == Approx(dist[i]));
        }
    }
    SECTION("one base") {
        map<Profile, int> profiles {{{10,0,0,0,10}, {1}}};
        array<double, 4> dist {1, 0, 0, 0};

        array<double, 4> result = computeNucleotideDistribution(profiles);
        for(int i = 0; i < 4; ++i) {
            CHECK(result[i] == Approx(dist[i]));
        }
    }
    SECTION("multiple bases") {
        map<Profile, int> profiles {
            {{1,0,0,0,1}, {4}},
            {{1,1,0,0,2}, {2}},
            {{0,0,0,1,1}, {2}},
        };
        array<double, 4> dist {0.6, 0.2, 0, 0.2};

        array<double, 4> result = computeNucleotideDistribution(profiles);
        for(int i = 0; i < 4; ++i) {
            CHECK(result[i] == Approx(dist[i]));
        }
    }
}
