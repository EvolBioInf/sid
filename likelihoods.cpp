#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

#include "simplex.hpp"
#include "likelihoods.hpp"

using namespace std;

const double DEFAULT_STEPSIZE = 1e-4;
const double DEFAULT_PI = 1e-3;
const double DEFAULT_EPSILON = 1e-3;

struct MemoizedLogGamma {
    // uncached values are indicated by -2, since the minimum of lngamma in the positive reals is about -0.12
    vector<double> cache = vector<double>(100, -2.0);

    double operator()(int x) {
        // should not happen, but just in case...
        if (x < 0) {
            return gsl_sf_lngamma(x);
        } else if (x == 0) {
	    return 0;
	} else {
            if ((unsigned)x >= cache.size()) {
                cache.resize(x+1, -2.0);
            }
            if (cache[x] < -1) {
                cache[x] = gsl_sf_lngamma(x);
            }
            return cache[x];
        }
    };
} lngamma {};

double inline multinomialCoefficient(int n, int n1, int n2, int n3, int n4) {
    return exp(lngamma(n+1) - lngamma(n1+1) - lngamma(n2+1) - lngamma(n3+1) - lngamma(n4+1));
}

double inline profileLikelihoodHomozygous(Profile p, array<double, 4>& nucleotide_dist, double p_error) {
    double l = 0.0;

    for (int i = 0; i < 4; ++i) {
	l += nucleotide_dist[i] * pow(1 - p_error, p[i]) * pow(p_error / 3, p[COV] - p[i]);
    }
    l *= multinomialCoefficient(p[COV], p[A], p[C], p[G], p[T]);

    return l; 
}

double inline profileLikelihoodHeterozygous(Profile p, array<double, 4>& nucleotide_dist, double p_error) {
    double l = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i+1; j < 4; ++j) {
	    l += nucleotide_dist[i] * nucleotide_dist[j] * pow((1.0 - 2.0*p_error/3.0)/2.0, p[i] + p[j]) * pow(p_error/3.0, p[COV] - p[i] - p[j]);
	}
    }
    l *= multinomialCoefficient(p[COV], p[A], p[C], p[G], p[T]);

    double s = 0.0;
    for (int i = 0; i < 4; ++i) {
       s += nucleotide_dist[i] * nucleotide_dist[i];
    }
    l /= (1 - s);

    return l;
}

struct LikelihoodParams {
    map<Profile, int>& profiles;
    array<double, 4> nucleotide_dist;

    LikelihoodParams(map<Profile, int>& p, array<double, 4> nd) : profiles{p}, nucleotide_dist{nd} {}
};

double logLikelihood(const gsl_vector* v, void *params_) {
    double pi = gsl_vector_get(v, 0);
    double epsilon = gsl_vector_get(v, 1);

    struct LikelihoodParams* params = (struct LikelihoodParams*)params_;
    map<Profile, int>& profiles = params->profiles;
    array<double, 4>& nd = params->nucleotide_dist;

    if (pi < 0 || pi > 1 || epsilon < 0 || epsilon > 1) {
        return numeric_limits<double>::max();
    }

    double likelihood = 0;
    for (const auto& [profile, count] : profiles) {
        double l1 = profileLikelihoodHomozygous(profile, nd, epsilon);
        double l2 = profileLikelihoodHeterozygous(profile, nd, epsilon);
        double l = (1.0 - pi) * l1 + pi * l2;
        if (l > 0) {
            likelihood += log(l) * count;
        }
    } 
    return -likelihood;
}

array<double, 4> computeNucleotideDistribution(map<Profile, int>& profiles) {
    Profile acc {0,0,0,0,0};
    acc = accumulate(profiles.begin(), profiles.end(), acc, 
            [](Profile acc, pair<Profile, int> p) {return acc + p.first * p.second;});

    int n = acc[COV];
    if (n != 0) {
        return {(double)acc[A]/n, (double)acc[C]/n, (double)acc[G]/n, (double)acc[T]/n};
    } else {
        return {0.25,0.25,0.25,0.25};
    }
}

void computeLikelihoods(map<Profile, int> profiles) {
    array<double, 4> nd = computeNucleotideDistribution(profiles);
    struct LikelihoodParams params {profiles, nd};

    Simplex2D simplex {2, DEFAULT_PI, DEFAULT_EPSILON, DEFAULT_STEPSIZE};
    Simplex2DResult result = simplex.run(logLikelihood, (void*)&params);
    double pi = result.x1;
    double epsilon = result.x2;
    double logL = result.fval;

    cout << scientific;
    cout << "# pi: " <<  pi << '\t';
    cout << "epsilon: " <<  epsilon << '\t';
    cout << "log likelihood: " << logL << endl;

    for (const auto& [profile, count] : profiles) {
        double l1 = profileLikelihoodHomozygous(profile, nd, epsilon);
        double l2 = profileLikelihoodHeterozygous(profile, nd, epsilon);
        cout << profile << '\t' << l1 << '\t' << l2 << endl;
    } 
}
