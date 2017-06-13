#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

#include "simplex.hpp"
#include "likelihoods.hpp"

//#define REFERENCE

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



// Source: https://arachnoid.com/binomial_probability/index.html 02.06.2017
double binom_probability_gamma(int n, int k, double p) {
    //return exp(gsl_sf_lngamma(n+1) - gsl_sf_lngamma(k+1) - gsl_sf_lngamma(n-k+1) + log(p)*k + log(1-p)*(n-k));
    return exp(lngamma(n+1) - lngamma(k+1) - lngamma(n-k+1) + log(p)*k + log(1-p)*(n-k));
}

double inline profileLikelihoodHomozygous(Profile p, array<double, 4>& nucleotide_dist, double p_error) {
#ifndef REFERENCE
    double l = 0.0;

    for (int i = 0; i < 4; ++i) {
        l += nucleotide_dist[i] * binom_probability_gamma(p[COV], p[COV] - p[i], p_error);
    }
#else
    double l = 0.0;
    double g = 0.0;

    for (int i = 0; i < 4; ++i) {
        l += nucleotide_dist[i] * pow(1-p_error, p[i]) * pow(p_error / 3.0, p[COV] - p[i]);
        g += gsl_sf_lngamma(p[i] + 1);
    }
    g = exp(gsl_sf_lngamma(p[COV] + 1) - g);
    l *= g;
#endif
    return l; 
}

double inline profileLikelihoodHeterozygous(Profile p, array<double, 4>& nucleotide_dist, double p_error) {
#ifndef REFERENCE
    double l = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i+1; j < 4; ++j) {
            double summand = 2 * nucleotide_dist[i] * nucleotide_dist[j];
            summand *= binom_probability_gamma(p[COV], p[COV] - p[i] - p[j], 2.0/3.0 * p_error);
            summand *= binom_probability_gamma(p[i] + p[j], p[i], 0.5);
            l += summand;
        }
    }
    double s = 0;
    for (int i = 0; i < 4; ++i) {
       s += nucleotide_dist[i] * nucleotide_dist[i];
    }
    l /= (1-s);
#else
    double g = 0.0;
    for (int i = 0; i < 4; ++i) {
        g += gsl_sf_lngamma(p[i] + 1);
    }
    g = exp(gsl_sf_lngamma(p[COV] + 1) - g);

    double S = 0.0;
    for (int i = 0; i < 4; ++i) {
        S += nucleotide_dist[i] * nucleotide_dist[i];
    }
    S = 1.0 - S;

    double l = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i+1; j < 4; ++j) {
            l += nucleotide_dist[i] * nucleotide_dist[j] / S * pow(0.5 - p_error / 3.0, p[i]+p[j]) * pow(p_error / 3.0, p[COV] - p[i] - p[j]);
        }
    }
    l *= g;
#endif
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
    Profile acc {0,0,0,0};
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
