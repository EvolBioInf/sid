#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <gsl/gsl_sf_gamma.h>

#include "likelihoods.hpp"
#include "optimization.hpp"

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

long double profileLikelihoodHomozygous(const Profile& p, double p_error, int i) {
    long double l = powl(1 - p_error, p[i]) * powl(p_error / 3, p[COV] - p[i]);

    // compute multiomial coefficient with logGamma trick
    l *= expl(lngamma(p[COV] + 1) - lngamma(p[A] + 1) - lngamma(p[C] + 1) - lngamma(p[G] + 1) - lngamma(p[T] + 1));
    return l;
}

long double profileLikelihoodHeterozygous(const Profile& p, double p_error, int i, int j) {
    long double l = powl((1.0 - 2.0*p_error/3.0)/2.0, p[i] + p[j]) * powl(p_error/3.0, p[COV] - p[i] - p[j]);

    // compute multiomial coefficient with logGamma trick
    l *= expl(lngamma(p[COV] + 1) - lngamma(p[A] + 1) - lngamma(p[C] + 1) - lngamma(p[G] + 1) - lngamma(p[T] + 1));
    return l;
}

long double profileLikelihoodHomozygous(const Profile& p, const array<double, 4>& nucleotide_dist, double p_error) {
    long double l = 0.0;

    for (int i = 0; i < 4; ++i) {
        l += nucleotide_dist[i] * powl(1 - p_error, p[i]) * powl(p_error / 3, p[COV] - p[i]);
    }
    // compute multiomial coefficient with logGamma trick
    l *= expl(lngamma(p[COV] + 1) - lngamma(p[A] + 1) - lngamma(p[C] + 1) - lngamma(p[G] + 1) - lngamma(p[T] + 1));

    return l;
}

long double profileLikelihoodHeterozygous(const Profile& p, const array<double, 4>& nucleotide_dist, double p_error) {
    long double l = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i+1; j < 4; ++j) {
            l += nucleotide_dist[i] * nucleotide_dist[j] * powl((1.0 - 2.0*p_error/3.0)/2.0, p[i] + p[j]) * powl(p_error/3.0, p[COV] - p[i] - p[j]);
        }
    }
    // compute multiomial coefficient with logGamma trick
    l *= expl(lngamma(p[COV] + 1) - lngamma(p[A] + 1) - lngamma(p[C] + 1) - lngamma(p[G] + 1) - lngamma(p[T] + 1));

    long double s = 0.0;
    for (int i = 0; i < 4; ++i) {
       s += nucleotide_dist[i] * nucleotide_dist[i];
    }
    l /= (1 - s);

    return l;
}

struct LikelihoodParams {
    const vector<Profile>& profiles;
    const vector<int>& counts;
    const array<double, 4>& nucleotide_dist;

    LikelihoodParams(const vector<Profile>& ps, const vector<int>& cs, const array<double, 4>& nd)
        : profiles{ps}, counts{cs}, nucleotide_dist{nd} {}
};

double logLikelihood(const gsl_vector* v, void *params_) {
    double pi = gsl_vector_get(v, 0);
    double epsilon = gsl_vector_get(v, 1);

    struct LikelihoodParams* params = (struct LikelihoodParams*)params_;
    const vector<Profile>& profiles = params->profiles;
    const vector<int>& counts = params->counts;
    const array<double, 4>& nd = params->nucleotide_dist;

    if (pi < 0 || pi > 1 || epsilon < 0 || epsilon > 1) {
        return numeric_limits<double>::max();
    }

    long double likelihood = 0;
    // for (const auto& [profile, count] : profiles) {
    for (int i = 0; i < profiles.size(); ++i) {
        long double l1 = profileLikelihoodHomozygous(profiles[i], nd, epsilon);
        long double l2 = profileLikelihoodHeterozygous(profiles[i], nd, epsilon);
        long double l = (1.0 - pi) * l1 + pi * l2;
        if (l > 0) {
            likelihood += log(l) * counts[i];
        }
    }
    return (double)-likelihood;
}

array<double, 4> computeNucleotideDistribution(const std::vector<Profile>& profiles, const std::vector<int>& counts) {
    array<unsigned long int, 5> acc {0,0,0,0,0};
    for (int i = 0; i < profiles.size(); ++i) {
        for (auto j : {A,C,G,T,COV}) {
            acc[j] += counts[i] * profiles[i][j];
        }
    }

    if (acc[COV] != 0) {
        return {
            (double)acc[A] / acc[COV],
            (double)acc[C] / acc[COV],
            (double)acc[G] / acc[COV],
            (double)acc[T] / acc[COV]};
    } else {
        return {0.25,0.25,0.25,0.25};
    }
}

GenomeParameters estimateGenomeParameters(const std::vector<Profile>& profiles, const std::vector<int>& nucleotide_counts) {
    array<double, 4> nd = computeNucleotideDistribution(profiles, nucleotide_counts);
    struct LikelihoodParams params {profiles, nucleotide_counts, nd};

    FunctionMinimizer<2> mini {{DEFAULT_PI, DEFAULT_EPSILON}, {DEFAULT_STEPSIZE, DEFAULT_STEPSIZE}};
    FunctionMinimizerResult<2> result = mini.run(logLikelihood, (void*)&params);
    double pi = result.x[0];
    double epsilon = result.x[1];
    double logL = result.fval;

    cerr << scientific;
    cerr << "# pi: " <<  pi << '\t';
    cerr << "epsilon: " <<  epsilon << '\t';
    cerr << "log likelihood: " << logL << endl;

    GenomeParameters gp = GenomeParameters(pi, epsilon);
    for (const Profile& profile : profiles) {
        long double l1 = profileLikelihoodHomozygous(profile, nd, epsilon);
        gp.hom_likelihoods.emplace_back(l1);
        long double l2 = profileLikelihoodHeterozygous(profile, nd, epsilon);
        gp.het_likelihoods.emplace_back(l2);
    }
    return gp;
}
