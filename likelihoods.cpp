#include <array>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

#include "likelihoods.hpp"

//#define REFERENCE

using namespace std;

const double DEFAULT_STEPSIZE = 1e-4;
const double DEFAULT_PI = 1e-3;
const double DEFAULT_EPSILON = 1e-3;

vector<double> lngamma_cache(100, -1);

double inline lngamma(int x) {
    // should not happen, but just in case...
    if (x < 0) {
        return gsl_sf_lngamma(x);
    } else {
        if ((unsigned)x >= lngamma_cache.size()) {
            lngamma_cache.resize(x+1, -1.0);
        }
        if (lngamma_cache[x] < 0) {
            lngamma_cache[x] = gsl_sf_lngamma(x);
        }
        return lngamma_cache[x];
    }
}

/*
 * Source: https://arachnoid.com/binomial_probability/index.html 02.06.2017
 */
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
    for(auto i = profiles.begin(); i != profiles.end(); ++i) {
        int count;
        Profile p;
        tie(p, count) = *i;
        double l1 = profileLikelihoodHomozygous(p, nd, epsilon);
        double l2 = profileLikelihoodHeterozygous(p, nd, epsilon);
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
    // starting point
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, DEFAULT_PI);
    gsl_vector_set(x, 1, DEFAULT_EPSILON);

    gsl_vector* step_sizes = gsl_vector_alloc(2);
    gsl_vector_set_all(step_sizes, DEFAULT_STEPSIZE);

    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
    gsl_multimin_function minimize;
    minimize.n = 2;
    minimize.f = &logLikelihood;

    array<double, 4> nd = computeNucleotideDistribution(profiles);
    struct LikelihoodParams p {profiles, nd};
    minimize.params = (void*)&p;

    gsl_multimin_fminimizer_set(minimizer, &minimize, x, step_sizes);

    int i = 0;
    int status = 0;
    do {
        ++i;
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if (status != 0) {
            break;
        }

        double size = gsl_multimin_fminimizer_size(minimizer);
        status = gsl_multimin_test_size(size, 1e-5);

        if (status == GSL_SUCCESS) {
            cerr << "Converged in " << i << " iterations." << endl;
        }
        //printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", i, gsl_vector_get (minimizer->x, 0), gsl_vector_get (minimizer->x, 1), minimizer->fval, size);
    } while (status == GSL_CONTINUE && i < 1000);
    if (status == GSL_CONTINUE) {
        cerr << "Error: max likelihood commputation did not converge in " << i << " iterations!" << endl;
    }

    double pi = gsl_vector_get(minimizer->x, 0);
    double epsilon = gsl_vector_get(minimizer->x, 1);

    cout << scientific;
    cout << "# pi: " <<  pi << '\t';
    cout << "epsilon: " <<  epsilon << '\t';
    cout << "log likelihood: " << minimizer->fval << endl;

    for(auto i = profiles.begin(); i != profiles.end(); ++i) {
        int count;
        Profile p;
        tie(p, count) = *i;
        double l1 = profileLikelihoodHomozygous(p, nd, epsilon);
        double l2 = profileLikelihoodHeterozygous(p, nd, epsilon);
        cout << p << '\t' << l1 << '\t' << l2 << endl;
    } 

    gsl_vector_free(x);
    gsl_vector_free(step_sizes);
    gsl_multimin_fminimizer_free(minimizer);
}
