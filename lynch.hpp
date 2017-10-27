#pragma once

#include <array>
#include <vector>

#include <math.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

static struct MemoizedLogGamma {
	// uncached values are indicated by -2, since the minimum of lngamma in the positive reals is about -0.12
	std::vector<double> cache = std::vector<double>(100, -2.0);

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
} log_gamma {};

typedef struct GenotypeLikelihood {
	long double L_homozygous;
	long double L_heterozygous;
} GenotypeLikelihood;

typedef struct {
	double heterozygosity;
	double error_rate;
	std::vector<GenotypeLikelihood> profile_likelihoods;
} ProfileGenotypeLikelihoods;

ProfileGenotypeLikelihoods estimateProfileGenotypeLikelihoods(const std::vector<UniqueProfile>&, const std::array<double, 4>);

double compoundLikelihood(const gsl_vector*, void*);

long double inline multinomialCoefficient(const UniqueProfile& p) {
	// compute with logGamma trick
	return expl(log_gamma(p.coverage + 1) 
		- log_gamma(p.profile[0] + 1)
		- log_gamma(p.profile[1] + 1)
		- log_gamma(p.profile[2] + 1)
		- log_gamma(p.profile[3] + 1));
}

long double inline heterozygousLikelihood(const UniqueProfile& p, double error_probability, const std::array<double, 4> nucleotide_distribution) {
        long double L = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
				L += nucleotide_distribution[i] * nucleotide_distribution[j]
					* powl((1 - 2./3. * error_probability) / 2., p.profile[i] + p.profile[j])
					* powl(error_probability / 3., p.coverage - p.profile[i] - p.profile[j]);
            }
        }
        // account for the fact that (i,i) is not a possible nucleotide configuration
        // in the heterozygous case
        long double s = 0;
        for (int i = 0; i < 4; ++i) {
            s += nucleotide_distribution[i] * nucleotide_distribution[i];
        }
        L /= (1 - s);
	return multinomialCoefficient(p) * L;
}

long double inline heterozygousLikelihood(const UniqueProfile& p, double error_probability, int ref0, int ref1) {
	return multinomialCoefficient(p)
		* powl((1 - 2./3. * error_probability) / 2., p.profile[ref0] + p.profile[ref1])
		* powl(error_probability / 3., p.coverage - p.profile[ref0] - p.profile[ref1]);
}

long double inline homozygousLikelihood(const UniqueProfile& p, double error_probability, const std::array<double, 4> nucleotide_distribution) {
    long double L = 0;
    for (int i = 0; i < 4; ++i) {
		L += nucleotide_distribution[i]
			* powl(1 - error_probability, p.profile[i])
			* powl(error_probability / 3., p.coverage - p.profile[i]);
    }
	return multinomialCoefficient(p) * L;
}

long double inline homozygousLikelihood(const UniqueProfile& p, double error_probability, int ref) {
	return multinomialCoefficient(p)
		* powl(1 - error_probability, p.profile[ref])
		* powl(error_probability / 3., p.coverage - p.profile[ref]);
}