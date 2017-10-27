#include <limits>

#include "optimization.hpp"
#include "pileup.hpp"

#include "lynch.hpp"

const double DEFAULT_STEPSIZE = 1e-4;
const double DEFAULT_PI = 1e-3;
const double DEFAULT_EPSILON = 1e-3;

typedef struct {
    const std::vector<UniqueProfile>& profiles;
    const std::array<double, 4> nucleotide_distribution;
} LikelihoodFunctionParameters;

ProfileGenotypeLikelihoods estimateProfileGenotypeLikelihoods(const std::vector<UniqueProfile>& profiles, const std::array<double, 4> nucleotide_distribution) {
    LikelihoodFunctionParameters params {profiles, nucleotide_distribution};

    FunctionMinimizer<2> minimizer {{DEFAULT_PI, DEFAULT_EPSILON}, {DEFAULT_STEPSIZE, DEFAULT_STEPSIZE}};
    FunctionMinimizerResult<2> result = minimizer.run(compoundLikelihood, (void*)&params);
    double pi = result.x[0];
    double epsilon = result.x[1];
    // double logL = result.fval;

    std::vector<GenotypeLikelihood> likelihoods;
    likelihoods.reserve(profiles.size());
    for (const auto& p : profiles) {
        likelihoods.push_back({
            homozygousLikelihood(p, epsilon, nucleotide_distribution),
            heterozygousLikelihood(p, epsilon, nucleotide_distribution)
        });
    }
    return {pi, epsilon, likelihoods};
}

double compoundLikelihood(const gsl_vector* v, void* parameter_data) {
    double pi = gsl_vector_get(v, 0);
    double epsilon = gsl_vector_get(v, 1);

    if (pi < 0 || pi > 1 || epsilon < 0 || epsilon > 1) {
        return numeric_limits<double>::max();
    }
    LikelihoodFunctionParameters* parameters = static_cast<LikelihoodFunctionParameters*>(parameter_data);
    long double logLikelihood = 0;
    for (const UniqueProfile& p : parameters->profiles) {
        long double L = (1. - pi) * homozygousLikelihood(p, epsilon, parameters->nucleotide_distribution)
            + pi * heterozygousLikelihood(p, epsilon, parameters->nucleotide_distribution);
        if (L > 0) {
            logLikelihood += log(L) * p.count;
        }
    }
    return static_cast<double>(-logLikelihood);
}