#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <map>
#include <vector>

#include "profiles.hpp"

double binom_probability(int n, int k, double p);
double binom_probability_gamma(int n, int k, double p);

std::array<double, 4> computeNucleotideDistribution(const std::vector<Profile>&, const std::vector<int>&);

typedef struct GenomeParameters {
    double heterozygosity;
    double error_rate;
    std::vector<double> hom_likelihoods;
    std::vector<double> het_likelihoods;
    GenomeParameters (double h, double e) : heterozygosity{h}, error_rate{e} {}
} GenomeParameters;

GenomeParameters estimateGenomeParameters(const std::vector<Profile>& profiles, const std::vector<int>& nucleotide_counts);

double profileLikelihoodHomozygous(const Profile& p, const std::array<double, 4>& nucleotide_dist, double p_error);
double profileLikelihoodHeterozygous(const Profile& p, const std::array<double, 4>& nucleotide_dist, double p_error);

double profileLikelihoodHomozygous(const Profile& p, double p_error, int i);
double profileLikelihoodHeterozygous(const Profile& p, double p_error, int i, int j);
#endif
