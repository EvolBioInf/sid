#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <map>
#include <vector>

#include "profiles.hpp"

double binom_probability(int n, int k, double p);
double binom_probability_gamma(int n, int k, double p);

std::array<double, 4> computeNucleotideDistribution(const std::vector<Profile>&, const std::vector<int>&);
std::vector<std::pair<double, double>> computeLikelihoods(const std::vector<Profile>&, const std::vector<int>&);

double profileLikelihoodHomozygous(const Profile& p, const std::array<double, 4>& nucleotide_dist, double p_error);
double profileLikelihoodHeterozygous(const Profile& p, const std::array<double, 4>& nucleotide_dist, double p_error);

#endif
