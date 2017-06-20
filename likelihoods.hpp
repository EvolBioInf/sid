#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <map>
#include <vector>

#include "profiles.hpp"

double binom_probability(int n, int k, double p);
double binom_probability_gamma(int n, int k, double p);

std::array<double, 4> computeNucleotideDistribution(const std::vector<Profile>&, const std::vector<int>&);
std::vector<std::pair<double, double>> computeLikelihoods(const std::vector<Profile>&, const std::vector<int>&);

#endif
