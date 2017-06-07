#ifndef LIKELIHOODS_H
#define LIKELIHOODS_H

#include <map>

#include "profiles.hpp"

double binom_probability(int n, int k, double p);
double binom_probability_gamma(int n, int k, double p);

std::array<double, 4> computeNucleotideDistribution(std::map<Profile, int>& profiles);
void computeLikelihoods(std::map<Profile, int> profiles);

#endif
