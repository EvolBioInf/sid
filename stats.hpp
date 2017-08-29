#ifndef STATS_H
#define STATS_H

#include <vector>

std::vector<std::pair<double, double>> relativeLikelihoods(const std::vector<std::pair<double, double>>&);

double likelihoodRatioTest(long double l_H0, long double l_H1);
std::vector<double> likelihoodRatioTest(const std::vector<std::pair<double, double>>&);
std::vector<double> adjustBonferroni(const std::vector<double>&, int n = 0);
std::vector<double> adjustBenjaminiHochberg(const std::vector<double>&);

std::vector<size_t> descending_sorted_indices(const std::vector<double>&);
#endif
