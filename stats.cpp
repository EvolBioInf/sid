#include <algorithm>
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <cfloat>

#include "stats.hpp"

using namespace std;

double inline aic(double likelihood, double num_params) {
    return 2 * num_params - 2 * log(likelihood);
}

vector<pair<double, double>> relativeLikelihoods(const vector<pair<double, double>>& likelihoods) {
    vector<pair<double, double>> relative_likelihoods (likelihoods.size());
    transform(likelihoods.begin(), likelihoods.end(), relative_likelihoods.begin(),
            [](pair<double, double> ls) -> pair<double, double> {
                double firstAIC = aic(ls.first, 2);
                double secondAIC = aic(ls.second, 2);
                if (firstAIC < secondAIC) {
                    return make_pair(1.0, exp((firstAIC - secondAIC) / 2.0));
                } else {
                    return make_pair(exp((secondAIC - firstAIC) / 2.0), 1.0);
                }
            });
    return  relative_likelihoods;
}

double likelihoodRatioTest(long double l_H0, long double l_H1) {
    if (l_H0 != 0) {
        long double chisq = -2 * (logl(l_H0) - logl(fmaxl(l_H0, l_H1)));
        // gsl works only with doubles, so we truncate
        return gsl_cdf_chisq_Q((double)chisq, 1);
    } else {
        return gsl_cdf_chisq_Q(DBL_MAX, 1);
    }
}

vector<double> likelihoodRatioTest(const vector<pair<double, double>>& likelihoods) {
    vector<double> p_values (likelihoods.size());
    transform(likelihoods.begin(), likelihoods.end(), p_values.begin(),
              [](pair<double, double> l) -> double {
                  return likelihoodRatioTest(l.first, l.second);
              });
    return p_values;
}

vector<double> adjustBonferroni(const vector<double>& p_values, size_t n) {
    if (n <= 0) {
        n = p_values.size();
    }
    vector<double> adjusted_p_values (p_values.size());
    transform(p_values.begin(), p_values.end(), adjusted_p_values.begin(),
              [&n](double p) { return p * double(n); });
    return adjusted_p_values;
}

vector<size_t> descending_sorted_indices(const vector<double>& v) {
    vector<size_t> indices (v.size(), 0);

    size_t i = 0;
    generate(indices.begin(), indices.end(), [&i](){ return i++; });

    sort(indices.begin(), indices.end(), [&v](size_t i, size_t j) { return v[i] > v[j]; });
    return indices;
}

vector<double> adjustBenjaminiHochberg(const vector<double>& p_values) {
    vector<double> adjusted_p_values (p_values.size());
    vector<size_t> sorted = descending_sorted_indices(p_values);

    size_t m = p_values.size();
    adjusted_p_values[sorted[0]] = p_values[sorted[0]];
    for(size_t i = 1; i < sorted.size(); ++i) {
        adjusted_p_values[sorted[i]] = min(adjusted_p_values[sorted[i-1]], p_values[sorted[i]] * double(m) / double(m-i));
    }
    replace_if(adjusted_p_values.begin(), adjusted_p_values.end(),
               [](double d) {return d > 1;}, 1.0);
    return adjusted_p_values;
}

