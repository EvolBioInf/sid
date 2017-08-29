#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include <gsl/gsl_multimin.h>
#include <functional>
#include <array>

template <int N>
struct FunctionMinimizerResult {
    std::array<double, N> x;
    double fval;

    FunctionMinimizerResult(std::array<double, N> x, double fval)
        : x{x}, fval{fval} {}
};

template <int N>
class FunctionMinimizer {
    private:
        gsl_vector* x;
        gsl_vector* step_sizes;
        gsl_multimin_fminimizer* minimizer;
        gsl_multimin_function minimize;
    public:
        FunctionMinimizer(std::array<double, N> init_x, std::array<double, N> init_stepsizes);
        ~FunctionMinimizer();
        FunctionMinimizerResult<N> run(double (*f)(const gsl_vector*, void*), void* params);
};

using namespace std;

template <int N>
FunctionMinimizer<N>::FunctionMinimizer(array<double, N> init_x, array<double, N> init_stepsizes) {
    this->x = gsl_vector_alloc(N);
    for (int i = 0; i < N; ++i) {
        gsl_vector_set(this->x, i, init_x[i]);
    }
    this->step_sizes = gsl_vector_alloc(N);
    for (int i = 0; i < N; ++i) {
        gsl_vector_set(this->step_sizes, i, init_stepsizes[i]);
    }

    this->minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, N);
    this->minimize.n = N;
}

template <int N>
FunctionMinimizerResult<N> FunctionMinimizer<N>::run(double (*f)(const gsl_vector*, void*), void* params) {
    minimize.f = f;
    minimize.params = params;

    gsl_multimin_fminimizer_set(this->minimizer, &this->minimize, this->x, this->step_sizes);

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
            std::cerr << "# Function minimization converged in " << i << " iterations." << std::endl;
        }
    } while (status == GSL_CONTINUE && i < 1000);
    if (status == GSL_CONTINUE) {
        std::cerr << "Error: function minimization did not converge in " << i << " iterations!" << std::endl;
    }

    array<double, N> result;
    for (int i = 0; i < N; ++i) {
        result[i] = gsl_vector_get(this->minimizer->x, i);
    }

    return FunctionMinimizerResult<N> {
        result,
        this->minimizer->fval
    };
}

template <int N>
FunctionMinimizer<N>::~FunctionMinimizer() {
    gsl_vector_free(x);
    gsl_vector_free(step_sizes);
    gsl_multimin_fminimizer_free(minimizer);
}

#endif
