#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>
#include <functional>

struct Simplex2DResult {
    double x1;
    double x2;
    double fval;

    Simplex2DResult(double x1, double x2, double fval) : x1{x1}, x2{x2}, fval{fval} {}
};

class Simplex2D {
    private:
        gsl_vector* x;
        gsl_vector* step_sizes;
        gsl_multimin_fminimizer* minimizer;
        gsl_multimin_function minimize;
    public:
        Simplex2D(int n, double init_x1, double init_x2, double init_stepsizes);
        ~Simplex2D();
        Simplex2DResult run(double (*f)(const gsl_vector*, void*), void* params);
};

#endif
