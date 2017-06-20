#include <iostream>

#include "simplex.hpp"

Simplex2D::Simplex2D(int n, double init_x1, double init_x2, double init_stepsizes) {
    this->x = gsl_vector_alloc(2);
    gsl_vector_set(this->x, 0, init_x1);
    gsl_vector_set(this->x, 1, init_x2);

    this->step_sizes = gsl_vector_alloc(2);
    gsl_vector_set_all(this->step_sizes, init_stepsizes);

    this->minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
    this->minimize.n = 2;
}

Simplex2DResult Simplex2D::run(double (*f)(const gsl_vector*, void*), void* params) {
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
            std::cerr << "# Simplex converged in " << i << " iterations." << std::endl;
        }
    } while (status == GSL_CONTINUE && i < 1000);
    if (status == GSL_CONTINUE) {
        std::cerr << "Error: max likelihood commputation did not converge in " << i << " iterations!" << std::endl;
    }

    return Simplex2DResult {
        gsl_vector_get(this->minimizer->x, 0), 
        gsl_vector_get(this->minimizer->x, 1),
        this->minimizer->fval
    };
}

Simplex2D::~Simplex2D() {
    gsl_vector_free(x);
    gsl_vector_free(step_sizes);
    gsl_multimin_fminimizer_free(minimizer);
}
