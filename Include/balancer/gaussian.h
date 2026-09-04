#pragma once
#include "matrix.h"

void gaussian_rref(Matrix* m);

int extract_solution {
    const Matrix* m,
    int* coefficients,
    int n_species
}