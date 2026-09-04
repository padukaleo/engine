#pragma once
#include "types.h"
#include "matrix.h"
#include "gaussian.h"

BalancerResult balance_equation(const char* equation);
BalancerResult balance_equation_parts(
    const char** reactants, int n_reactants,
    const char** products, int n_products
);