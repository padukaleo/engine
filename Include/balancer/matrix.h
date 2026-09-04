#pragma once
#include "types.h"
#include "parser.h"

#define MAX_SPECIES 20
#define MAX_ELEMENTS 20

typedef struct {
    chem_float data[MAX_ELEMENTS][MAX_SPECIES +1];
    int rows; //jumlah element unik
    int cols; //jumlah senyawa
} Matrix;

int build_matrix(
    const char** reactants, int n_reactants,
    const char** products, int n_products,
    Matrix* out
);

//debug