#pragma once
#include "types.h"

ChemResult calc_kc(
    const chem_float* conc_products,
    const chem_int* coef_products,
    int n_products,
    const chem_float* conc_reactants,
    const chem_int* coef_reactants,
    int n_reactants
);

ChemResult calc_kp(
    const chem_float* press_products,
    const chem_int* coef_products,
    int n_products,
    const chem_float* press_reactants,
    const chem_int* coef_reactants,
    int n_reactants
);

ChemResult Kp_from_Kc(chem_float Kc, chem_float temp_K, chem_int delta_n);
ChemResult Kc_from_Kp(chem_float Kp, chem_float tempt_K, chem_int delta_n);

ChemResult calc_Qc(
    const chem_float* conc_products,
    const chem_int* coef_products,
    int n_products,
    const chem_float* conc_reactants,
    const chem_int* coef_reactants,
    int n_reactants
);

int reaction_directions(chem_float Q, chem_float K);

ChemResult ionization_degree_eq(chem_float Ka, chem_float conc);
int preciptation_check(chem_float Q, chem_float Ksp);