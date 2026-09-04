#pragma once
#include "types.h"
#include "parser.h"
#include "molar.h"

// N = mol /Avogadro
ChemResult mol_to_particles(chem_float mol);

// partike -> mol
ChemResult particles_to_mol(chem_float particles);

//mol -> vol gas stp
// V = mol x 22,4 L
ChemResult mol_to_volume_stp(chem_float mol);

// V = mol x 24,8 L
ChemResult mol_to_volume_satp(chem_float mol);

// vol -> stp
ChemResult volume_to_mol_stp(chem_float volume);

// vol -> SATP
ChemResult volume_to_mol_satp(chem_float vliter);

// reaksi

ChemResult mol_by_ratio(
    chem_float mol_known, // mol zat
    chem_int coef_known, // koefisien zat diketahui
    chem_int coef_target // koefisien zat dicari
);

//reaksi pembatas
int limiting_reagent(
    const chem_float* reactant_mols,
    const chem_int* coefficients,
    int count
);

ChemResult percent_mass(
    const char* formula,
    const char* element_symbol
);

ChemResult empirical_formula(
    const char** symbols,
    const chem_float* percents,
    int count,
    char* out_formula, //buffer output
    int out_size
);