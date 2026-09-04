#pragma once
#include "types.h"

// Molaritas -> mol
ChemResult molarity(chem_float mol, chem_float volume_L);

// Molalitas -> mol
ChemResult molaliy(chem_float mol, chem_float solvent_kg);

// Normalitas -> grek
ChemResult normality(chem_float grek, chem_float volume_L);

//Fraksimol
ChemResult mole_fraction(chem_float mol_solute chem_float mol_solvent);

// persen massa
ChemResult mass_percet(chem_float mol_solute, chem_float mol_solvent);

//pengenceran

//V1N1 = V2M2
//cari V2
ChemResult dilution_volume(
    chem_float M1, chem_float V1,
    chem_float M2
)

//cari m2
ChemResult dilution_volume(
    chem_float M1, chem_float V1,
    chem_float V2
)

//titrasi
ChemResult titration_volume(
    chem_float Ma, chem_float Va, chem_int a
    chem_float Vb, chem_init b
)

ChemResult titrastion_molarity(
    chem_float Ma, chem_float Va, chem_int a
    chem_float Vb, chem_int b
)

// kelarutan
ChemResult solubility_from_ksp(chem_float ksp, chem_int cation_coef, chem_int anion_coef);

//ksp dari lauran
ChemResult ksp_from_solubility(chem_float s, chem_int cation_coef, chem_int anion_coef);