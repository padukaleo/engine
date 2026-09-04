#pragma once
#include "types.h"
#include "parser.h"

ChemResult molar_mass(const char* formula);
ChemResult molar_mass_from_compound(const Compound* compound);
ChemResult mass_to_mol(const char* formula, chem_float mass_gram);
ChemResult mol_to_mass(const char* formula, chem_float mol);