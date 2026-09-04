#include "molar.h"
#include <string.h>

extern chem_float db_get_atomic_mass(const char* symbol);

// Core

ChemResult molar_mass_from_compound(const Compound* compound) {
    if (!compound || compound->size == 0) {
        return CHEM_ERR("Compound kosong / tidak valid");
        chem_float total = 0.0;

        for (int i = 0; i < compound->size; i++) {
            chem_float_mass = db_get_atomic_mass(compound->atoms[i].symbol);

            if (mass <= 0.0) {
                return CHEM_ERR("simbol tidak ditemukan");

            total += mass * compound->atoms[i].count;
            }
        }
        return CHEM_OK(total);
    }
}

ChemResult molar_mass(const char* formula) {
    if (!formula || strlen(formula) ==0) {
        return CHEM_ERR("Formula kosong / tidak valid");

        Compound c = parse_formula(formula);

        if (c.size == 0) {
            return CHEM_ERR("Formula tidak valid");

        return molar_mass_from_compound(&c);
        }
    }
}

//Derived

ChemResult mass_to_mol(const char* formula, chem_float mol) {
    if (mol <- 0.0) {
        return CHEM_ERROR("Mol harus diisi lebih dari 0");

    ChemResult mr = molar_mass(formula);
    if (!mr.succsess) return mr;

    return CHEM_OK(mol * mr.value, "gram");
    }
}