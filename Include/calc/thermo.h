#pragma once
#include "types.h"

#define R_THERMO    8.314   // J/mol·K

// ── Entalpi ───────────────────────────────────────────

// ΔH = ΣH produk - ΣH reaktan
ChemResult delta_enthalpy(
    chem_float h_products,
    chem_float h_reactants
);

// Hukum Hess → ΔH total dari array ΔH tiap reaksi
// signs = +1 atau -1 (reaksi dibalik)
ChemResult hess_law(
    const chem_float* delta_h,
    const chem_int*   signs,
    int               count
);

// Entalpi dari energi ikatan
// ΔH = ΣE ikatan putus - ΣE ikatan terbentuk
ChemResult enthalpy_bond_energy(
    chem_float bonds_broken,    // total energi ikatan putus (kJ)
    chem_float bonds_formed     // total energi ikatan terbentuk (kJ)
);

// q = m × c × ΔT
ChemResult heat_transfer(
    chem_float mass_g,
    chem_float specific_heat,   // J/g·K
    chem_float delta_T          // Kelvin atau Celsius (sama saja)
);

// Cari c dari q, m, ΔT
ChemResult specific_heat(
    chem_float q,
    chem_float mass_g,
    chem_float delta_T
);

// ── Entropi ───────────────────────────────────────────

// ΔS = ΣS produk - ΣS reaktan
ChemResult delta_entropy(
    chem_float s_products,
    chem_float s_reactants
);

// ── Energi Gibbs ──────────────────────────────────────

// ΔG = ΔH - TΔS
// deltaH dalam kJ/mol, deltaS dalam J/mol·K
ChemResult gibbs_energy(
    chem_float delta_H,     // kJ/mol
    chem_float delta_S,     // J/mol·K
    chem_float temp_K
);

// ΔG dari ΔG° dan Q
// ΔG = ΔG° + RT ln(Q)
ChemResult gibbs_from_Q(
    chem_float delta_G0,    // kJ/mol
    chem_float Q,           // reaction quotient
    chem_float temp_K
);

// ΔG° dari Keq
// ΔG° = -RT ln(Keq)
ChemResult gibbs_from_Keq(chem_float Keq, chem_float temp_K);

// Keq dari ΔG°
// Keq = e^(-ΔG° / RT)
ChemResult keq_from_gibbs(chem_float delta_G0, chem_float temp_K);

// ── Spontanitas ───────────────────────────────────────

// return:
//  1  = spontan (ΔG < 0)
//  0  = setimbang (ΔG = 0)
// -1  = tidak spontan (ΔG > 0)
int spontaneity(chem_float delta_G);

// Suhu kesetimbangan → T = ΔH / ΔS
ChemResult equilibrium_temp(
    chem_float delta_H,     // kJ/mol
    chem_float delta_S      // J/mol·K
);