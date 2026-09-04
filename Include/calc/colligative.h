#pragma once
#include "types.h"

// ── Sifat Koligatif ───────────────────────────────────

// Kenaikan titik didih → ΔTb = i × Kb × m
ChemResult boiling_point_elevation(
    chem_float Kb,          // tetapan ebullioskopi pelarut
    chem_float molality,    // mol/kg
    chem_float i            // faktor Van't Hoff (default 1)
);

// Penurunan titik beku → ΔTf = i × Kf × m
ChemResult freezing_point_depression(
    chem_float Kf,
    chem_float molality,
    chem_float i
);

// Tekanan osmotik → π = i × M × R × T
ChemResult osmotic_pressure(
    chem_float molarity,
    chem_float temp_K,
    chem_float i
);

// Penurunan tekanan uap → ΔP = Xsolute × P°
ChemResult vapor_pressure_lowering(
    chem_float P_pure,          // tekanan uap pelarut murni
    chem_float mole_frac_solute // fraksi mol solute
);

// Tekanan uap larutan → P = Xsolvent × P°
ChemResult vapor_pressure_solution(
    chem_float P_pure,
    chem_float mole_frac_solvent
);

// ── Faktor Van't Hoff ─────────────────────────────────

// i = 1 + α(n-1)
// α = derajat ionisasi, n = jumlah ion
ChemResult vant_hoff_factor(chem_float alpha, chem_int n_ions);

// Molalitas dari ΔTb atau ΔTf
ChemResult molality_from_deltaT(chem_float deltaT, chem_float K, chem_float i);

// Mr solute dari ΔTf
// Mr = (Kf × massa_solute × 1000) / (ΔTf × massa_solvent_g)
ChemResult molar_mass_from_freezing(
    chem_float Kf,
    chem_float mass_solute_g,
    chem_float mass_solvent_g,
    chem_float delta_Tf
);