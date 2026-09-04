#pragma once
#include "types.h"

#define FARADAY     96485.3321  // C/mol
#define R_ELECTRO   8.314       // J/mol·K

// ── Sel Elektrokimia ──────────────────────────────────

// E°cell = E°katoda - E°anoda
ChemResult cell_potential(chem_float E_cathode, chem_float E_anode);

// ΔG dari E°cell
// ΔG = -nFE
ChemResult gibbs_from_cell(chem_float n, chem_float E_cell);

// E°cell dari ΔG°
ChemResult cell_from_gibbs(chem_float delta_G0, chem_float n);

// Keq dari E°cell
// ln(Keq) = nFE° / RT
ChemResult keq_from_cell(chem_float E_cell, chem_float n, chem_float temp_K);

// ── Persamaan Nernst ──────────────────────────────────

// E = E° - (RT/nF) × ln(Q)
ChemResult nernst(
    chem_float E_standard,
    chem_float n,
    chem_float Q,
    chem_float temp_K       // default 298.15 K
);

// Nernst disederhanakan (25°C)
// E = E° - (0.0592/n) × log(Q)
ChemResult nernst_25(chem_float E_standard, chem_float n, chem_float Q);

// ── Elektrolisis (Hukum Faraday) ──────────────────────

// Massa endapan → m = (I × t × Ar) / (n × F)
ChemResult faraday_mass(
    chem_float I,       // Ampere
    chem_float t,       // detik
    chem_float Ar,      // massa atom relatif
    chem_int   n        // jumlah elektron
);

// Waktu dari massa target
// t = (m × n × F) / (I × Ar)
ChemResult faraday_time(
    chem_float mass_g,
    chem_float I,
    chem_float Ar,
    chem_int   n
);

// Volume gas di elektroda (STP)
// V = (I × t × 22400) / (n × F)  → mL
ChemResult faraday_gas_volume(
    chem_float I,
    chem_float t,
    chem_int   n
);

// Muatan listrik → Q = I × t
ChemResult electric_charge(chem_float I, chem_float t);