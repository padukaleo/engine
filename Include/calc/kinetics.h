#pragma once
#include "types.h"

// ── Konstanta Laju (Arrhenius) ────────────────────────

// k = A × e^(-Ea/RT)
ChemResult rate_constant(
    chem_float A,       // faktor frekuensi
    chem_float Ea,      // energi aktivasi J/mol
    chem_float temp_K
);

// Ea dari dua k dan dua T
// ln(k2/k1) = (Ea/R) × (1/T1 - 1/T2)
ChemResult activation_energy(
    chem_float k1, chem_float T1,
    chem_float k2, chem_float T2
);

// ── Waktu Paruh ───────────────────────────────────────

// Orde 0 → t½ = [A]0 / 2k
ChemResult half_life_zero(chem_float conc_initial, chem_float k);

// Orde 1 → t½ = ln2 / k
ChemResult half_life_first(chem_float k);

// Orde 2 → t½ = 1 / (k[A]0)
ChemResult half_life_second(chem_float k, chem_float conc_initial);

// ── Konsentrasi vs Waktu ──────────────────────────────

// Orde 0 → [A] = [A]0 - kt
ChemResult conc_zero_order(chem_float conc_initial, chem_float k, chem_float t);

// Orde 1 → [A] = [A]0 × e^(-kt)
ChemResult conc_first_order(chem_float conc_initial, chem_float k, chem_float t);

// Orde 2 → 1/[A] = 1/[A]0 + kt
ChemResult conc_second_order(chem_float conc_initial, chem_float k, chem_float t);

// ── Laju Reaksi ───────────────────────────────────────

// v = k × [A]^m × [B]^n
// Untuk 1-2 reaktan
ChemResult reaction_rate(
    chem_float k,
    chem_float conc_A, chem_int order_A,
    chem_float conc_B, chem_int order_B  // order_B = 0 kalau 1 reaktan
);