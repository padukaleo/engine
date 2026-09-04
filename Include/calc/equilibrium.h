#pragma once
#include "types.h"

// ── Tetapan Kesetimbangan ─────────────────────────────

// Kc dari konsentrasi produk & reaktan
// Kc = [produk]^koef / [reaktan]^koef
ChemResult calc_Kc(
    const chem_float* conc_products,
    const chem_int*   coef_products,
    int               n_products,
    const chem_float* conc_reactants,
    const chem_int*   coef_reactants,
    int               n_reactants
);

// Kp dari tekanan parsial
ChemResult calc_Kp(
    const chem_float* press_products,
    const chem_int*   coef_products,
    int               n_products,
    const chem_float* press_reactants,
    const chem_int*   coef_reactants,
    int               n_reactants
);

// Konversi Kc ↔ Kp
// Kp = Kc × (RT)^Δn
ChemResult Kp_from_Kc(chem_float Kc, chem_float temp_K, chem_int delta_n);
ChemResult Kc_from_Kp(chem_float Kp, chem_float temp_K, chem_int delta_n);

// ── Reaction Quotient ─────────────────────────────────

// Qc — sama struktur dengan Kc
ChemResult calc_Qc(
    const chem_float* conc_products,
    const chem_int*   coef_products,
    int               n_products,
    const chem_float* conc_reactants,
    const chem_int*   coef_reactants,
    int               n_reactants
);

// Bandingkan Q dengan K
// return  1 = reaksi maju (Q < K)
// return  0 = setimbang   (Q = K)
// return -1 = reaksi balik (Q > K)
int reaction_direction(chem_float Q, chem_float K);

// ── Derajat Ionisasi ──────────────────────────────────

// α = [H+] / Ca → dari Ka dan konsentrasi
ChemResult ionization_degree_eq(chem_float Ka, chem_float conc);

// ── Kesetimbangan Heterogen ───────────────────────────

// Kesetimbangan kelarutan → sudah di solution.h (Ksp)
// Di sini: ion product vs Ksp
// return  1 = endapan terbentuk (Q > Ksp)
// return  0 = tepat jenuh
// return -1 = belum jenuh
int precipitation_check(chem_float Q, chem_float Ksp);