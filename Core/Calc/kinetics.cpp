#include "kinetics.h"
#include <math.h>

#define R_KIN   8.314   // J/mol·K
#define LN2     0.693147180559945

// ── Konstanta Laju ────────────────────────────────────

ChemResult rate_constant(chem_float A, chem_float Ea, chem_float temp_K) {
    if (A      <= 0.0) return CHEM_ERR("Faktor frekuensi harus lebih dari 0");
    if (Ea     <= 0.0) return CHEM_ERR("Energi aktivasi harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    return CHEM_OK(A * exp(-Ea / (R_KIN * temp_K)), "");
}

ChemResult activation_energy(
    chem_float k1, chem_float T1,
    chem_float k2, chem_float T2
) {
    if (k1 <= 0.0 || k2 <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (T1 <= 0.0 || T2 <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    if (T1 == T2) return CHEM_ERR("T1 dan T2 tidak boleh sama");
    // Ea = R × ln(k2/k1) / (1/T1 - 1/T2)
    chem_float Ea = R_KIN * log(k2 / k1) / ((1.0 / T1) - (1.0 / T2));
    if (Ea <= 0.0) return CHEM_ERR("Ea tidak valid, periksa input");
    return CHEM_OK(Ea, "J/mol");
}

// ── Waktu Paruh ───────────────────────────────────────

ChemResult half_life_zero(chem_float conc_initial, chem_float k) {
    if (conc_initial <= 0.0) return CHEM_ERR("Konsentrasi awal harus lebih dari 0");
    if (k            <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    return CHEM_OK(conc_initial / (2.0 * k), "s");
}

ChemResult half_life_first(chem_float k) {
    if (k <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    return CHEM_OK(LN2 / k, "s");
}

ChemResult half_life_second(chem_float k, chem_float conc_initial) {
    if (k            <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (conc_initial <= 0.0) return CHEM_ERR("Konsentrasi awal harus lebih dari 0");
    return CHEM_OK(1.0 / (k * conc_initial), "s");
}

// ── Konsentrasi vs Waktu ──────────────────────────────

ChemResult conc_zero_order(chem_float c0, chem_float k, chem_float t) {
    if (c0 <= 0.0) return CHEM_ERR("Konsentrasi awal harus lebih dari 0");
    if (k  <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (t  <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    chem_float c = c0 - k * t;
    if (c < 0.0) return CHEM_ERR("Reaktan sudah habis sebelum waktu t");
    return CHEM_OK(c, "mol/L");
}

ChemResult conc_first_order(chem_float c0, chem_float k, chem_float t) {
    if (c0 <= 0.0) return CHEM_ERR("Konsentrasi awal harus lebih dari 0");
    if (k  <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (t  <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    return CHEM_OK(c0 * exp(-k * t), "mol/L");
}

ChemResult conc_second_order(chem_float c0, chem_float k, chem_float t) {
    if (c0 <= 0.0) return CHEM_ERR("Konsentrasi awal harus lebih dari 0");
    if (k  <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (t  <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    // 1/[A] = 1/[A]0 + kt
    return CHEM_OK(1.0 / ((1.0 / c0) + k * t), "mol/L");
}

// ── Laju Reaksi ───────────────────────────────────────

ChemResult reaction_rate(
    chem_float k,
    chem_float conc_A, chem_int order_A,
    chem_float conc_B, chem_int order_B
) {
    if (k      <= 0.0) return CHEM_ERR("k harus lebih dari 0");
    if (conc_A <= 0.0) return CHEM_ERR("Konsentrasi A harus lebih dari 0");
    if (order_A < 0)   return CHEM_ERR("Orde reaksi tidak boleh negatif");

    chem_float rate = k * pow(conc_A, order_A);

    if (order_B > 0) {
        if (conc_B <= 0.0) return CHEM_ERR("Konsentrasi B harus lebih dari 0");
        rate *= pow(conc_B, order_B);
    }

    return CHEM_OK(rate, "mol/L·s");
}