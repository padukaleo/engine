#include "electrochem.h"
#include <math.h>

// ── Sel Elektrokimia ──────────────────────────────────

ChemResult cell_potential(chem_float E_cathode, chem_float E_anode) {
    return CHEM_OK(E_cathode - E_anode, "V");
}

ChemResult gibbs_from_cell(chem_float n, chem_float E_cell) {
    if (n <= 0.0) return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    // ΔG = -nFE → konversi J ke kJ
    return CHEM_OK((-n * FARADAY * E_cell) / 1000.0, "kJ/mol");
}

ChemResult cell_from_gibbs(chem_float delta_G0, chem_float n) {
    if (n <= 0.0) return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    // E = -ΔG / nF → konversi kJ ke J
    return CHEM_OK(-(delta_G0 * 1000.0) / (n * FARADAY), "V");
}

ChemResult keq_from_cell(chem_float E_cell, chem_float n, chem_float temp_K) {
    if (n      <= 0.0) return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // Keq = e^(nFE / RT)
    chem_float keq = exp((n * FARADAY * E_cell) / (R_ELECTRO * temp_K));
    return CHEM_OK(keq, "");
}

// ── Persamaan Nernst ──────────────────────────────────

ChemResult nernst(
    chem_float E_standard,
    chem_float n,
    chem_float Q,
    chem_float temp_K
) {
    if (n      <= 0.0) return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    if (Q      <= 0.0) return CHEM_ERR("Q harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // E = E° - (RT/nF) × ln(Q)
    chem_float E = E_standard - ((R_ELECTRO * temp_K) / (n * FARADAY)) * log(Q);
    return CHEM_OK(E, "V");
}

ChemResult nernst_25(chem_float E_standard, chem_float n, chem_float Q) {
    if (n <= 0.0) return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    if (Q <= 0.0) return CHEM_ERR("Q harus lebih dari 0");
    // E = E° - (0.0592/n) × log10(Q)
    chem_float E = E_standard - (0.0592 / n) * log10(Q);
    return CHEM_OK(E, "V");
}

// ── Elektrolisis ──────────────────────────────────────

ChemResult faraday_mass(
    chem_float I,
    chem_float t,
    chem_float Ar,
    chem_int   n
) {
    if (I  <= 0.0) return CHEM_ERR("Arus harus lebih dari 0");
    if (t  <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    if (Ar <= 0.0) return CHEM_ERR("Massa atom harus lebih dari 0");
    if (n  <= 0)   return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    // m = (I × t × Ar) / (n × F)
    return CHEM_OK((I * t * Ar) / (n * FARADAY), "gram");
}

ChemResult faraday_time(
    chem_float mass_g,
    chem_float I,
    chem_float Ar,
    chem_int   n
) {
    if (mass_g <= 0.0) return CHEM_ERR("Massa harus lebih dari 0");
    if (I      <= 0.0) return CHEM_ERR("Arus harus lebih dari 0");
    if (Ar     <= 0.0) return CHEM_ERR("Massa atom harus lebih dari 0");
    if (n      <= 0)   return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    // t = (m × n × F) / (I × Ar)
    return CHEM_OK((mass_g * n * FARADAY) / (I * Ar), "detik");
}

ChemResult faraday_gas_volume(chem_float I, chem_float t, chem_int n) {
    if (I <= 0.0) return CHEM_ERR("Arus harus lebih dari 0");
    if (t <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    if (n <= 0)   return CHEM_ERR("Jumlah elektron harus lebih dari 0");
    // V = (I × t × 22400) / (n × F) → mL di STP
    return CHEM_OK((I * t * 22400.0) / (n * FARADAY), "mL");
}

ChemResult electric_charge(chem_float I, chem_float t) {
    if (I <= 0.0) return CHEM_ERR("Arus harus lebih dari 0");
    if (t <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    return CHEM_OK(I * t, "C");
}