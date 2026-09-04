#include "thermo.h"
#include <math.h>

// ── Entalpi ───────────────────────────────────────────

ChemResult delta_enthalpy(chem_float h_products, chem_float h_reactants) {
    return CHEM_OK(h_products - h_reactants, "kJ/mol");
}

ChemResult hess_law(
    const chem_float* delta_h,
    const chem_int*   signs,
    int               count
) {
    if (!delta_h || !signs || count <= 0)
        return CHEM_ERR("Input tidak valid");

    chem_float total = 0.0;
    for (int i = 0; i < count; i++) {
        if (signs[i] != 1 && signs[i] != -1)
            return CHEM_ERR("Sign harus +1 atau -1");
        total += signs[i] * delta_h[i];
    }
    return CHEM_OK(total, "kJ/mol");
}

ChemResult enthalpy_bond_energy(
    chem_float bonds_broken,
    chem_float bonds_formed
) {
    if (bonds_broken < 0.0) return CHEM_ERR("Energi ikatan tidak boleh negatif");
    if (bonds_formed < 0.0) return CHEM_ERR("Energi ikatan tidak boleh negatif");
    // ΔH = putus - terbentuk
    return CHEM_OK(bonds_broken - bonds_formed, "kJ/mol");
}

ChemResult heat_transfer(
    chem_float mass_g,
    chem_float specific_heat_val,
    chem_float delta_T
) {
    if (mass_g <= 0.0)             return CHEM_ERR("Massa harus lebih dari 0");
    if (specific_heat_val <= 0.0)  return CHEM_ERR("Kalor jenis harus lebih dari 0");
    if (delta_T == 0.0)            return CHEM_ERR("ΔT tidak boleh 0");
    return CHEM_OK(mass_g * specific_heat_val * delta_T, "J");
}

ChemResult specific_heat(
    chem_float q,
    chem_float mass_g,
    chem_float delta_T
) {
    if (mass_g  <= 0.0) return CHEM_ERR("Massa harus lebih dari 0");
    if (delta_T == 0.0) return CHEM_ERR("ΔT tidak boleh 0");
    // c = q / (m × ΔT)
    return CHEM_OK(q / (mass_g * delta_T), "J/g·K");
}

// ── Entropi ───────────────────────────────────────────

ChemResult delta_entropy(chem_float s_products, chem_float s_reactants) {
    return CHEM_OK(s_products - s_reactants, "J/mol·K");
}

// ── Energi Gibbs ──────────────────────────────────────

ChemResult gibbs_energy(
    chem_float delta_H,
    chem_float delta_S,
    chem_float temp_K
) {
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // ΔG = ΔH - T(ΔS/1000) → konversi J ke kJ
    chem_float delta_G = delta_H - (temp_K * (delta_S / 1000.0));
    return CHEM_OK(delta_G, "kJ/mol");
}

ChemResult gibbs_from_Q(
    chem_float delta_G0,
    chem_float Q,
    chem_float temp_K
) {
    if (Q      <= 0.0) return CHEM_ERR("Q harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // ΔG = ΔG° + RT ln(Q) → konversi J ke kJ
    chem_float delta_G = delta_G0 + ((R_THERMO * temp_K * log(Q)) / 1000.0);
    return CHEM_OK(delta_G, "kJ/mol");
}

ChemResult gibbs_from_Keq(chem_float Keq, chem_float temp_K) {
    if (Keq    <= 0.0) return CHEM_ERR("Keq harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // ΔG° = -RT ln(Keq) → konversi J ke kJ
    chem_float delta_G0 = -(R_THERMO * temp_K * log(Keq)) / 1000.0;
    return CHEM_OK(delta_G0, "kJ/mol");
}

ChemResult keq_from_gibbs(chem_float delta_G0, chem_float temp_K) {
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // Keq = e^(-ΔG° × 1000 / RT) → konversi kJ ke J
    chem_float keq = exp(-(delta_G0 * 1000.0) / (R_THERMO * temp_K));
    return CHEM_OK(keq, "");
}

// ── Spontanitas ───────────────────────────────────────

int spontaneity(chem_float delta_G) {
    if (delta_G < -1e-9) return  1;   // spontan
    if (delta_G >  1e-9) return -1;   // tidak spontan
    return 0;                          // setimbang
}

ChemResult equilibrium_temp(chem_float delta_H, chem_float delta_S) {
    if (delta_S == 0.0)
        return CHEM_ERR("ΔS tidak boleh 0");
    // T = ΔH(kJ) × 1000 / ΔS(J/mol·K)
    chem_float T = (delta_H * 1000.0) / delta_S;
    if (T < 0.0)
        return CHEM_ERR("Tidak ada suhu kesetimbangan (ΔH dan ΔS berlawanan tanda)");
    return CHEM_OK(T, "K");
}