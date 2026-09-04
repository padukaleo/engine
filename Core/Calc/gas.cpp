#include "gas.h"

// ── Hukum Gas Ideal ───────────────────────────────────

ChemResult ideal_gas_pressure(chem_float n, chem_float T, chem_float V) {
    if (n <= 0.0) return CHEM_ERR("Mol harus lebih dari 0");
    if (T <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    if (V <= 0.0) return CHEM_ERR("Volume harus lebih dari 0");
    // P = nRT / V
    return CHEM_OK((n * R_LATM * T) / V, "atm");
}

ChemResult ideal_gas_volume(chem_float n, chem_float T, chem_float P) {
    if (n <= 0.0) return CHEM_ERR("Mol harus lebih dari 0");
    if (T <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    if (P <= 0.0) return CHEM_ERR("Tekanan harus lebih dari 0");
    // V = nRT / P
    return CHEM_OK((n * R_LATM * T) / P, "L");
}

ChemResult ideal_gas_mol(chem_float P, chem_float V, chem_float T) {
    if (P <= 0.0) return CHEM_ERR("Tekanan harus lebih dari 0");
    if (V <= 0.0) return CHEM_ERR("Volume harus lebih dari 0");
    if (T <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // n = PV / RT
    return CHEM_OK((P * V) / (R_LATM * T), "mol");
}

ChemResult ideal_gas_temp(chem_float P, chem_float V, chem_float n) {
    if (P <= 0.0) return CHEM_ERR("Tekanan harus lebih dari 0");
    if (V <= 0.0) return CHEM_ERR("Volume harus lebih dari 0");
    if (n <= 0.0) return CHEM_ERR("Mol harus lebih dari 0");
    // T = PV / nR
    return CHEM_OK((P * V) / (n * R_LATM), "K");
}

// ── Hukum Gas Gabungan ────────────────────────────────

ChemResult combined_gas_P2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float V2, chem_float T2
) {
    if (P1 <= 0.0 || V1 <= 0.0 || T1 <= 0.0 || V2 <= 0.0 || T2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    // P2 = P1V1T2 / T1V2
    return CHEM_OK((P1 * V1 * T2) / (T1 * V2), "atm");
}

ChemResult combined_gas_V2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float P2, chem_float T2
) {
    if (P1 <= 0.0 || V1 <= 0.0 || T1 <= 0.0 || P2 <= 0.0 || T2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    // V2 = P1V1T2 / T1P2
    return CHEM_OK((P1 * V1 * T2) / (T1 * P2), "L");
}

ChemResult combined_gas_T2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float P2, chem_float V2
) {
    if (P1 <= 0.0 || V1 <= 0.0 || T1 <= 0.0 || P2 <= 0.0 || V2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    // T2 = P2V2T1 / P1V1
    return CHEM_OK((P2 * V2 * T1) / (P1 * V1), "K");
}

// ── Hukum Boyle, Charles, Gay-Lussac ─────────────────

ChemResult boyle_P2(chem_float P1, chem_float V1, chem_float V2) {
    if (P1 <= 0.0 || V1 <= 0.0 || V2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((P1 * V1) / V2, "atm");
}

ChemResult boyle_V2(chem_float P1, chem_float V1, chem_float P2) {
    if (P1 <= 0.0 || V1 <= 0.0 || P2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((P1 * V1) / P2, "L");
}

ChemResult charles_V2(chem_float V1, chem_float T1, chem_float T2) {
    if (V1 <= 0.0 || T1 <= 0.0 || T2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((V1 * T2) / T1, "L");
}

ChemResult charles_T2(chem_float V1, chem_float T1, chem_float V2) {
    if (V1 <= 0.0 || T1 <= 0.0 || V2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((V2 * T1) / V1, "K");
}

ChemResult gay_lussac_P2(chem_float P1, chem_float T1, chem_float T2) {
    if (P1 <= 0.0 || T1 <= 0.0 || T2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((P1 * T2) / T1, "atm");
}

ChemResult gay_lussac_T2(chem_float P1, chem_float T1, chem_float P2) {
    if (P1 <= 0.0 || T1 <= 0.0 || P2 <= 0.0)
        return CHEM_ERR("Semua nilai harus lebih dari 0");
    return CHEM_OK((P2 * T1) / P1, "K");
}

// ── Campuran Gas (Dalton) ─────────────────────────────

ChemResult partial_pressure(
    chem_float mol_i,
    chem_float mol_total,
    chem_float P_total
) {
    if (mol_i     <= 0.0) return CHEM_ERR("Mol gas harus lebih dari 0");
    if (mol_total <= 0.0) return CHEM_ERR("Total mol harus lebih dari 0");
    if (P_total   <= 0.0) return CHEM_ERR("Tekanan total harus lebih dari 0");
    if (mol_i > mol_total)
        return CHEM_ERR("Mol gas tidak boleh melebihi total mol");
    // Pi = (ni / ntotal) × Ptotal
    return CHEM_OK((mol_i / mol_total) * P_total, "atm");
}

// ── Konversi Suhu ─────────────────────────────────────

ChemResult celsius_to_kelvin(chem_float celsius) {
    if (celsius < -273.15)
        return CHEM_ERR("Suhu tidak boleh di bawah nol mutlak (-273.15°C)");
    return CHEM_OK(celsius + 273.15, "K");
}

ChemResult kelvin_to_celsius(chem_float kelvin) {
    if (kelvin < 0.0)
        return CHEM_ERR("Suhu Kelvin tidak boleh negatif");
    return CHEM_OK(kelvin - 273.15, "°C");
}