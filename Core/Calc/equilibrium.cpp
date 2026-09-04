#include "equilibrium.h"
#include <math.h>

#define R_LATM  0.08206

// Helper → hitung produk [C]^n dari array
static chem_float power_product(
    const chem_float* values,
    const chem_int*   exponents,
    int               count
) {
    chem_float result = 1.0;
    for (int i = 0; i < count; i++)
        result *= pow(values[i], exponents[i]);
    return result;
}

// ── Tetapan Kesetimbangan ─────────────────────────────

ChemResult calc_Kc(
    const chem_float* cp, const chem_int* ep, int np,
    const chem_float* cr, const chem_int* er, int nr
) {
    if (!cp || !ep || !cr || !er || np <= 0 || nr <= 0)
        return CHEM_ERR("Input tidak valid");

    chem_float num = power_product(cp, ep, np);
    chem_float den = power_product(cr, er, nr);

    if (den == 0.0) return CHEM_ERR("Konsentrasi reaktan tidak boleh 0");
    return CHEM_OK(num / den, "");
}

ChemResult calc_Kp(
    const chem_float* pp, const chem_int* ep, int np,
    const chem_float* pr, const chem_int* er, int nr
) {
    if (!pp || !ep || !pr || !er || np <= 0 || nr <= 0)
        return CHEM_ERR("Input tidak valid");

    chem_float num = power_product(pp, ep, np);
    chem_float den = power_product(pr, er, nr);

    if (den == 0.0) return CHEM_ERR("Tekanan parsial reaktan tidak boleh 0");
    return CHEM_OK(num / den, "");
}

ChemResult Kp_from_Kc(chem_float Kc, chem_float temp_K, chem_int delta_n) {
    if (Kc     <= 0.0) return CHEM_ERR("Kc harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // Kp = Kc × (RT)^Δn
    return CHEM_OK(Kc * pow(R_LATM * temp_K, delta_n), "");
}

ChemResult Kc_from_Kp(chem_float Kp, chem_float temp_K, chem_int delta_n) {
    if (Kp     <= 0.0) return CHEM_ERR("Kp harus lebih dari 0");
    if (temp_K <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    // Kc = Kp / (RT)^Δn
    chem_float rt = pow(R_LATM * temp_K, delta_n);
    if (rt == 0.0) return CHEM_ERR("RT^Δn tidak boleh 0");
    return CHEM_OK(Kp / rt, "");
}

// ── Reaction Quotient ─────────────────────────────────

ChemResult calc_Qc(
    const chem_float* cp, const chem_int* ep, int np,
    const chem_float* cr, const chem_int* er, int nr
) {
    // Sama persis dengan Kc, bedanya ini kondisi tidak setimbang
    return calc_Kc(cp, ep, np, cr, er, nr);
}

int reaction_direction(chem_float Q, chem_float K) {
    chem_float diff = Q - K;
    if (diff < -1e-12) return  1;   // maju
    if (diff >  1e-12) return -1;   // balik
    return 0;                        // setimbang
}

// ── Derajat Ionisasi ──────────────────────────────────

ChemResult ionization_degree_eq(chem_float Ka, chem_float conc) {
    if (Ka   <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // α = sqrt(Ka / C)
    chem_float alpha = sqrt(Ka / conc);
    if (alpha > 1.0) alpha = 1.0; // tidak bisa lebih dari 100%
    return CHEM_OK(alpha * 100.0, "%");
}

// ── Kesetimbangan Heterogen ───────────────────────────

int precipitation_check(chem_float Q, chem_float Ksp) {
    chem_float diff = Q - Ksp;
    if (diff >  1e-20) return  1;   // endapan terbentuk
    if (diff < -1e-20) return -1;   // belum jenuh
    return 0;                        // tepat jenuh
}