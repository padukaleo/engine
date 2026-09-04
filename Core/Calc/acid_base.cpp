#include "acid_base.h"
#include <math.h>

#define KW      1e-14       // tetapan air 25°C
#define PKW     14.0

// ── pH & pOH dasar ────────────────────────────────────

ChemResult calc_pH(chem_float h_plus) {
    if (h_plus <= 0.0) return CHEM_ERR("[H+] harus lebih dari 0");
    return CHEM_OK(-log10(h_plus), "");
}

ChemResult calc_pOH(chem_float oh_minus) {
    if (oh_minus <= 0.0) return CHEM_ERR("[OH-] harus lebih dari 0");
    return CHEM_OK(-log10(oh_minus), "");
}

ChemResult h_plus_from_pH(chem_float pH) {
    if (pH < 0.0 || pH > 14.0) return CHEM_ERR("pH harus antara 0-14");
    return CHEM_OK(pow(10.0, -pH), "mol/L");
}

ChemResult oh_minus_from_pOH(chem_float pOH) {
    if (pOH < 0.0 || pOH > 14.0) return CHEM_ERR("pOH harus antara 0-14");
    return CHEM_OK(pow(10.0, -pOH), "mol/L");
}

ChemResult pH_from_pOH(chem_float pOH) {
    if (pOH < 0.0 || pOH > 14.0) return CHEM_ERR("pOH harus antara 0-14");
    return CHEM_OK(PKW - pOH, "");
}

ChemResult pOH_from_pH(chem_float pH) {
    if (pH < 0.0 || pH > 14.0) return CHEM_ERR("pH harus antara 0-14");
    return CHEM_OK(PKW - pH, "");
}

// ── Asam & Basa ───────────────────────────────────────

ChemResult pH_strong_acid(chem_float conc) {
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    return CHEM_OK(-log10(conc), "");
}

ChemResult pH_strong_base(chem_float conc) {
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    chem_float pOH = -log10(conc);
    return CHEM_OK(PKW - pOH, "");
}

ChemResult pH_weak_acid(chem_float Ka, chem_float conc) {
    if (Ka   <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // [H+] = sqrt(Ka × C)
    chem_float h_plus = sqrt(Ka * conc);
    return CHEM_OK(-log10(h_plus), "");
}

ChemResult pH_weak_base(chem_float Kb, chem_float conc) {
    if (Kb   <= 0.0) return CHEM_ERR("Kb harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // [OH-] = sqrt(Kb × C)
    chem_float oh_minus = sqrt(Kb * conc);
    chem_float pOH      = -log10(oh_minus);
    return CHEM_OK(PKW - pOH, "");
}

ChemResult ionization_degree(chem_float Ka, chem_float conc) {
    if (Ka   <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // α = sqrt(Ka / C)
    return CHEM_OK(sqrt(Ka / conc) * 100.0, "%");
}

// ── Buffer ────────────────────────────────────────────

ChemResult pH_buffer_acid(
    chem_float Ka,
    chem_float mol_acid,
    chem_float mol_base
) {
    if (Ka       <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (mol_acid <= 0.0) return CHEM_ERR("Mol asam harus lebih dari 0");
    if (mol_base <= 0.0) return CHEM_ERR("Mol basa konjugat harus lebih dari 0");
    // pH = pKa + log([A-]/[HA])
    chem_float pKa = -log10(Ka);
    chem_float pH  = pKa + log10(mol_base / mol_acid);
    return CHEM_OK(pH, "");
}

ChemResult pH_buffer_base(
    chem_float Kb,
    chem_float mol_base,
    chem_float mol_acid
) {
    if (Kb       <= 0.0) return CHEM_ERR("Kb harus lebih dari 0");
    if (mol_base <= 0.0) return CHEM_ERR("Mol basa harus lebih dari 0");
    if (mol_acid <= 0.0) return CHEM_ERR("Mol asam konjugat harus lebih dari 0");
    // pOH = pKb + log([B+]/[BOH])
    chem_float pKb = -log10(Kb);
    chem_float pOH = pKb + log10(mol_acid / mol_base);
    return CHEM_OK(PKW - pOH, "");
}

// ── Hidrolisis ────────────────────────────────────────

ChemResult pH_hydrolysis_acid(chem_float Kb, chem_float conc) {
    if (Kb   <= 0.0) return CHEM_ERR("Kb harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // Kh = Kw / Kb
    // [H+] = sqrt(Kh × C)
    chem_float Kh     = KW / Kb;
    chem_float h_plus = sqrt(Kh * conc);
    return CHEM_OK(-log10(h_plus), "");
}

ChemResult pH_hydrolysis_base(chem_float Ka, chem_float conc) {
    if (Ka   <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (conc <= 0.0) return CHEM_ERR("Konsentrasi harus lebih dari 0");
    // Kh = Kw / Ka
    // [OH-] = sqrt(Kh × C)
    chem_float Kh      = KW / Ka;
    chem_float oh_minus = sqrt(Kh * conc);
    chem_float pOH     = -log10(oh_minus);
    return CHEM_OK(PKW - pOH, "");
}

ChemResult pH_hydrolysis_both(chem_float Ka, chem_float Kb) {
    if (Ka <= 0.0) return CHEM_ERR("Ka harus lebih dari 0");
    if (Kb <= 0.0) return CHEM_ERR("Kb harus lebih dari 0");
    // pH = 7 + 0.5 × (pKa - pKb)
    chem_float pKa = -log10(Ka);
    chem_float pKb = -log10(Kb);
    return CHEM_OK(7.0 + 0.5 * (pKa - pKb), "");
}