#include "colligative.h"
#include <math.h>

#define R_LATM  0.08206

ChemResult boiling_point_elevation(
    chem_float Kb,
    chem_float molality,
    chem_float i
) {
    if (Kb       <= 0.0) return CHEM_ERR("Kb harus lebih dari 0");
    if (molality <= 0.0) return CHEM_ERR("Molalitas harus lebih dari 0");
    if (i        <= 0.0) return CHEM_ERR("Faktor Van't Hoff harus lebih dari 0");
    return CHEM_OK(i * Kb * molality, "°C");
}

ChemResult freezing_point_depression(
    chem_float Kf,
    chem_float molality,
    chem_float i
) {
    if (Kf       <= 0.0) return CHEM_ERR("Kf harus lebih dari 0");
    if (molality <= 0.0) return CHEM_ERR("Molalitas harus lebih dari 0");
    if (i        <= 0.0) return CHEM_ERR("Faktor Van't Hoff harus lebih dari 0");
    return CHEM_OK(i * Kf * molality, "°C");
}

ChemResult osmotic_pressure(
    chem_float molarity,
    chem_float temp_K,
    chem_float i
) {
    if (molarity <= 0.0) return CHEM_ERR("Molaritas harus lebih dari 0");
    if (temp_K   <= 0.0) return CHEM_ERR("Suhu harus lebih dari 0 Kelvin");
    if (i        <= 0.0) return CHEM_ERR("Faktor Van't Hoff harus lebih dari 0");
    // π = iMRT
    return CHEM_OK(i * molarity * R_LATM * temp_K, "atm");
}

ChemResult vapor_pressure_lowering(
    chem_float P_pure,
    chem_float mole_frac_solute
) {
    if (P_pure          <= 0.0) return CHEM_ERR("Tekanan uap murni harus lebih dari 0");
    if (mole_frac_solute < 0.0 || mole_frac_solute > 1.0)
        return CHEM_ERR("Fraksi mol harus antara 0-1");
    return CHEM_OK(mole_frac_solute * P_pure, "mmHg");
}

ChemResult vapor_pressure_solution(
    chem_float P_pure,
    chem_float mole_frac_solvent
) {
    if (P_pure            <= 0.0) return CHEM_ERR("Tekanan uap murni harus lebih dari 0");
    if (mole_frac_solvent < 0.0 || mole_frac_solvent > 1.0)
        return CHEM_ERR("Fraksi mol harus antara 0-1");
    // P = Xsolvent × P°
    return CHEM_OK(mole_frac_solvent * P_pure, "mmHg");
}

ChemResult vant_hoff_factor(chem_float alpha, chem_int n_ions) {
    if (alpha < 0.0 || alpha > 1.0)
        return CHEM_ERR("Derajat ionisasi harus antara 0-1");
    if (n_ions <= 1)
        return CHEM_ERR("Jumlah ion harus lebih dari 1");
    // i = 1 + α(n-1)
    return CHEM_OK(1.0 + alpha * (n_ions - 1), "");
}

ChemResult molality_from_deltaT(
    chem_float deltaT,
    chem_float K,
    chem_float i
) {
    if (deltaT <= 0.0) return CHEM_ERR("ΔT harus lebih dari 0");
    if (K      <= 0.0) return CHEM_ERR("Tetapan K harus lebih dari 0");
    if (i      <= 0.0) return CHEM_ERR("Faktor Van't Hoff harus lebih dari 0");
    // m = ΔT / (i × K)
    return CHEM_OK(deltaT / (i * K), "mol/kg");
}

ChemResult molar_mass_from_freezing(
    chem_float Kf,
    chem_float mass_solute_g,
    chem_float mass_solvent_g,
    chem_float delta_Tf
) {
    if (Kf            <= 0.0) return CHEM_ERR("Kf harus lebih dari 0");
    if (mass_solute_g <= 0.0) return CHEM_ERR("Massa solute harus lebih dari 0");
    if (mass_solvent_g<= 0.0) return CHEM_ERR("Massa solvent harus lebih dari 0");
    if (delta_Tf      <= 0.0) return CHEM_ERR("ΔTf harus lebih dari 0");
    // Mr = (Kf × w × 1000) / (ΔTf × W)
    chem_float Mr = (Kf * mass_solute_g * 1000.0) / (delta_Tf * mass_solvent_g);
    return CHEM_OK(Mr, "g/mol");
}