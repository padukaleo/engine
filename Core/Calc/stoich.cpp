#include "stoich.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#define AVOGADRO    6.02214076e23
#define MOLAR_STP   22.4        // L/mol di 0°C
#define MOLAR_SATP  24.8        // L/mol di 25°C

// ── Konversi dasar ────────────────────────────────────

ChemResult mol_to_particles(chem_float mol) {
    if (mol <= 0.0)
        return CHEM_ERR("Mol harus lebih dari 0");
    return CHEM_OK(mol * AVOGADRO, "partikel");
}

ChemResult particles_to_mol(chem_float particles) {
    if (particles <= 0.0)
        return CHEM_ERR("Jumlah partikel harus lebih dari 0");
    return CHEM_OK(particles / AVOGADRO, "mol");
}

ChemResult mol_to_volume_stp(chem_float mol) {
    if (mol <= 0.0)
        return CHEM_ERR("Mol harus lebih dari 0");
    return CHEM_OK(mol * MOLAR_STP, "L");
}

ChemResult mol_to_volume_satp(chem_float mol) {
    if (mol <= 0.0)
        return CHEM_ERR("Mol harus lebih dari 0");
    return CHEM_OK(mol * MOLAR_SATP, "L");
}

ChemResult volume_to_mol_stp(chem_float liter) {
    if (liter <= 0.0)
        return CHEM_ERR("Volume harus lebih dari 0");
    return CHEM_OK(liter / MOLAR_STP, "mol");
}

ChemResult volume_to_mol_satp(chem_float liter) {
    if (liter <= 0.0)
        return CHEM_ERR("Volume harus lebih dari 0");
    return CHEM_OK(liter / MOLAR_SATP, "mol");
}

// ── Stoikiometri reaksi ───────────────────────────────

ChemResult mol_by_ratio(
    chem_float mol_known,
    chem_int   coef_known,
    chem_int   coef_target
) {
    if (coef_known <= 0 || coef_target <= 0)
        return CHEM_ERR("Koefisien harus lebih dari 0");
    if (mol_known <= 0.0)
        return CHEM_ERR("Mol harus lebih dari 0");

    return CHEM_OK(mol_known * coef_target / coef_known, "mol");
}

int limiting_reagent(
    const chem_float* mols,
    const chem_int*   coefs,
    int               count
) {
    if (!mols || !coefs || count <= 0) return -1;

    int    limiting = 0;
    double min_ratio = mols[0] / coefs[0];

    for (int i = 1; i < count; i++) {
        double ratio = mols[i] / coefs[i];
        if (ratio < min_ratio) {
            min_ratio = ratio;
            limiting  = i;
        }
    }

    return limiting;
}

ChemResult percent_yield(
    chem_float actual_gram,
    chem_float theoretical_gram
) {
    if (actual_gram <= 0.0)
        return CHEM_ERR("Massa aktual harus lebih dari 0");
    if (theoretical_gram <= 0.0)
        return CHEM_ERR("Massa teoritis harus lebih dari 0");
    if (actual_gram > theoretical_gram)
        return CHEM_ERR("Massa aktual tidak boleh melebihi teoritis");

    return CHEM_OK((actual_gram / theoretical_gram) * 100.0, "%");
}

// ── Komposisi ─────────────────────────────────────────

ChemResult percent_mass(const char* formula, const char* element_symbol) {
    if (!formula || !element_symbol)
        return CHEM_ERR("Input tidak valid");

    Compound c = parse_formula(formula);
    if (c.size == 0)
        return CHEM_ERR("Formula tidak valid");

    ChemResult mr = molar_mass_from_compound(&c);
    if (!mr.success) return mr;

    // Cari elemen yang diminta
    for (int i = 0; i < c.size; i++) {
        if (strcmp(c.atoms[i].symbol, element_symbol) == 0) {
            chem_float elem_mass = db_get_atomic_mass(element_symbol)
                                 * c.atoms[i].count;
            return CHEM_OK((elem_mass / mr.value) * 100.0, "%");
        }
    }

    return CHEM_ERR("Elemen tidak ditemukan dalam senyawa");
}

ChemResult empirical_formula(
    const char**      symbols,
    const chem_float* percents,
    int               count,
    char*             out_formula,
    int               out_size
) {
    if (!symbols || !percents || count <= 0)
        return CHEM_ERR("Input tidak valid");

    // 1. Bagi persen dengan massa atom → mol ratio
    chem_float mols[20] = {0};
    for (int i = 0; i < count; i++) {
        chem_float mass = db_get_atomic_mass(symbols[i]);
        if (mass <= 0.0)
            return CHEM_ERR("Simbol elemen tidak valid");
        mols[i] = percents[i] / mass;
    }

    // 2. Bagi semua dengan nilai terkecil
    chem_float min_mol = mols[0];
    for (int i = 1; i < count; i++)
        if (mols[i] < min_mol) min_mol = mols[i];

    chem_float ratios[20] = {0};
    for (int i = 0; i < count; i++)
        ratios[i] = mols[i] / min_mol;

    // 3. Bulatkan ke integer terdekat (toleransi 0.1)
    chem_int int_ratios[20] = {0};
    for (int i = 0; i < count; i++) {
        // Cek apakah perlu kalikan (contoh: 1.5 → kalikan 2)
        chem_float r = ratios[i];
        if (fabs(r - round(r)) < 0.1) {
            int_ratios[i] = (chem_int)round(r);
        } else if (fabs(r * 2 - round(r * 2)) < 0.1) {
            // Semua × 2
            for (int j = 0; j < count; j++)
                int_ratios[j] = (chem_int)round(ratios[j] * 2);
            break;
        } else if (fabs(r * 3 - round(r * 3)) < 0.1) {
            // Semua × 3
            for (int j = 0; j < count; j++)
                int_ratios[j] = (chem_int)round(ratios[j] * 3);
            break;
        } else {
            int_ratios[i] = (chem_int)round(r);
        }
    }

    // 4. Build string output "C6H12O6" dll
    int pos = 0;
    for (int i = 0; i < count; i++) {
        int written = 0;
        if (int_ratios[i] == 1)
            written = snprintf(out_formula + pos, out_size - pos,
                               "%s", symbols[i]);
        else
            written = snprintf(out_formula + pos, out_size - pos,
                               "%s%d", symbols[i], int_ratios[i]);

        if (written < 0 || pos + written >= out_size)
            return CHEM_ERR("Buffer output terlalu kecil");
        pos += written;
    }

    return CHEM_OK(1.0, "empirical_formula");
}