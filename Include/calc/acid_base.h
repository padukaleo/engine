#pragma once
#include "types.h"

// ── pH & pOH dasar ────────────────────────────────────

// pH = -log[H+]
ChemResult calc_pH(chem_float h_plus);

// pOH = -log[OH-]
ChemResult calc_pOH(chem_float oh_minus);

// [H+] dari pH
ChemResult h_plus_from_pH(chem_float pH);

// [OH-] dari pOH
ChemResult oh_minus_from_pOH(chem_float pOH);

// pH + pOH = 14 (25°C)
ChemResult pH_from_pOH(chem_float pOH);
ChemResult pOH_from_pH(chem_float pH);

// ── Asam & Basa ───────────────────────────────────────

// Asam kuat → [H+] = Ca
ChemResult pH_strong_acid(chem_float conc);

// Basa kuat → [OH-] = Cb
ChemResult pH_strong_base(chem_float conc);

// Asam lemah → [H+] = sqrt(Ka × Ca)
ChemResult pH_weak_acid(chem_float Ka, chem_float conc);

// Basa lemah → [OH-] = sqrt(Kb × Cb)
ChemResult pH_weak_base(chem_float Kb, chem_float conc);

// Derajat ionisasi asam lemah → α = [H+] / Ca
ChemResult ionization_degree(chem_float Ka, chem_float conc);

// ── Buffer ────────────────────────────────────────────

// Henderson-Hasselbalch
// pH = pKa + log([A-] / [HA])
ChemResult pH_buffer_acid(
    chem_float Ka,
    chem_float mol_weak_acid,
    chem_float mol_conjugate_base
);

// pOH = pKb + log([B+] / [BOH])
ChemResult pH_buffer_base(
    chem_float Kb,
    chem_float mol_weak_base,
    chem_float mol_conjugate_acid
);

// ── Hidrolisis ────────────────────────────────────────

// Garam dari asam kuat + basa lemah → pH < 7
// [H+] = sqrt(Kw / Kb × C)
ChemResult pH_hydrolysis_acid(chem_float Kb, chem_float conc);

// Garam dari asam lemah + basa kuat → pH > 7
// [OH-] = sqrt(Kw / Ka × C)
ChemResult pH_hydrolysis_base(chem_float Ka, chem_float conc);

// Garam dari asam lemah + basa lemah
// pH = 7 + 0.5 × (pKa - pKb)
ChemResult pH_hydrolysis_both(chem_float Ka, chem_float Kb);