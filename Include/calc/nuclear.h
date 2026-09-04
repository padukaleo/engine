#pragma once
#include "types.h"

// ── Peluruhan Radioaktif ──────────────────────────────

// λ = ln2 / t½
ChemResult decay_constant(chem_float half_life);

// N = N0 × (1/2)^(t/t½)
ChemResult remaining_amount(
    chem_float initial,
    chem_float half_life,
    chem_float time
);

// Aktivitas → A = λN
ChemResult activity(chem_float lambda, chem_float N);

// Waktu yang sudah berlalu
// t = t½ × log(N0/N) / log2
ChemResult time_elapsed(
    chem_float initial,
    chem_float remaining,
    chem_float half_life
);

// Persentase yang tersisa
ChemResult percent_remaining(
    chem_float initial,
    chem_float half_life,
    chem_float time
);

// ── Energi Nuklir ─────────────────────────────────────

// E = mc² → dalam MeV
// defect_mass dalam satuan u (unified atomic mass unit)
ChemResult binding_energy(chem_float mass_defect_u);

// Energi ikat per nukleon
ChemResult binding_energy_per_nucleon(
    chem_float mass_defect_u,
    chem_int   mass_number      // A
);