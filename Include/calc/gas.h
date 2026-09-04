#pragma once
#include "types.h"

#define R_IDEAL     8.314       // J/mol·K
#define R_LATM      0.08206     // L·atm/mol·K
#define STP_TEMP    273.15      // K
#define STP_PRESS   1.0         // atm
#define SATP_TEMP   298.15      // K

// ── Hukum Gas Ideal ───────────────────────────────────

// PV = nRT → cari P
ChemResult ideal_gas_pressure(
    chem_float n,       // mol
    chem_float T,       // Kelvin
    chem_float V        // Liter
);

// PV = nRT → cari V
ChemResult ideal_gas_volume(
    chem_float n,
    chem_float T,
    chem_float P        // atm
);

// PV = nRT → cari n
ChemResult ideal_gas_mol(
    chem_float P,
    chem_float V,
    chem_float T
);

// PV = nRT → cari T
ChemResult ideal_gas_temp(
    chem_float P,
    chem_float V,
    chem_float n
);

// ── Hukum Gas Gabungan ────────────────────────────────

// P1V1/T1 = P2V2/T2 → cari P2
ChemResult combined_gas_P2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float V2, chem_float T2
);

// → cari V2
ChemResult combined_gas_V2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float P2, chem_float T2
);

// → cari T2
ChemResult combined_gas_T2(
    chem_float P1, chem_float V1, chem_float T1,
    chem_float P2, chem_float V2
);

// ── Hukum Boyle, Charles, Gay-Lussac ─────────────────

// Boyle → P1V1 = P2V2 → cari P2
ChemResult boyle_P2(chem_float P1, chem_float V1, chem_float V2);

// Boyle → cari V2
ChemResult boyle_V2(chem_float P1, chem_float V1, chem_float P2);

// Charles → V1/T1 = V2/T2 → cari V2
ChemResult charles_V2(chem_float V1, chem_float T1, chem_float T2);

// Charles → cari T2
ChemResult charles_T2(chem_float V1, chem_float T1, chem_float V2);

// Gay-Lussac → P1/T1 = P2/T2 → cari P2
ChemResult gay_lussac_P2(chem_float P1, chem_float T1, chem_float T2);

// Gay-Lussac → cari T2
ChemResult gay_lussac_T2(chem_float P1, chem_float T1, chem_float P2);

// ── Campuran Gas (Dalton) ─────────────────────────────

// Tekanan parsial → Pi = Xi × Ptotal
ChemResult partial_pressure(
    chem_float mol_i,       // mol gas ke-i
    chem_float mol_total,   // total mol semua gas
    chem_float P_total      // tekanan total
);

// ── Konversi Suhu ─────────────────────────────────────

ChemResult celsius_to_kelvin(chem_float celsius);
ChemResult kelvin_to_celsius(chem_float kelvin);