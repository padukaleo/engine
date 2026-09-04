#include "nuclear.h"
#include <math.h>

#define LN2             0.693147180559945
#define U_TO_MEV        931.494          // 1 u = 931.494 MeV/c²

ChemResult decay_constant(chem_float half_life) {
    if (half_life <= 0.0) return CHEM_ERR("Waktu paruh harus lebih dari 0");
    return CHEM_OK(LN2 / half_life, "1/s");
}

ChemResult remaining_amount(
    chem_float initial,
    chem_float half_life,
    chem_float time
) {
    if (initial   <= 0.0) return CHEM_ERR("Jumlah awal harus lebih dari 0");
    if (half_life <= 0.0) return CHEM_ERR("Waktu paruh harus lebih dari 0");
    if (time      <= 0.0) return CHEM_ERR("Waktu harus lebih dari 0");
    return CHEM_OK(initial * pow(0.5, time / half_life), "");
}

ChemResult activity(chem_float lambda, chem_float N) {
    if (lambda <= 0.0) return CHEM_ERR("Konstanta peluruhan harus lebih dari 0");
    if (N      <= 0.0) return CHEM_ERR("Jumlah inti harus lebih dari 0");
    return CHEM_OK(lambda * N, "Bq");
}

ChemResult time_elapsed(
    chem_float initial,
    chem_float remaining,
    chem_float half_life
) {
    if (initial   <= 0.0) return CHEM_ERR("Jumlah awal harus lebih dari 0");
    if (remaining <= 0.0) return CHEM_ERR("Jumlah sisa harus lebih dari 0");
    if (half_life <= 0.0) return CHEM_ERR("Waktu paruh harus lebih dari 0");
    if (remaining > initial)
        return CHEM_ERR("Sisa tidak boleh melebihi jumlah awal");
    // t = t½ × ln(N0/N) / ln2
    return CHEM_OK(half_life * log(initial / remaining) / LN2, "s");
}

ChemResult percent_remaining(
    chem_float initial,
    chem_float half_life,
    chem_float time
) {
    ChemResult r = remaining_amount(initial, half_life, time);
    if (!r.success) return r;
    return CHEM_OK((r.value / initial) * 100.0, "%");
}

ChemResult binding_energy(chem_float mass_defect_u) {
    if (mass_defect_u <= 0.0)
        return CHEM_ERR("Defek massa harus lebih dari 0");
    return CHEM_OK(mass_defect_u * U_TO_MEV, "MeV");
}

ChemResult binding_energy_per_nucleon(
    chem_float mass_defect_u,
    chem_int   mass_number
) {
    if (mass_defect_u <= 0.0) return CHEM_ERR("Defek massa harus lebih dari 0");
    if (mass_number   <= 0)   return CHEM_ERR("Nomor massa harus lebih dari 0");
    chem_float total = mass_defect_u * U_TO_MEV;
    return CHEM_OK(total / mass_number, "MeV/nukleon");
}