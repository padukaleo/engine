#include "balancer.h"
#include <string.h>
#include <stdlib.h>

#define MAX_FORMULA_LEN 32

// Split string "A + B + C" → array ["A","B","C"]
static int split_species(
    const char* str,
    char        out[][MAX_FORMULA_LEN],
    int         max_count
) {
    int count = 0;
    int i = 0, len = (int)strlen(str);

    while (i < len && count < max_count) {
        // Skip spasi
        while (i < len && str[i] == ' ') i++;
        if (i >= len) break;

        // Skip "+" separator
        if (str[i] == '+') { i++; continue; }

        // Baca formula sampai spasi atau +
        int j = 0;
        while (i < len && str[i] != ' ' && str[i] != '+') {
            if (j < MAX_FORMULA_LEN - 1)
                out[count][j++] = str[i];
            i++;
        }
        out[count][j] = '\0';
        if (j > 0) count++;
    }
    return count;
}

BalancerResult balance_equation(const char* equation) {
    BalancerResult result;
    memset(&result, 0, sizeof(BalancerResult));

    if (!equation || strlen(equation) == 0) {
        result.error = "Persamaan kosong";
        return result;
    }

    // Cari separator "->" atau "→" atau "="
    const char* sep = strstr(equation, "->");
    if (!sep) sep = strstr(equation, "=");
    if (!sep) {
        result.error = "Separator reaksi tidak ditemukan (gunakan -> atau =)";
        return result;
    }

    int sep_len = (strncmp(sep, "->", 2) == 0) ? 2 : 1;

    // Split reaktan dan produk
    char lhs[128] = {0};
    char rhs[128] = {0};
    strncpy(lhs, equation, sep - equation);
    strncpy(rhs, sep + sep_len, sizeof(rhs) - 1);

    char reactant_arr[MAX_SPECIES][MAX_FORMULA_LEN];
    char product_arr[MAX_SPECIES][MAX_FORMULA_LEN];

    int n_r = split_species(lhs, reactant_arr, MAX_SPECIES);
    int n_p = split_species(rhs, product_arr,  MAX_SPECIES);

    if (n_r == 0 || n_p == 0) {
        result.error = "Reaktan atau produk tidak valid";
        return result;
    }

    // Buat array pointer untuk diteruskan
    const char* r_ptrs[MAX_SPECIES];
    const char* p_ptrs[MAX_SPECIES];
    for (int i = 0; i < n_r; i++) r_ptrs[i] = reactant_arr[i];
    for (int i = 0; i < n_p; i++) p_ptrs[i] = product_arr[i];

    return balance_equation_parts(r_ptrs, n_r, p_ptrs, n_p);
}

BalancerResult balance_equation_parts(
    const char** reactants, int n_reactants,
    const char** products,  int n_products
) {
    BalancerResult result;
    memset(&result, 0, sizeof(BalancerResult));

    int n_species = n_reactants + n_products;
    if (n_species < 2 || n_species > MAX_SPECIES) {
        result.error = "Jumlah senyawa tidak valid";
        return result;
    }

    // Build matrix
    Matrix m;
    if (!build_matrix(reactants, n_reactants, products, n_products, &m)) {
        result.error = "Gagal membangun matrix";
        return result;
    }

    // RREF
    gaussian_rref(&m);

    // Ekstrak solusi
    int coefficients[MAX_SPECIES] = {0};
    if (!extract_solution(&m, coefficients, n_species)) {
        result.error = "Persamaan tidak dapat diseimbangkan";
        return result;
    }

    // Isi result
    result.size = n_species;
    for (int i = 0; i < n_species; i++)
        result.coefficients[i] = coefficients[i];
    result.success = 1;

    return result;
}