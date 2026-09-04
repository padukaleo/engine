#include "matrix.h"
#include <string.h>

// Cari index elemen di daftar, return -1 kalau tidak ada
static int find_element(
    const char  elements[][4],
    int         count,
    const char* symbol
) {
    for (int i = 0; i < count; i++)
        if (strcmp(elements[i], symbol) == 0) return i;
    return -1;
}

int build_matrix(
    const char** reactants, int n_reactants,
    const char** products,  int n_products,
    Matrix*      out
) {
    memset(out, 0, sizeof(Matrix));

    int   n_species  = n_reactants + n_products;
    char  elements[MAX_ELEMENTS][4];
    int   n_elements = 0;

    // Parse semua senyawa, kumpulkan elemen unik
    Compound compounds[MAX_SPECIES];
    for (int i = 0; i < n_reactants; i++) {
        compounds[i] = parse_formula(reactants[i]);
        for (int j = 0; j < compounds[i].size; j++) {
            const char* sym = compounds[i].atoms[j].symbol;
            if (find_element(elements, n_elements, sym) < 0) {
                if (n_elements >= MAX_ELEMENTS) return 0;
                strncpy(elements[n_elements++], sym, 3);
            }
        }
    }
    for (int i = 0; i < n_products; i++) {
        int idx = n_reactants + i;
        compounds[idx] = parse_formula(products[i]);
        for (int j = 0; j < compounds[idx].size; j++) {
            const char* sym = compounds[idx].atoms[j].symbol;
            if (find_element(elements, n_elements, sym) < 0) {
                if (n_elements >= MAX_ELEMENTS) return 0;
                strncpy(elements[n_elements++], sym, 3);
            }
        }
    }

    out->rows = n_elements;
    out->cols = n_species;

    // Isi matrix
    // reaktan → positif, produk → negatif
    for (int i = 0; i < n_reactants; i++) {
        for (int j = 0; j < compounds[i].size; j++) {
            int row = find_element(elements, n_elements,
                                   compounds[i].atoms[j].symbol);
            if (row >= 0)
                out->data[row][i] += compounds[i].atoms[j].count;
        }
    }
    for (int i = 0; i < n_products; i++) {
        int idx = n_reactants + i;
        for (int j = 0; j < compounds[idx].size; j++) {
            int row = find_element(elements, n_elements,
                                   compounds[idx].atoms[j].symbol);
            if (row >= 0)
                out->data[row][n_reactants + i] -=
                    compounds[idx].atoms[j].count;
        }
    }

    return 1;
}