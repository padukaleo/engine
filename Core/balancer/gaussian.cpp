#include "gaussian.h"
#include <math.h>
#include <string.h>

#define EPSILON 1e-10

// GCD dua bilangan bulat
static int gcd(int a, int b) {
    a = abs(a);
    b = abs(b);
    while (b) { int t = b; b = a % b; a = t; }
    return a;
}

// LCM dua bilangan bulat
static int lcm(int a, int b) {
    return (a / gcd(a, b)) * b;
}

void gaussian_rref(Matrix* m) {
    int pivot_row = 0;

    for (int col = 0; col < m->cols && pivot_row < m->rows; col++) {
        // Cari pivot di kolom ini
        int pivot = -1;
        for (int row = pivot_row; row < m->rows; row++) {
            if (fabs(m->data[row][col]) > EPSILON) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) continue;

        // Swap baris pivot ke posisi pivot_row
        if (pivot != pivot_row) {
            for (int j = 0; j <= m->cols; j++) {
                chem_float tmp         = m->data[pivot_row][j];
                m->data[pivot_row][j]  = m->data[pivot][j];
                m->data[pivot][j]      = tmp;
            }
        }

        // Normalisasi baris pivot
        chem_float scale = m->data[pivot_row][col];
        for (int j = 0; j <= m->cols; j++)
            m->data[pivot_row][j] /= scale;

        // Eliminasi kolom ini di semua baris lain
        for (int row = 0; row < m->rows; row++) {
            if (row == pivot_row) continue;
            chem_float factor = m->data[row][col];
            for (int j = 0; j <= m->cols; j++)
                m->data[row][j] -= factor * m->data[pivot_row][j];
        }

        pivot_row++;
    }
}

int extract_solution(
    const Matrix* m,
    int*          coefficients,
    int           n_species
) {
    // Set variabel bebas terakhir = 1, lalu solve
    // Semua koefisien harus positif dan bulat
    chem_float sol[MAX_SPECIES] = {0};
    sol[n_species - 1] = 1.0;  // variabel bebas

    // Back-substitution dari RREF
    for (int row = m->rows - 1; row >= 0; row--) {
        // Cari pivot di baris ini
        int pivot_col = -1;
        for (int col = 0; col < n_species; col++) {
            if (fabs(m->data[row][col] - 1.0) < EPSILON) {
                pivot_col = col;
                break;
            }
        }
        if (pivot_col < 0) continue;

        sol[pivot_col] = 0.0;
        for (int col = pivot_col + 1; col < n_species; col++)
            sol[pivot_col] -= m->data[row][col] * sol[col];
    }

    // Konversi ke integer dengan LCM
    // Cari denominator terbesar (pendekatan sederhana dengan skala)
    int scale = 1;
    for (int i = 0; i < n_species; i++) {
        // Kalikan sampai semua jadi integer (max 12 iterasi)
        for (int k = 1; k <= 12; k++) {
            chem_float val = sol[i] * scale;
            if (fabs(val - round(val)) < 1e-6) break;
            scale = lcm(scale, k + 1);
        }
    }

    // Isi output dan validasi semua positif
    int g = 0;
    for (int i = 0; i < n_species; i++) {
        coefficients[i] = (int)round(sol[i] * scale);
        if (coefficients[i] <= 0) return 0; // gagal
        g = (g == 0) ? coefficients[i] : gcd(g, coefficients[i]);
    }

    // Sederhanakan dengan GCD
    if (g > 1)
        for (int i = 0; i < n_species; i++)
            coefficients[i] /= g;

    return 1;
}