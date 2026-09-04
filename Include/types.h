#pragma once
#include <stdint.h>

//Tipe dasar, biar konsisten
typedef double chem_float;
typedef int32_t chem_int;
typedef uint8_t chem_u8;

// hasil kalkulasi
typedef struct {
    chem_float value;
    int succsess; // 1 = ok, 0 = gagal
    const char*error; // error message
    const char*unit; // "g/mol", mol, dan lain lain
} ChemResult;

typedef struct {
    chem_u8 number; // nomor atom
    const char* symbol; // "H", "Fe", dll
    const char* name; // "Hydrogen"
    chem_float atomic_mass; // gram/mol
    chem_u8 valence; // elektronvalensi
    chem_u8 period;
    chem_u8 group;
} Element;

// compound hasil pairsing
typedef struct {
    const char* symbol;
    chem_int count;
} AtomCount;

typedef struct {
    AtomCount atoms[20]; // max 20 atom
    int size;
    chem_float molar_mass; 
} Compound;

//hasil balancing
typedef struct {
    int coefficients[20]; // max 20 senyawa
    int size;
    int succsess; 
    const char* error;
} BalanceResult;

// Macro helper
#define CHEM_OK(val, u) (ChemResult){ val, 1, NULL, u}
#define CHEM_ERR(msg) (ChemResult){0.0, 0, msg, NULL}