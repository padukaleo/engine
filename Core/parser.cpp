#include "parser.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

int is_upper(char c) {
    return c >= 'A' && c <= 'Z';
}

int is_lower(char c) {
    return c >= 'a' && c <= 'z';
}

int is_digit(char c) {
    return c >= '0' && c <= '9';
}

static void add_atom(Compound* c, const char* symbol, chem_int count) {
    // Check if the atom already exists in the compound
    for (int i = 0; i < c->size; i++) {
        if (strcmp(c->atoms[i].symbol, symbol) == 0) {
            c->atoms[i].count += count;
            return;
        }
    }
    // If the atom doesn't exist, add it
    if (c->size < 20) {
        c->atoms[c->size].symbol = symbol;
        c->atoms[c->size].count = count;
        c->size++;
    }
}

// core parser

//multiplier

static int parse_segment(
    const char*f, //formula str
    int i, //index sekarang
    int len, //panjan str
    int mult, //multiplier bracket
    Compound* out //output
) {
    while ( i < len) {
        char c = f[i];

        //buka braket
        if ( c == '(' || c == '[') {
            char close = ( c == '(') ? ')' : ']';
            Compound inner;
            memset(&inner, 0, sizeof(Compound));

            //parse
            i = parse_segment(f, i + 1, len, 1, &inner);

            //pastikan isi braket tertutup
            if (i >= len || f[i] != close) {
                return -1;
            }
            i++; //skip close bracket

            //baca angka setelah braket
            int num = 0; 
            while (i < len && is_digit(f[i])) {
                num = num * 10 + (f[i] - '0');
                i++;
            }
            if ( num == 0) num = 1;

            //gabungkan bracket, num, mult ke input
            for (int j = 0; j < inner.size; j++) {
                add_atom(out, inner.atoms[j].symbol, inner.atoms[j].count * num * mult);
            }
        } continue;
    }
}

// tutup bracket
if (c == ')' || c == ']') {
    break;
}

// huruf kapital
if (is_upper(c)) {
    char symbol[3] = {0};
    symbol[0] = f[i++];
    if (i < len && is_lower(f[i])) {
        symbol[1] = f[i++];
    }
    //baca angka setelah simbol
    int num = 0;
    while (i < len && is_digit(f[i])) {
        num = num * 10 + (f[i] - '0');
        i++;
    }
    if (num == 0) num = 1;
    add_atom(out, symbol, num * mult);
    continue;
}
return i;

// public api
int validate_formula(const char* formula) {
    if (!formula || strlen(formula) == 0) return 0;

    int depth = 0;
    for (int i = 0; formula[i]; i++) {
        char c = formula[i];
        if (c == '(' || c == '[') {
            depth++;
        }
        else if (c == ')' || c == ']') {
            depth--;
            if (depth < 0) return 0; 
    }
    else if (!is_upper(c) && !is_lower(c) && !is_digit(c)) {
        return 0;
    }
  }
  return depth == 0;
}

Compound parse_formula(const char* formula) {
    Compound result;
    memset(&result, 0, sizeof(Compound));

    if (!validate_formula(formula)) {
        return result; // return empty compound on invalid formula
    }

    int len = (int)strlen(formula);
    int end = parse_segment(formula, 0, len, 1, &result);

    if (end < 0 || end != len) {
        memset(&result,0,sizeof(Compound)); 
    }
    return result;
}



