#pragma once
#include "types.h"

Compound parse_formula(const char* formula);

//validate
int validate_formula(const char* formula);

//helper
int is_upper(char c);
int is_lower(char c);
int is_digit(char c);