#ifndef GGSW_SECRET_KEY_H
#define GGSW_DECRET_KEY_H

#include <cstdint>
//#include something for INT_POL

typedef struct gsw_secret_key{
    INT_POL* values; // The key itself.
} GSW_SECRET_KEY;

#endif