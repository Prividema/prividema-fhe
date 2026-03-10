#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw_key.h"

#define NBASE            1024
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 2
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 2
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3

//! bivGGSW KEY PART (begin)

