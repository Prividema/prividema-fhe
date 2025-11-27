#include "distributions_first.h"

RNG_UNI* cpp_new_rng_uni(int uni_min, int uni_max){
    return new RNG_UNI(uni_min, uni_max);
}

void cpp_delete_rng_uni(RNG_UNI* p) {
    delete p;
}

int cpp_draw_uniform(RNG_UNI* uni){
    return uni->drawUniformInt();
}