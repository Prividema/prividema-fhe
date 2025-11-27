#include "distributions.h"
#include <iostream>
#include <random>
 
int main()
{
    RNG_UNI* uni = cpp_new_rng_uni(0,9);
 
    std::map<long, unsigned> histogram{};
    for (auto n{10000}; n; --n)
        ++histogram[cpp_draw_uniform(uni)];
 
    for (const auto [k, v] : histogram)
        std::cout << std::setw(2) << k << ' ' << std::string(v / 200, '*') << '\n';
}