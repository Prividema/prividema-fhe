#ifndef DISTRIBUTIONS_H
#define DISTRIBTUIONS_H


#ifdef __cplusplus
    #include <random>

    class RNG_UNI {
    public:
        RNG_UNI(int uni_min, int uni_max)
        : gen(std::random_device{}()),                     // initialise le générateur
          uniform(uni_min, uni_max)                    // initialise la distribution normale
        { }

        void setUniform(int a, int b) {
        uniform.param(std::uniform_int_distribution<>::param_type(a, b));
        }

        int drawUniformInt() {
            return uniform(gen);
        }
    private:
        std::mt19937 gen;
        std::uniform_int_distribution<> uniform;
    };
#else 
    typedef struct RNG_UNI RNG_UNI;
#endif

#ifdef __cplusplus
extern "C" {
#endif

    extern RNG_UNI* cpp_new_rng_uni(int uni_min, int uni_max);

    extern void cpp_delete_rng_uni(RNG_UNI* p) ;

    extern int cpp_draw_uniform(RNG_UNI* uni);

#ifdef __cplusplus
}
#endif

#endif /*DISTRIBUTIONS*/
