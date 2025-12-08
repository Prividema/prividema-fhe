#include "common/utils.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NB_BOXES 100
#define NB_SAMPLES 100000

// Chi-square critical value for alpha = 0.05 depending on the degrees of freedom
double chi_critical_05[100] = {
    3.841, 5.991, 7.815, 9.488, 11.070, 12.592, 14.067, 15.507, 16.919, 18.307, 19.675, 21.026, 22.362,
    23.685, 24.996, 26.296, 27.587, 28.869, 30.144, 31.410, 32.671, 33.924, 35.172, 36.415, 37.652, 38.885, 40.113, 
    41.337, 42.557, 43.773, 44.985, 46.194, 47.400, 48.602, 49.802, 50.998, 52.192, 53.384, 54.572, 55.758, 56.942, 
    58.124, 59.304, 60.481, 61.656, 62.830, 64.001, 65.171, 66.339, 67.505, 68.669, 69.832, 70.993, 72.153, 73.311, 
    74.468, 75.624, 76.778, 77.931, 79.082, 80.232, 81.381, 82.529, 83.675, 84.821, 85.965, 87.108, 88.250, 89.391, 
    90.531, 91.670, 92.808, 93.945, 95.081, 96.217, 97.351, 98.484, 99.617, 100.749, 101.879, 103.010, 104.139, 105.267, 
    06.395, 107.522, 108.648, 109.773, 110.898, 112.022, 113.145, 114.268, 115.390, 116.511, 117.632, 118.752, 119.871, 
    120.990, 122.108, 123.225, 124.342
};

int rand_uniform_aux()
{
    // boxes = {0, ..., 0}
    int64_t boxes[NB_BOXES];
    for(int i = 0; i < NB_BOXES; i++)
        boxes[i] = 0;
    
    // Fill the boxes
    double step = ((double)(INT64_MAX) - (double)INT64_MIN) / (double)NB_BOXES;
    for(size_t i = 0; i < NB_SAMPLES; i++) {
        int64_t sample = 0;
        if(rand_uniform(&sample) < 0)
            return 1;

       
        // This loop only fill NB_BOXES - 1 elements of boxes.
        int j = 0;
        for(double v = (double)INT64_MIN; v < (double)INT64_MAX; v += step) {
            if(sample < v + step) {
                boxes[j]++;
                break;
            }
            j++;
        }
    }

    // Apply Chi square test
    double expected = (double) NB_SAMPLES / (double) NB_BOXES; // For an uniform distribution
    double T = 0.0;
    for(int i = 0; i < NB_BOXES; i++) {
        double num = (((double)boxes[i]) - expected);
        T += (num * num) / expected;
    }

    // The test has 5% chance of failing
    if (T < chi_critical_05[NB_BOXES - 2])
        return 0;
    return 1;
}

// Test rand_uniform
Test(common, test_rand_uniform) 
{
    int count = 0;
    for(int i = 0; i < 100; i++)
        count += rand_uniform_aux();

    // As Chi square has 5% chance of failing, we count the number of times it fails.
    // On 100 iterations it fails 5 times in average.
    // We determined experimentally the standard deviation sigma = 4
    // So the number of failures should be in range 5 +- 4.
    cr_assert(le(int, 1, count), "The number of errors should be between in range 5 +- 4");
    cr_assert(ge(int, 9, count), "The number of errors should be between in range 5 +- 4");
}

// Test rand_normal
Test(common, test_rand_normal) 
{

}