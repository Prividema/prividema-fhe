#include <stdio.h>
#if GGSW 
    #include "ggsw.h"
#else
    #include "glwe.h"
#endif

int main() {
    int result;

    #if GGSW
        // Test addition
        result = add(1, 10);
        printf("1+10= %d\n", result);  // Expected output: 7

        // Test subtraction
        result = subtract(10, 5);
        printf("10 - 5 = %d\n", result);  // Expected output: 5
    #else
        printf("GGSW Not built.\n");
    #endif

    #if DEBUG
        printf("Hello Debug !\n");
    #endif
    return 0;
}