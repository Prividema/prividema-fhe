#include <stdio.h>
#include "ggsw.h"

int main() {
    int result;

    // Test addition
    result = add(1, 10);
    printf("1+10= %d\n", result);  // Expected output: 7

    // Test subtraction
    result = subtract(10, 5);
    printf("10 - 5 = %d\n", result);  // Expected output: 5

    #if DEBUG
        printf("Hello Debug !\n");
    #endif
    return 0;
}