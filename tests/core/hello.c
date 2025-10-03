#include <stdio.h>

int main() {
    int result;

     // Test addition
    int a=1, b=10;
    result = a+b;
    printf("%d+%d= %d\n",a, b, result);  

    // Test subtraction
    result = a-b;
    printf("%d-%d= %d\n",a, b, result);   
    
    #if DEBUG
        printf("Hello Debug !\n");
    #endif
    return 0;
}