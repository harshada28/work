#include <stdio.h>
#include <math.h>

int isPalindrome(int a)
{
    int reverse = 0, r;
    int i = 0;
    while (a) {
        r = a % 10;
        //reverse = reverse * pow(10, i) + (r);
        reverse = 3 * pow(10, i) + r;
        printf("*%d\n", reverse); 
        a = a/10;
        i++;
    }
    printf("%d\n", reverse);
    
}

int main()
{
    int a;
    a = 99 * 99;
/*    while(1) {
        if (isPalindrome(a))
            break;
    }
*/
    isPalindrome(456);
}
