#include <stdio.h>

int main()
{
    unsigned long int a, sum, ans, pr;
    sum = ans = 0;
    pr = 1;
    printf("1 2 ");
    for (a = 2; a < 4000000;)
    {
//        if (a%2 == 0)
  //          sum += a;
        sum = a + pr;
        printf("%ld ", sum);
        if ((sum % 2 == 0) && sum < 4000000) {
            ans += sum;
        }
        pr = a;
        a = sum;
    }
    printf("!%ld\n", ans);
}
