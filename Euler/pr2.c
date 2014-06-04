/*
Each new term in the Fibonacci sequence is generated by adding the previous two terms. 
By starting with 1 and 2, the first 10 terms will be:

1, 2, 3, 5, 8, 13, 21, 34, 55, 89, ...

By considering the terms in the Fibonacci sequence whose values do not exceed four million, find the sum of the even-valued terms.
*/
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
