/*
The sum of the primes below 10 is 2 + 3 + 5 + 7 = 17.

Find the sum of all the primes below two million.
*/

#include <stdio.h>
#include <math.h>

int isPrime(int x)
{
    int y, i;
    double sq;
    sq = sqrt(x);
    y = sq + 0.9;
   
    if (x%2 == 0)
        return 0;
    if (x%5 == 0)
        return 0; 
    for (i = 2; i <= y; i++)
    {
        if (x%i == 0)
            return 0;
    }
    printf("*%d ", x);
    return 1;

}

int main()
{
    int cnt, x;
    unsigned long long int sum;
    sum = 0; cnt = 11;
    x = 0;
    while (cnt < 2000000)
    {
        if(isPrime(cnt))
        {
            x++;
            sum += cnt;
            printf("%lld \n", sum);
        }
        cnt++;
    }
    sum += 17;
    printf("%d %lld %d\n", cnt, sum, x);
    return 0;
}
