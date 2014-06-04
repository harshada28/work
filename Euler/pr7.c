/*
By listing the first six prime numbers: 2, 3, 5, 7, 11, and 13,
we can see that the 6th prime is 13.

What is the 10 001st prime number? */

#include <stdio.h>
#include <math.h>

int isPrime(int x, int cnt)
{
  double sq;
  int y, i;

  if (x%2 == 0)
    return 0;
  if (x%5 == 0)
    return 0;
  sq = sqrt(x);
  y = sq + 0.9;
  for (i = 3; i < y; i++)
  {
    if (x%i == 0)
      return 0;
  }
  printf("%d: %d \n", cnt, x);
  return 1;

}

int main()
{
  int x = 17, y, cnt;

  cnt = 7;
  while (1) {
    if (cnt == 10002)
      break;
    if(isPrime(x, cnt))
      cnt++;
    x++;
  }
}
