/* The prime factors of 13195 are 5, 7, 13 and 29.

What is the largest prime factor of the number 600851475143 ?*/

#include <stdio.h>

int main()
{
  unsigned long int a, x, r, max_prime;
  a = 2;
  x = 600851475143;
  max_prime = 1;
  while (a <= x) {
    r = x % a;
    if (r == 0 && max_prime < a)
      max_prime = a;
    a++;
  }
  printf("%ld\n", max_prime);
}

