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

