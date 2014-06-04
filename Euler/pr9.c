/*
A Pythagorean triplet is a set of three natural numbers, a  b  c, for which,

a2 + b2 = c2
For example, 32 + 42 = 9 + 16 = 25 = 52.

There exists exactly one Pythagorean triplet for which a + b + c = 1000.
Find the product abc.*/

#include <stdio.h>

int main()
{
  int i, j, x;
  for (i = 200; i < 1000; i++) {
    for (j = i+1; j < 1000; j++) {
      x = i * j - 1000 * j - 1000 * i + 500000;
      if (!x) {
        printf("\n\n%d %d\n", i, j);
        return;
      }
    }

  }
}
