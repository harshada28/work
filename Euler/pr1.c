#include <stdio.h>

int main()
{
  int i, a1, a2, a3, sum;
  sum = 0;
  for (i = 0; i < 1000; i++)
  {
    a1 = i%3;
    a2 = i%5;
    a3 = i%15;
    if (!a1 || !a2 || !a3)
      sum += i;
  }
  printf("%d\n", sum);
}
