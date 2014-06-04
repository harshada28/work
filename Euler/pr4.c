/*
A palindromic number reads the same both ways.
The largest palindrome made from the product of two 2-digit numbers is 9009 = 91 99.

Find the largest palindrome made from the product of two 3-digit numbers.
*/

#include <stdio.h>
#include <math.h>

int isPalindrome(int a)
{
    int reverse = 0, r, num;
    int i = 0;
    num = a;
    while (a) {
        r = a % 10;
        reverse = reverse * 10 + r;
        a = a/10;
        i++;
    }
    if (reverse == num) 
      return 1;
    else
      return 0;
}

int main()
{
    int a, x, y, max;
    max = 0;
    for (x = 999; x > 99; x--)
    {
        for (y = 999; y > 99; y--)
        {
          if (isPalindrome(x*y))
            {
              if (max < x*y)
                max = x*y;
              else
                break;
//              return 0;
            }
        }
    }
              printf("%d\n", max);
}
