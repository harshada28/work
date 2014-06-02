#include <stdio.h>

int fact(int no)
{
  unsigned long long int ans;
  ans = 1;
  while (no) {
    ans = ans * no;
    no--;
  printf("%lld\n", ans);
  }
  printf("%lld\n", ans);
}

int main()
{
  int t, no;
  scanf("%d", &t);
  while (t) {
    scanf("%d", &no);
    fact(no);
    t--;
  }

}
