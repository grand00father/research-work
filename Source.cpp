#include <iostream>
using namespace std;

int main()
{
	int n=5;
	int i;
	int res;

	res = 1;
	for (i = 1; i <= n; i++) {
		res = res * i;
	}
	cout << res;
	return 0;
}