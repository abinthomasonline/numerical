#include <iostream>
#include <omp.h>

using namespace std;

int main() {
	int a[13] = {1,2,3,9,8,3,5,2,86,992,885,193,123};
	int l=0;
	#pragma omp parallel for
	for(int i=0;i<13;++i) {
		#pragma omp critical
		if(l<a[i]) {
			l=a[i];
		}
	}
	cout<<l;
	return 0;
}