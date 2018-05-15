#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main() {
	ifstream file3d;
	ofstream file2d;

	file3d.open("vel.dat");
	file2d.open("velxy.dat");

	int x,y,z;
	double v1, v2, v3;

	do {
		file3d>>x;
		file3d>>y;
		file3d>>z;
		file3d>>v1;
		file3d>>v2;
		file3d>>v3;
		if(z==8) {
			file2d<<x<<" "<<y<<" "<<v1/sqrt(v1*v1+v2*v2)<<" "<<v2/sqrt(v1*v1+v2*v2)<<endl;
		}
		cout<<x<<" "<<y<<" "<<z<<" "<<endl;
	} while(x<=199);

	file3d.close();
	file2d.close();

	return 0;
}