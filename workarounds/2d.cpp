#include <iostream>
#include <fstream>

using namespace std;

int main() {
	ifstream file3d;
	ofstream file2d;

	file3d.open("stfn.dat");
	file2d.open("stfnxy.dat");

	int x,y,z;
	double st;

	do {
		file3d>>x;
		file3d>>y;
		file3d>>z;
		file3d>>st;
		if(z==8) {
			file2d<<st<<" ";
			if(y==15) {
				file2d<<endl;
			}
		}
		cout<<x<<" "<<y<<" "<<z<<" "<<st<<endl;
	} while(x<=199);

	file3d.close();
	file2d.close();

	return 0;
}