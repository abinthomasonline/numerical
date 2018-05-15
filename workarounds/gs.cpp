#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void poisson(double*** phi, double*** f, int l, int b, int h, double dx, double* bc, double error_max) {
	double**** p = new double***[2];

	#pragma omp parallel
	{
		#pragma omp for
		for(int i1=0;i1<2;++i1) {
			p[i1] = new double**[l];
		}
		#pragma omp for collapse(2)
		for(int i1=0;i1<2;++i1) {
			for(int i2=0;i2<l;++i2) {
				p[i1][i2] = new double*[b];
			}
		}
		#pragma omp for collapse(3)
		for(int i1=0;i1<2;++i1) {
			for(int i2=0;i2<l;++i2) {
				for(int i3=0;i3<b;++i3) {
					p[i1][i2][i3] = new double[h];
				}
			}
		}
		#pragma omp for collapse(4)
		for(int i1=0;i1<2;++i1) {
			for(int i2=0;i2<l;++i2) {
				for(int i3=0;i3<b;++i3) {
					for(int i4=0;i4<h;++i4) {
						p[i1][i2][i3][i4]=phi[i2][i3][i4];
					}
				}
			}
		}
	}

	int i0=1;
	double error;
	do {
		error=0;
		#pragma omp parallel
		{
			#pragma omp for collapse(3)
			for(int i1=1;i1<l-1;++i1) {
				for(int i2=1;i2<b-1;++i2) {
					for(int i3=1;i3<h-1;++i3) {
						p[(i0+1)/2][i1][i2][i3]=(p[(-i0+1)/2][i1+1][i2][i3]+
												   p[(-i0+1)/2][i1-1][i2][i3]+
												   p[(-i0+1)/2][i1][i2+1][i3]+
												   p[(-i0+1)/2][i1][i2-1][i3]+
												   p[(-i0+1)/2][i1][i2][i3+1]+
												   p[(-i0+1)/2][i1][i2][i3-1]-
												   f[i1][i2][i3]*dx*dx)/6;
						double error_ijk=abs(p[(i0+1)/2][i1][i2][i3]-p[(-i0+1)/2][i1][i2][i3]);
						#pragma omp critical
						if(error_ijk>error) {
							error=error_ijk;
						}
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<b;++i1) {
				for(int i2=0;i2<h;++i2) {
					p[(i0+1)/2][0][i1][i2]=p[(i0+1)/2][1][i1][i2]-bc[0]*dx;
					double error_ijk=abs(p[(i0+1)/2][0][i1][i2]-p[(-i0+1)/2][0][i1][i2]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<b;++i1) {
				for(int i2=0;i2<h;++i2) {
					p[(i0+1)/2][l-1][i1][i2]=p[(i0+1)/2][l-2][i1][i2]+bc[3]*dx;
					double error_ijk=abs(p[(i0+1)/2][l-1][i1][i2]-p[(-i0+1)/2][l-1][i1][i2]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<h;++i1) {
				for(int i2=0;i2<l;++i2) {
					p[(i0+1)/2][i2][0][i1]=p[(i0+1)/2][i2][1][i1]-bc[1]*dx;
					double error_ijk=abs(p[(i0+1)/2][i2][0][i1]-p[(-i0+1)/2][i2][0][i1]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<h;++i1) {
				for(int i2=0;i2<l;++i2) {
					p[(i0+1)/2][i2][b-1][i1]=p[(i0+1)/2][i2][b-2][i1]+bc[4]*dx;
					double error_ijk=abs(p[(i0+1)/2][i2][b-1][i1]-p[(-i0+1)/2][i2][b-1][i1]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<l;++i1) {
				for(int i2=0;i2<b;++i2) {
					p[(i0+1)/2][i1][i2][0]=p[(i0+1)/2][i1][i2][1]-bc[2]*dx;
					double error_ijk=abs(p[(i0+1)/2][i1][i2][0]-p[(-i0+1)/2][i1][i2][0]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
			#pragma omp for collapse(2) nowait
			for(int i1=0;i1<l;++i1) {
				for(int i2=0;i2<b;++i2) {
					p[(i0+1)/2][i1][i2][h-1]=p[(i0+1)/2][i1][i2][h-2]+bc[5]*dx;
					double error_ijk=abs(p[(i0+1)/2][i1][i2][h-1]-p[(-i0+1)/2][i1][i2][h-1]);
					#pragma omp critical
					if(error_ijk>error) {
						error=error_ijk;
					}
				}
			}
		}
		i0*=-1;
		cout<<"error: "<<error<<endl;
	} while(error>error_max);

	#pragma omp parallel
	{
		#pragma omp for collapse(3)
		for(int i1=0;i1<l;++i1) {
			for(int i2=0;i2<b;++i2) {
				for(int i3=0;i3<h;++i3) {
					phi[i1][i2][i3]=p[(-i0+1)/2][i1][i2][i3];
				}
			}
		}

		#pragma omp for collapse(3)
		for(int i1=0;i1<2;++i1) {
			for(int i2=0;i2<l;++i2) {
				for(int i3=0;i3<b;++i3) {
					delete [] p[i1][i2][i3];
				}
			}
		}
		#pragma omp for collapse(2)
		for(int i1=0;i1<2;++i1) {
			for(int i2=0;i2<l;++i2) {
				delete [] p[i1][i2];
			}
		}
		#pragma omp for
		for(int i1=0;i1<2;++i1) {
			delete [] p[i1];
		}
	}
	delete [] p;
}

void get_velocity(double*** phi, double**** vel, int l, int b, int h, double dx) {
	#pragma omp parallel
	{
		#pragma omp for collapse(3) nowait
		for(int i1=1;i1<l-1;++i1) {
			for(int i2=1;i2<b-1;++i2) {
				for(int i3=1;i3<h-1;++i3) {
					vel[0][i1][i2][i3]=(phi[i1+1][i2][i3]-phi[i1-1][i2][i3])/2;
				}
			}
		}
		#pragma omp for collapse(3) nowait
		for(int i1=1;i1<l-1;++i1) {
			for(int i2=1;i2<b-1;++i2) {
				for(int i3=1;i3<h-1;++i3) {
					vel[1][i1][i2][i3]=(phi[i1][i2+1][i3]-phi[i1][i2-1][i3])/2;
				}
			}
		}
		#pragma omp for collapse(3) nowait
		for(int i1=1;i1<l-1;++i1) {
			for(int i2=1;i2<b-1;++i2) {
				for(int i3=1;i3<h-1;++i3) {
					vel[2][i1][i2][i3]=(phi[i1][i2][i3+1]-phi[i1][i2][i3-1])/2;
				}
			}
		}
	}
}

int main() {
	const int length=200, breadth=15, height=15;
	const double grid_size=0.1;
	const double error_max=0.1;

	double bc[6] = {};
	double*** phi = new double**[length];
	double*** source_sink = new double**[length];
	double**** velocity = new double***[3];
	#pragma omp parallel
	{
		#pragma omp for
		for(int i1=0;i1<length;++i1) {
			phi[i1] = new double*[breadth];
			source_sink[i1] = new double*[breadth];
		}
		#pragma omp for collapse(2)
		for(int i1=0;i1<length;++i1) {
			for(int i2=0;i2<breadth;++i2) {
				phi[i1][i2] = new double[height];
				source_sink[i1][i2] = new double[height];
			}
		}
		#pragma omp for collapse(3) nowait
		for(int i1=0;i1<length;++i1) {
			for(int i2=0;i2<breadth;++i2) {
				for(int i3=0;i3<height;++i3) {
					phi[i1][i2][i3] = 0;
				}
			}
		}
		#pragma omp for
		for(int i1=0;i1<3;++i1) {
			velocity[i1] = new double**[length];
		}
		#pragma omp for collapse(2)
		for(int i1=0;i1<3;++i1) {
			for(int i2=0;i2<length;++i2) {
				velocity[i1][i2] = new double*[breadth];
			}
		}
		#pragma omp for collapse(3)
		for(int i1=0;i1<3;++i1) {
			for(int i2=0;i2<length;++i2) {
				for(int i3=0;i3<breadth;++i3) {
					velocity[i1][i2][i3] = new double[height];
				}
			}
		}
	}

	source_sink[40][7][7]=100;
	source_sink[159][7][7]=-100;
	bc[0]=bc[3]=100;

	poisson(phi, source_sink, length, breadth, height, grid_size, bc, error_max);
	get_velocity(phi, velocity, length, breadth, height, grid_size);

	ofstream stfn_file;
	stfn_file.open("stfn.dat");
	for(int i1=0;i1<length;++i1) {
		for(int i2=0;i2<breadth;++i2) {
			for(int i3=0;i3<height;++i3) {
				stfn_file<<i1+1<<" "<<i2+1<<" "<<i3+1<<" "<<phi[i1][i2][i3]<<endl;
			}
		}
		stfn_file<<endl;
	}
	stfn_file.close();

	ofstream vel_file;
	vel_file.open("vel.dat");
	for(int i1=1;i1<length-1;++i1) {
		for(int i2=1;i2<breadth-1;++i2) {
			for(int i3=1;i3<height-1;++i3) {
				vel_file<<i1+1<<" "<<i2+1<<" "<<i3+1<<" "<<velocity[0][i1][i2][i3]<<" "<<velocity[1][i1][i2][i3]<<" "<<velocity[2][i1][i2][i3]<<endl;
			}
		}
	}
	vel_file.close();

	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for(int i1=0;i1<length;++i1) {
			for(int i2=0;i2<breadth;++i2) {
				delete [] phi[i1][i2];
				delete [] source_sink[i1][i2];
			}
		}
		#pragma omp for nowait
		for(int i1=0;i1<length;++i1) {
			delete [] phi[i1];
			delete [] source_sink[i1];
		}
		#pragma omp for collapse(3)
		for(int i1=0;i1<3;++i1) {
			for(int i2=0;i2<length;++i2) {
				for(int i3=0;i3<breadth;++i3) {
					delete [] velocity[i1][i2][i3];
				}
			}
		}
		#pragma omp for collapse(2)
		for(int i1=0;i1<3;++i1) {
			for(int i2=0;i2<length;++i2) {
				delete [] velocity[i1][i2];
			}
		}
		#pragma omp for
		for(int i1=0;i1<3;++i1) {
			delete [] velocity[i1];
		}
	}
	delete [] phi;
	delete [] source_sink;
	delete [] velocity;

	return 0;
}