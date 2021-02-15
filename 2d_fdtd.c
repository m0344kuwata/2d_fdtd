#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//constants//
#define NX 150				//X_cell_number
#define NY 150		//Y_cell_number
#define NT 500			//time_number
#define DEL_T 0.00000000001473	 //delta_time
#define DEL_X 0.006246			//delta_x
#define DEL_Y 0.006246		//delta_y
#define E0 0.000000000008854	//epsilon0
#define M0 0.000001257			//mu0
#define TAU 0.000000001				//timeconstant
#define E 0						//flag of e-field
#define H 1						//flag of h-field
#define NM 5					//material_number
#define FREE 0					//free space

// DEL定数の決め方　⇒　1.5.1節，（1.79）（1.80）参照

//functions//
void set_material(void);	//今回は使用せず．Vaccume以外の材質の場合に使用
void cal_ezfld(double ez[NX][NY],  double hx[NX][NY], double hy[NX][NY], int step, double cez[NX][NY], double cezlx[NX][NY], double cezly[NX][NY]) ;	//calculate ez-field
void Mur_surface(double ez[NX][NY],  double ezb[NX][NY], double cxdx, double cxdy);	//boundary cordition
void cal_hfld(double ez[NX][NY], double hy[NX][NY], double hx[NX][NY], double chxly[NX][NY], double chylx[NX][NY]);				//calculate h-field
double current_source(double);					//current_source
void print_field(double f[NX][NY], int flag, int step) ;			//output


//main//
int main() {

	//parameters//
	double ez[NX][NY] = { 0.0 }, hx[NX][NY] = { 0.0 }, hy[NX][NY] = { 0.0 };			//ez,hx,hy fields
	int n, i;												//number of timestep
	double ezb[NX][NY] = { 0.0 };								//before ez-filed
	double v0 = 1 / sqrt(E0*M0);							//velocity0
	double cxdx = (v0*DEL_T - DEL_X) / (v0*DEL_T + DEL_X);	//constant of velocity0
	double cxdy = (v0*DEL_T - DEL_Y) / (v0*DEL_T + DEL_Y);	//constant of velocity0


	int idex[NX][NY] = { FREE };									//material id
	double eps[NM] = { E0 };
	double sigm[NM] = { 0.0 };
	double mu[NM] = { M0 };

	double cez[NX][NY];
	double cezlx[NX][NY];										//constant of (1.26)(1.30)
	double cezly[NX][NY];
	double chxly[NX][NY];
	double chylx[NX][NY];

	double cez_0 = 1;													//constant where sigma = 0
	double cezlx_0 = DEL_T / DEL_X / E0;
	double cezly_0 = DEL_T / DEL_Y / E0;
	double chxly_0 = DEL_T / M0 / DEL_Y;
  double chylx_0 = DEL_T / M0 / DEL_X;


	set_material();

	//calculate constants//
	for (int i = 0; i < NX; i++) {
		for(int j = 0; j < NY; j++)
		if (idex[i][j] == FREE) {
			cez[i][j] = cez_0;
			cezlx[i][j] = cezlx_0;
			cezly[i][j] = cezly_0;
			chxly[i][j] = chxly_0;
			chylx[i][j] = chylx_0;
		}
		//材質がVaccume以外であればここで処理
	}

	//main loop//
	for (n = 1; n < NT; n++) {

		for (int i = 0; i < NX ; i++) {
			for (int j = 0; j < NY ; j++) {
				ezb[i][j] =ez[i][j];
			}
		}

		cal_ezfld(ez, hx, hy, n, cez, cezlx, cezly);								//calculate ez-field
		printf("n=%d\n",n);

		Mur_surface(ez ,ezb ,cxdx, cxdy);													//calculate ez surface

		cal_hfld(ez, hy, hx, chxly, chylx);											//calculate h-field

		print_field(ez, E, n);
																																//output e-field
		printf("finish\n");
	}

	return 0;
}


void set_material(void){

	printf("All field is Vaccume\n");

}



void cal_ezfld(double ez[NX][NY],  double hx[NX][NY], double hy[NX][NY], int step, double cez[NX][NY], double cezlx[NX][NY], double cezly[NX][NY]) {

	//calculate ex-field//
	for (int i = 1; i < NX - 1; i++) {
		for (int j = 1; j < NY - 1; j++) {
			ez[i][j] = (cez[i][j] * ez[i][j]) + (cezlx[i][j] * (hy[i][j] - hy[i-1][j])) - (cezly[i][j] * (hx[i][j] - hx[i][j-1]));
		}
	}

	//current_source//
	ez[NX / 2][NY / 2] = current_source((double)step*DEL_T);

}

void Mur_surface(double ez[NX][NY],  double ezb[NX][NY], double cxdx, double cxdy){

	// x = 0 面
	for(int j = 1; j < NY - 1; j++){
		ez[0][j] = ezb[1][j] + cxdx * (ez[1][j] - ezb[0][j]);
		ez[NX-1][j] = ezb[NX-2][j] + cxdx * (ez[NX-2][j] - ezb[NX-1][j]);
	}
	// y = 0 面
	for(int i = 1; i < NX - 1; i++){
		ez[i][0] = ezb[i][1] + cxdy * (ez[i][1] - ezb[i][0]);
		ez[i][NY-1] = ezb[i][NY-2] + cxdy * (ez[i][NY-2] - ezb[i][NY-1]);
	}

}


//calculate h-field//
void cal_hfld(double ez[NX][NY], double hy[NX][NY], double hx[NX][NY], double chxly[NX][NY], double chylx[NX][NY]) {

	//calculate hx-field//
	for (int i = 1; i < NX - 1; i++) {
		for (int j = 0; j < NY - 1; j++) {
		hx[i][j] = hx[i][j] - chxly[i][j] * (ez[i][j+1]-ez[i][j]);
		}
	}
	//calculate hy-field//
	for (int i = 0; i < NX - 1; i++) {
		for (int j = 1; j < NY - 1; j++) {
		hy[i][j] = hy[i][j] + chylx[i][j] * (ez[i+1][j]-ez[i][j]);
		}
	}

}

//current_source//
double current_source(double t) {

	//gaussian pulse//
	if (t >= 0 && t <= 2 * TAU) {

		return (exp(-1 * (4 / TAU)*(4 / TAU)*(t - TAU)*(t - TAU)));
	}
	else {
		return 0;
	}

}

//output//
void print_field(double f[NX][NY], int flag, int step) {

	//parameters//
	int i,j;
	FILE *fp;							//file pointer
	char fname[] = "e_field_10000.txt";	//file name

	//label filed//
	if (flag == E) {
		fname[0] = 'e';
	}
	else {
		fname[0] = 'h';
	}

	//numbering//
	fname[9] = step % 10000 / 1000 + 48;
	fname[10] = step % 1000 / 100 + 48;
	fname[11] = step % 100 / 10 + 48;
	fname[12] = step % 10 / 1 + 48;

	//file open//
	fp = fopen(fname, "w");
	if (fp == NULL) {
		printf("%s file not open.\n", fname);
		exit(-1);
	}
	else {
		printf("%s file opened.\n", fname);
	}

	//output//
	for (i = 0; i < NX; i++) {
		for(j = 0; j < NY; j++){
			fprintf(fp, "%e\t%e\t%e\n", (double)i*DEL_X, (double)j*DEL_Y, f[i][j]);	//z and field
		}
	}

	//file close//
	fclose(fp);

}
