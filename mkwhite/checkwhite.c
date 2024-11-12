#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>


int main(int argc, char **argv){
	int i,j,k;
	int nx,ny,nz;
	int ichip, jchip,iseed;

	FILE *fp = fopen("white.out","r");
	fread(&ichip, sizeof(int), 1, fp);
	fread(&nx, sizeof(int), 1, fp);
	fread(&ny, sizeof(int), 1, fp);
	fread(&nz, sizeof(int), 1, fp);
	fread(&iseed, sizeof(int), 1, fp);
	fread(&ichip, sizeof(int), 1, fp);
	float *den = (float*)malloc(sizeof(float)*nx*ny);

	for(k=0;k<nz;k++){
		fread(&ichip, sizeof(int), 1, fp);
		fread(den,sizeof(float),nx*ny,fp);
		fread(&jchip, sizeof(int), 1, fp);
		float mean,std;
		mean = std = 0;
		for(i=0;i<nx*ny;i++){
			mean += den[i];
			std += den[i]*den[i];
		}
		mean = mean/nx/ny;
		std = std/nx/ny - mean*mean;

		printf("k=%d; chips= %d %d with mean/std = %g %g\n", k,ichip, jchip, mean,std);
	}
}
