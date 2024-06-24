#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void info_out (FILE * fp, 
	      double * rx, double * ry, double * rz,  //position coords
	      double * vx, double * vy, double * vz,  //velocity coords
	      int * ix, int * iy, int * iz, double L, //arrays of integers that keep track of the no. of times each particle has crossed PBC
	      int N, int put_vel, int unfold) {
    int i;

    //Print number of particles
    fprintf(fp,"%i\n",N);

    //Print the simulation box dimensions
    fprintf(fp,"BOX %.5lf %.5lf %.5lf\n",L,L,L);

    //Loop over each particle
    for (i=0;i<N;i++) {
        //Print particle data with coordinates with unfolding
        if (unfold) {
            fprintf(fp,"%i % 10.6lf % 10.6lf % 10.6lf ",
                rx[i] + (ix[i] * L),
                ry[i] + (iy[i] * L),
                rz[i] + (iz[i] * L));
        }
        else{
            fprintf(fp,"%i % 10.6lf % 10.6lf % 10.6lf ", rx[i], ry[i], rz[i]);
        }
        //If velocities should be included, print them
        if (put_vel)
        fprintf(fp,"% 10.6lf % 10.6lf % 10.6lf",vx[i],vy[i],vz[i]);
        
        //End the line for the current particle
        fprintf(fp,"\n");
    }
}

int info_in (FILE * fp, double * rx, double * ry, double * rz, 
	     double * vx, double * vy, double * vz, double * L,
	     int * N) {
    int i, l;
    int has_vel = 0, dum;
    char dummy[4];
    double Lx, Ly, Lz;

    //Read number of particles and velocity flag from file
    if (fscanf(fp, "%i %i\n", N, &has_vel) != EOF) {
        //Read box dimensions
        //dummy to discard the BOX string
        i = fscanf(fp, "%s %lf %lf %lf\n", dummy, &Lx, &Ly, &Lz);
        *L = Lx;

        //Loop over each particle
        for (i = 0; i < (*N); i++) {
        //Read particle identifier and coordinates
        //dum used to discard atomic number
        l = fscanf(fp, "%i %lf %lf %lf ", &dum, &rx[i], &ry[i], &rz[i]);

        //If velocities are included, read them
        if (has_vel) {
            l = fscanf(fp, "%lf %lf %lf", &vx[i], &vy[i], &vz[i]);
        }
        }
    }

    //Return whether velocities were included
    return has_vel;
}

