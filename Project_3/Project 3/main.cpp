//
//  main.cpp
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 21.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#define _USE_MATH_DEFINES // Mathematical constants

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>


// #include "lib.h"

using namespace std;

ofstream ofile;


void Eulers_Method(int, double, double, double, double*, double*, double*);
void Velocity_Verlet(int, double, double, double, double*, double*, double*);



int main(int argc, const char * argv[]) {
    // Set variable to decide which method should be implemented:
    // Eulers method = 0, Velocity Verlet = 1
    int method = 1;
    
    string outfile_name;
    
    int N = 1E4;
    double Final_time = 50.;
    double h = Final_time/N;
    
    // Inital values (given in [AU]):
    double x_initial = 1;
    double y_initial = 0;
    
    double *x = new double[N];
    double *y = new double[N];
    x[0] = x_initial; y[0] = y_initial;
    double v_x_0 = 0; double v_y_0 = 6.3;
    
    double *t = new double[N];
    
    
    if (method == 0){
        Eulers_Method(N, h, v_x_0, v_y_0, x, y, t);
        outfile_name = "eulers_method.txt";
    }
    else if (method == 1){
        Velocity_Verlet(N, h, v_x_0, v_y_0, x, y, t);
        outfile_name = "velocity_verlet.txt";
    }
    

    // Print to file
    
    ofile.open(outfile_name);
    
    ofile << "Results of program main.cpp, t - x - y" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0; i<N; i++)
    {
        ofile << setprecision(8) << t[i] << "  ";
        ofile << setprecision(8) << x[i] << "  ";
        ofile << setprecision(8) << y[i] << endl;
    }
    ofile.close();
    
    
    
    delete [] x;
    delete [] y;
    delete [] t;
}


// Function to solve an ordinary diffential equation using Eulers Method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void Eulers_Method(int N, double h, double initial_velocity_x, double initial_velocity_y, double *x, double *y, double *t) {
    double v_x_i = initial_velocity_x; double v_y_i = initial_velocity_y;
    double r_i; double velocity = 0;
    
    for (int i=0; i<N; i++){
        r_i = sqrt(x[i]*x[i] + y[i]*y[i]);
        x[i+1] = x[i] + h*v_x_i;
        y[i+1] = y[i] + h*v_y_i;
        t[i+1] = t[i] + h;
        // Update velocities
        velocity = (4*M_PI*M_PI)/(r_i*r_i*r_i);
        v_x_i -= h*velocity*x[i];
        v_y_i -= h*velocity*y[i];
    }
}



// Function to solve an ordinary diffential equation using the velocity Verlet method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void Velocity_Verlet(int N, double h, double initial_velocity_x, double initial_velocity_y, double *x, double *y, double *t){
    double hh = h*h;
    
    double G = 4*M_PI*M_PI;
    double M_sun = 1; //double M_earth = 0.000003;
    
    double a_x_new, a_y_new;
    double v_x_i = initial_velocity_x; double v_y_i = initial_velocity_y;
    
    double r = sqrt(x[0]*x[0] + y[0]*y[0]);
    double F = -G*M_sun/(r*r*r);
    double a_x_old = F*x[0]; double a_y_old = F*y[0];
    
    
    
    for (int i=0; i<N; i++){
        r = sqrt(x[i]*x[i] + y[i]*y[i]);

        x[i+1] = x[i] + h*v_x_i + hh/2*a_x_old;
        y[i+1] = y[i] + h*v_y_i + hh/2*a_y_old;
        t[i+1] = t[i] + h;
        
        F = -G*M_sun/(r*r*r); //Modified Force
        
        a_x_new = F*x[i+1]; a_y_new = F*y[i+1];
        
        v_x_i += h/2*(a_x_new + a_x_old);
        a_x_old = a_x_new;
        v_y_i += h/2*(a_y_new + a_y_old);
        a_y_old = a_y_new;
    }
}








