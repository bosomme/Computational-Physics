//
//  ODEsolver.cpp
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 21.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include "ODEsolver.h"

// Initializers
ODEsolver::ODEsolver(){
    G = 4*M_PI*M_PI;
    v_x_i = 0; v_y_i = 6.3;
}

ODEsolver::ODEsolver(double initial_velocity_x, double initial_velocity_y){
    G = 4*M_PI*M_PI;
    v_x_i = initial_velocity_x; v_y_i = initial_velocity_y;
}



// Function to solve an ordinary diffential equation using Eulers Method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Eulers_Method(int N, double h, double *x, double *y, double *t){
    double r_i; double velocity = 0;
    
    for (int i=0; i<N; i++){
        r_i = sqrt(x[i]*x[i] + y[i]*y[i]);
        x[i+1] = x[i] + h*v_x_i;
        y[i+1] = y[i] + h*v_y_i;
        t[i+1] = t[i] + h;
        // Update velocities
        velocity = G/(r_i*r_i*r_i);
        v_x_i -= h*velocity*x[i];
        v_y_i -= h*velocity*y[i];
    }
}



// Function to solve an ordinary diffential equation using the velocity Verlet method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Velocity_Verlet(int N, double h, double *x, double *y, double *t){
    double hh = h*h;

    double M_sun = 1; //double M_earth = 0.000003;
    
    double r = sqrt(x[0]*x[0] + y[0]*y[0]);
    double F = -G*M_sun/(r*r*r);
    double a_x_old = F*x[0]; double a_y_old = F*y[0];
    double a_x_new, a_y_new;
    

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





