//
//  ODEsolver.cpp
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 21.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include "ODEsolver.h"
#include "planet.h"

using namespace std;

// Initializers
ODEsolver::ODEsolver(){
    G = 4*M_PI*M_PI;
    v_x_i = 0; v_y_i = 6.3;
    Kinetic_Energy = 0.0;
    Potential_Energy = 0.0;
    Total_Energy = 0.0;
}

ODEsolver::ODEsolver(double initial_velocity_x, double initial_velocity_y){
    G = 4*M_PI*M_PI;
    v_x_i = initial_velocity_x; v_y_i = initial_velocity_y;
    Kinetic_Energy = 0.0;
    Potential_Energy = 0.0;
    Total_Energy = 0.0;
}



// Function to solve an ordinary diffential equation using Eulers Method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Eulers_Method(int N, double h, double *x, double *y, double *t){
    double r = sqrt(x[0]*x[0] + y[0]*y[0]);
    double velocity = 0;
    
    double M_earth = 0.000003;
    
    
    cout << setprecision(4) << scientific;
    Energy(M_earth, r);
    cout << "Start: " << endl;
    print_energy();
    
    clock_t start, end;
    start=clock();
    
    for (int i=0; i<N; i++){
        r = sqrt(x[i]*x[i] + y[i]*y[i]);
        x[i+1] = x[i] + h*v_x_i;
        y[i+1] = y[i] + h*v_y_i;
        t[i+1] = t[i] + h;
        // Update velocities
        velocity = G/(r*r*r);
        v_x_i -= h*velocity*x[i];
        v_y_i -= h*velocity*y[i];
    }
    
    end=clock();
    
    Energy(M_earth, r);
    cout << "End:   " << endl;
    print_energy();
    
    cout << scientific << "CPU time for Eulers method (sec): " << ((double)end-(double)start)/CLOCKS_PER_SEC << endl;
}



// Function to solve an ordinary diffential equation using the velocity Verlet method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Velocity_Verlet(int N, double h, double *x, double *y, double *t){
    double hh = h*h;

    double M_sun = 1; double M_earth = 0.000003;
    
    double r = sqrt(x[0]*x[0] + y[0]*y[0]);
    double F = -G*M_sun/(r*r*r);                        // Modified (F = -G*M_sun*M_earth/(r*r) and a_x = (F/M_earth) * x/r)
    double a_x_old = F*x[0]; double a_y_old = F*y[0];
    double a_x_new, a_y_new;
    
    cout << setprecision(4) << scientific;
    Energy(M_earth, r);
    cout << "Start: " << endl;
    print_energy();
    
    clock_t start, end;
    start=clock();

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
    
    end=clock();
    
    Energy(M_earth, r);
    cout << "End:   " << endl;
    print_energy();
    
    cout << scientific << "CPU time for the Verlet Velocity method (sec): " << ((double)end-(double)start)/CLOCKS_PER_SEC << endl;
}


// Function to find kinetic and potential energy
void ODEsolver::Energy(double mass_of_planet, double r){
    planet planet(mass_of_planet);
    Kinetic_Energy = planet.Kinetic_Energy(v_x_i, v_y_i);
    Potential_Energy = planet.Potential_Energy(G, r);
    Total_Energy = Kinetic_Energy+Potential_Energy;
}

// Print energy from ODEsolver::Energy
void ODEsolver::print_energy(){
    cout << "Kinetic Energy = " << Kinetic_Energy << " --- ";
    cout << "Potential Energy = " << Potential_Energy << endl;
    cout << "Total Energy = " << Total_Energy << endl;
}





