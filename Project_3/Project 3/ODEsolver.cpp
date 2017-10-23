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

ofstream ofile;

// Initializers
ODEsolver::ODEsolver(){
    G = 4*M_PI*M_PI;
    number_of_planets = 0;
    Kinetic_Energy = 0.0;
    Potential_Energy = 0.0;
    Total_Energy = 0.0;
}


void ODEsolver::Add_Planet(planet new_planet){
    number_of_planets += 1;
    all_planets.push_back(new_planet);
}


// Function to solve an ordinary diffential equation using Eulers Method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Eulers_Method(int N, double h, double *x, double *y, double *t){
    /*
    double r = sqrt(x[0]*x[0] + y[0]*y[0]);
    double velocity = 0;
    
    //double M_earth = 0.000003;
    
    
    cout << setprecision(4) << scientific;
    Energy(M_earth, x[N], y[N], v_x_i, v_y_i);
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
    
    
    Energy(M_earth, x[N], y[N], v_x_i, v_y_i);
    cout << "End:   " << endl;
    print_energy();
 

    cout << scientific << "CPU time for Eulers method (sec): " << ((double)end-(double)start)/CLOCKS_PER_SEC << endl;
     */
 }



// Function to solve an ordinary diffential equation using the velocity Verlet method
// Takes initial values for x and y as first values of arrays, and fills arrays x and y with corresponing positions.
void ODEsolver::Velocity_Verlet(int N, double h, string outfile_name){
    ofile.open(outfile_name);
    
    ofile << "Positions for planets (x, y)" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    
    double hh = h*h;

    double F, Fx, Fy; double X, Y;
    double a_x_old, a_y_old, a_x_new, a_y_new;

    
    clock_t start, end;
    start=clock();
    
    
    for (int i=0; i<N; i++){
        
        for (int nr_1=0; nr_1<number_of_planets; nr_1++){
            planet &current = all_planets[nr_1];
            Fx = 0.0; Fy = 0.0;
            
            ofile << setprecision(8) << current.x << "  ";
            ofile << setprecision(8) << current.y << "  ";
            
            for (int nr_2=nr_1+1; nr_2<number_of_planets; nr_2++){
                planet &other = all_planets[nr_2];
                F = current.Gravitational_Force(other, G);
                X = current.x - other.x; Y = current.y - other.y;
                Fx -= F*X/current.distance(other); Fy -= F*Y/current.distance(other);
            }
            
            a_x_old = Fx/current.mass_of_planet; a_y_old = Fy/current.mass_of_planet;
            current.x += h*current.v_x + hh/2*a_x_old;
            current.y += h*current.v_y + hh/2*a_y_old;
            
            Fx = 0.0; Fy = 0.0;
            
            for (int nr_2=nr_1+1; nr_2<number_of_planets; nr_2++){
                planet &other = all_planets[nr_2];
                F = current.Gravitational_Force(other, G);
                X = current.x - other.x; Y = current.y - other.x;
                Fx -= F*X/current.distance(other); Fy -= F*Y/current.distance(other);
            }
            
            a_x_new = Fx/current.mass_of_planet; a_y_new = Fy/current.mass_of_planet;
            current.v_x += h/2*(a_x_new + a_x_old);
            current.v_y += h/2*(a_y_new + a_y_old);
            
        }
        ofile << endl;
        
    }
    end=clock();
    
    ofile.close();
    
}


// Function to find kinetic and potential energy
void ODEsolver::Energy(planet planet, double x, double y, double vx, double vy){
    Kinetic_Energy = planet.Kinetic_Energy(vx, vy);
    double r = sqrt(x*x + y*y);
    Potential_Energy = planet.Potential_Energy(G, r);
    Total_Energy = Kinetic_Energy+Potential_Energy;
}

// Print energy from ODEsolver::Energy
void ODEsolver::print_energy(){
    cout << "Kinetic Energy = " << Kinetic_Energy << " --- ";
    cout << "Potential Energy = " << Potential_Energy << endl;
    cout << "Total Energy = " << Total_Energy << endl;
}





