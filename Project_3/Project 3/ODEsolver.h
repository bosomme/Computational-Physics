//
//  ODEsolver.h
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 21.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#ifndef ODEsolver_h
#define ODEsolver_h

#define _USE_MATH_DEFINES // Mathematical constants

#include "planet.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>


class ODEsolver
{
public:
    friend class planet;
    
    
    double G;
    double v_x_i; double v_y_i;
    double Kinetic_Energy; double Potential_Energy;
    double Total_Energy;
    
    // Initializers:
    ODEsolver();
    ODEsolver(double initial_velocity_x, double initial_velocity_y);
    
    // Functions:
    void Eulers_Method(int N, double h, double *x, double *y, double *t);
    void Velocity_Verlet(int N, double h,  double *x, double *y, double *t);
    
    void Energy(double mass_of_planet, double r);
    void print_energy();
};



#endif /* ODEsolver.h */
