//
//  ODEsolver.hpp
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 21.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#ifndef ODEsolver_h
#define ODEsolver_h

#define _USE_MATH_DEFINES // Mathematical constants

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>

using namespace std;

class ODEsolver
{
public:
    double G;
    double v_x_i; double v_y_i;
    
    // Initializers:
    ODEsolver();
    ODEsolver(double initial_velocity_x, double initial_velocity_y);
    
    // Functions:
    void Eulers_Method(int N, double h, double *x, double *y, double *t);
    void Velocity_Verlet(int N, double h,  double *x, double *y, double *t);
    
};



#endif //ODEsolver.h
