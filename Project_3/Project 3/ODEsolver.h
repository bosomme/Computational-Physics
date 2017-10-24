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

#include <vector>


class ODEsolver
{
public:
    friend class planet;
    
    
    double G;
    int number_of_planets;
    std::vector<planet> all_planets;
    double Kinetic_Energy; double Potential_Energy;
    double Total_Energy;
    
    // Initializers:
    ODEsolver();
    
    // Functions:
    void Add_Planet(planet new_planet);
    
    void Eulers_Method(int N, double h, double *x, double *y, double *t);
    void Velocity_Verlet(int N, double h, std::string outfile_name);
    
    void Gravitational_Force(planet &current, planet &other, double &Fx, double &Fy);
    
    void Energy(planet planet, double x, double y, double vx, double vy);
    void print_energy();
};



#endif /* ODEsolver.h */
