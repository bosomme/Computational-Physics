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



#include "ODEsolver.h"
#include "planet.h"



using namespace std;


int main(int argc, const char * argv[]) {
    int N = 1E6;
    double Final_time = 50.;
    double h = Final_time/N;
    
    // Initalialize planets+sun:
    // Using numbers from Planets.txt
    // Position precision: 6 decimals
    planet sun(1.0, 0.0, 0.0, 0.0, 0.0);
    planet earth(0.000003, 8.849192E-01, 4.661114E-01, -8.245182E-03, 1.520145E-02);
    planet jupiter(0.0009546, -4.568782E+00, -2.945080E+00, 3.999444E-03, -5.983215E-03);
    
    
    ODEsolver solve;
    solve.Add_Planet(earth);
    solve.Add_Planet(jupiter);
    solve.Add_Planet(sun);
    
    
    solve.Velocity_Verlet(N, h, "Earth_Jupiter.txt");

    
    return 0;
}














