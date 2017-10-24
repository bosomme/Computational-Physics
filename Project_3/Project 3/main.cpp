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

    // Choose to observe the entire solar system, or studying the perihelion precession of Mercury
    // perihelion = 1 means that we study Mercury.
    int perihelion = 0;
    
    if (perihelion == 0){
        int N = 1E4;
        double Final_time = 50.;
        double h = Final_time/((double) N);
        
        // Initalialize planets+sun:
        // Using numbers from Planets.txt
        // Position precision: 6 decimals

        // planet sun(1.0, 0.0, 0.0, 0.0, 0.0);
        planet sun(1.0, 2.202805E-03, 5.751760E-03, -5.254272E-06, 5.476856E-06);
        planet earth(0.000003, 8.849192E-01, 4.661114E-01, -8.245182E-03, 1.520145E-02);
        planet jupiter(0.0009546, -4.568782E+00, -2.945080E+00, 3.999444E-03, -5.983215E-03);
        
        planet mercury(0.0000002, -2.683268E-01, -3.621236E-01, 1.696680E-02, -1.534282E-02);
        planet venus(0.0000024, -6.731236E-01, 2.489269E-01, -6.953139E-03, -1.911807E-02);
        planet mars(0.0000003, -1.579421E+00, 5.244129E-01, -3.841610E-03, -1.209597E-02);
        planet saturn(0.0002859, -3.316577E-01, -1.005006E+01, 5.270171E-03, -2.022489E-04);
        planet uranus(0.0000427, 1.785256E+01, 8.823199E+00, -1.771305E-03, 3.342609E-03);
        planet neptune(0.0000515, 2.861745E+01, -8.809267E+00, 9.029314E-04, 3.019282E-03);
        
        
        
        ODEsolver solve;
        
        solve.Add_Planet(sun);
        solve.Add_Planet(mercury);
        solve.Add_Planet(venus);
        solve.Add_Planet(earth);
        solve.Add_Planet(mars);
        solve.Add_Planet(jupiter);
        solve.Add_Planet(saturn);
        solve.Add_Planet(uranus);
        solve.Add_Planet(neptune);
        
        
        solve.Velocity_Verlet(N, h, "Solar_system.txt");
    }
    
    else if (perihelion == 1){
        int N = 1E5;
        double Final_time = 100.;
        double h = Final_time/((double) N);
        
        
        planet sun(1, 0.0, 0.0, 0.0, 0.0);
        planet mercury(0.0000002, 0.3075, 0.0, 0.0, 12.44/365.0);
        
        
        ODEsolver solve;
        
        
        solve.Add_Planet(sun);
        solve.Add_Planet(mercury);
        
        solve.Velocity_Verlet(N, h, "Mercury_Perihelion.txt")
    }
    
    return 0;
}














