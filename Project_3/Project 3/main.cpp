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



using namespace std;

ofstream ofile;


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
        ODEsolver solve;
        solve.Eulers_Method(N, h, x, y, t);
        outfile_name = "eulers_method.txt";
    }
    else if (method == 1){
        ODEsolver solve(v_x_0, v_y_0);
        solve.Velocity_Verlet(N, h, x, y, t);
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














