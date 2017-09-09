//
//  main.cpp
//  LU_decomposition
//
//  Created by Børge Olsen-Hagen Sømme on 07.09.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>

#include "lib.h"

using namespace std;

ofstream ofile;

double func(double x) {return 100*exp(-10*x);}

double solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}


int main(int argc, const char * argv[])
{
    // LU-decomposition: A = LU
    // Solve exuation Ax = b
    
    //cout << "Please give n: ";
    //int n;
    //cin >> n;
    int n = 10000;
    
    /*
    cout << "Please give wanted output filename: "; // filename for printing to file
    string outfile_name;
    cin >> outfile_name;
    */

    double h = 1.0/n;
    double *x = new double[n+1];
    double *b = new double[n+1]; double *u = new double[n+1];
    b[0] = 0; u[0] = 0;
    
    // compose tridiagonal matrix:
    double **A;
    A = (double **) matrix(n, n, sizeof(double));
    A[0][0] = 2;
    
    clock_t start, finish;  //  declare start and final time
    start = clock();
    
    for(int i = 1; i < n; i++)
    {
        for(int j = 1; j < n; j++)
        {
            A[i][j] = 0;
        }
        A[i][i] = 2;
        A[i][i-1] = -1;
        A[i-1][i] = -1;
    }
    
    for (int i=0; i<=n; i++)
    {
        x[i] = i*h;
    }
    
    // constructing right hand side:
    
    for (int i=0; i<=n; i++)
    {
        b[i] = h*h*func(x[i]);
        
        u[i] = solution(x[i]);
    }
    
    int *indx = new int[n]; double d;

    
    ludcmp(A, n, indx, &d);
    
    lubksb(A, n, indx, b);
    
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC);
    
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setw(20) << setprecision(10) << "Time used for LU-decomposition:" << timeused << endl;
    
    /*
    // Print to file
    
    ofile.open(outfile_name);
    
    ofile << "Results of program main.cpp, numerical - analytical" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0; i<=n; i++)
    {
        ofile << setprecision(8) << x[i] << "  ";
        ofile << setprecision(8) << b[i+1] << "  ";
        ofile << setprecision(8) << u[i] << endl;
    }
    ofile.close();
    */

    delete [] x;
    delete [] b;
    delete [] u;
    delete [] indx;
    free_matrix((void **) A);

    
    return 0;
}
