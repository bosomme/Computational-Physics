/*
  General_matrix.cpp
  Project 1

  Created by Børge Olsen-Hagen Sømme on 06.09.2017.
  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.

  Program used for solving general matrix, and print results to file
  for plotting with python
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

ofstream ofile;


double func(double x) {return 100*exp(-10*x);}

double solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}


int main(int argc, const char * argv[])
{
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = f(i)
    // u(i) is analytical solution
    
    cout << "Please give n: ";
    int n;
    cin >> n;
    
    /*
    cout << "Please give wanted output filename: "; // filename for printing to file
    string outfile_name;
    cin >> outfile_name;
    */
    
    double h = 1.0/(n+1);
    double *x = new double[n+2];
    
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];
    
    double *u = new double[n+2];
    u[0] = 0;
    double *v = new double[n+2];
    v[0] = 0;
    double *f = new double[n+1];
    f[0] = 0;
    
    double *b_tilde = new double[n+1];
    double *f_tilde = new double[n+1];
    
    double temp_value;   // Value used in Gaussian elimination (expression for b_tilde)
    double f_temp;
    
    clock_t start, finish;  //  declare start and final time
    start = clock();
    
    for (int i=0; i<=n+1; i++)
    {
        x[i] = i*h;
    }
    
    // constructing right hand side:
    
    for (int i=1; i<=n; i++)
    {
        f[i] = h*h*func(x[i]);
        
        u[i] = solution(x[i]);
        
        // Also fill up a, b and c:
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    a[0] = 0;
    c[n] = 0;
    

    b_tilde[1] = b[1];
    f_tilde[1] = f[1];
    v[1] = f[1]/b[1];

    
    // Forward substitution
    for (int i=2; i<=n; i++)
    {
        temp_value = c[i-1]/b_tilde[i-1];
        f_temp = f_tilde[i-1]/b_tilde[i-1];
        
        b_tilde[i] = b[i] - a[i]*temp_value;
        
        f_tilde[i] = f[i] - a[i]*f_temp;
    }
    
    // Backward substitution
    for (int i=n-1; i>=1; i--)
    {
        v[i] = (f_tilde[i] - c[i]*v[i+1])/b_tilde[i];
    }
    
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC);
    
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setw(20) << setprecision(10) << "Time used:" << timeused << endl;
    
    /*
    // Print to file
 
    ofile.open(outfile_name);
    
    ofile << "Results of program main.cpp, numerical - analytical" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0; i<=n+1; i++)
    {
        ofile << setprecision(8) << x[i] << "  ";
        ofile << setprecision(8) << v[i] << "  ";
        ofile << setprecision(8) << u[i] << endl;
    }
    */
    
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
    delete [] x;
    delete [] b_tilde;
    delete [] f_tilde;
    delete [] v;
    delete [] u;
    
    ofile.close();
    
    return 0;
}

