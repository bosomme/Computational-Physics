//
//  main.cpp
//  Project 1
//
//  Created by Børge Olsen-Hagen Sømme on 01.09.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

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
    double h = 1.0/(n+1);
    double *x = new double[n+2];
    
    double *u = new double[n+2];
    u[0] = 0;
    double *v = new double[n+2];
    v[0] = 0;
    double *f = new double[n+1];
    f[0] = 0;
    
    for (int i=0; i<=n+1; i++)
    {
        x[i] = i*h;
    }
    
    // constructing right hand side:
    
    for (int i=1; i<=n; i++)
    {
        f[i] = h*h*func(x[i]);
        
        u[i] = solution(x[i]);
    }
    
    double *b_tilde = new double[n+1];
    double *f_tilde = new double[n+1];
    b_tilde[1] = 2;
    f_tilde[1] = f[1];
    
    //double temp_value = 2;   // Value used in Gaussian elimination (expression for b_tilde)
    
    // Forward substitution
    for (int i=2; i<=n; i++)
    {
        
        b_tilde[i] = float(i+1)/float(i);
        //temp_value = 2 - float(1)/float(temp_value);
        f_tilde[i] = f[i] + float((i-1)*f_tilde[i-1])/float(i);
    }
    
    v[n] = float(f_tilde[n])/float(b_tilde[n]);
    
    // Backward substitution
    for (int i=n-1; i>=1; i--)
    {
        v[i] = float(i)/float(i+1)*(f_tilde[i]-v[i+1]);
    }
    
    
    // Print to file
    cout << "Please give wanted output filename: ";
    string outfile_name;
    cin >> outfile_name;
    
    ofile.open(outfile_name);
    
    ofile << "Results of program main.cpp, numerical - analytical" << endl;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0; i<=n+1; i++)
    {
        ofile << setprecision(8) << x[i] << "  ";
        ofile << setprecision(8) << v[i] << "  ";
        ofile << setprecision(8) << u[i] << endl;
    }
    
    delete [] f;
    delete [] x;
    delete [] b_tilde;
    delete [] f_tilde;
    delete [] v;
    delete [] u;
    
    ofile.close();
    
    return 0;
}

