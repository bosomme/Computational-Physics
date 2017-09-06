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
    
    cout << "Please give n: ";
    int n;
    cin >> n;
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
    
    double *b_tilde = new double[n+1];
    double *f_tilde = new double[n+1];
    b_tilde[1] = b[1];
    f_tilde[1] = f[1];
    v[1] = f[1]/b[1];
    
    double temp_value;   // Value used in Gaussian elimination (expression for b_tilde)
    double f_temp;
    
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

