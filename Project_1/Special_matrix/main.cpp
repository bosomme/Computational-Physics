/*
  Special_matrix
  Project 1

  Created by Børge Olsen-Hagen Sømme on 01.09.2017.
  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h> // Computing time

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
    //int n = 100000;
    
    /*
    cout << "Please give wanted output filename: "; // filename for printing to file
    string outfile_name;
    cin >> outfile_name;
    */
     
    double h = 1.0/(n+1);
    double *x = new double[n+2];
    
    double *u = new double[n+2];
    u[0] = 0;
    double *v = new double[n+2];
    v[0] = 0;
    double *f = new double[n+1];
    f[0] = 0;
    
    double *b_tilde = new double[n+1];
    double *f_tilde = new double[n+1];
    
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
    }
    
    b_tilde[1] = 2;
    f_tilde[1] = f[1];
    
    // Compute b_tilde in advance
    for (int i=2; i<=n; i++)
    {
        b_tilde[i] = float(i+1)/float(i);
    }
    
    
    // Forward substitution
    for (int i=2; i<=n; i++)
    {
        f_tilde[i] = f[i] + float(f_tilde[i-1])/b_tilde[i-1];
    }
    
    v[n] = float(f_tilde[n])/float(b_tilde[n]);
    
    // Backward substitution
    for (int i=n-1; i>=1; i--)
    {
        v[i] = (f_tilde[i] + v[i+1])/b_tilde[i];
    }
    
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC);
    
    //cout << setiosflags(ios::showpoint | ios::uppercase);
    //cout << setw(20) << setprecision(10) << "Time used for special matrix:" << timeused << endl;
    
    // Computing and extracting maximum relative error
    double Relative_error; double new_error;
    Relative_error = log10(fabs((v[1]-u[1])/u[1]));
    for (int i=1; i<=n; i++)
    {
        new_error = log10(fabs((v[i]-u[i])/u[i]));
        if (new_error > Relative_error)
            Relative_error = new_error;
    }
    
    cout << "Maximum relative error: " << Relative_error << endl;
    cout << "h^2 = " << h*h << endl;
    
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
    ofile << "Time used by algorithm = " << timeused;
    ofile.close();
    */
    
    delete [] f;
    delete [] x;
    delete [] b_tilde;
    delete [] f_tilde;
    delete [] v;
    delete [] u;
    
    
    
    return 0;
}

