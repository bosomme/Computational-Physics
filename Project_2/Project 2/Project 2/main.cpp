/*
    main.cpp
    Project 2
 
    Created by Børge Olsen-Hagen Sømme on 29.09.2017.
    Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <time.h>

#include "lib.h"

using namespace std;



void initialize(int, double, double*, double**, double**);
void find_max_off_diag(int, double**, double*, int*, int*);



int main(int argc, const char * argv[])
{
    int n = 5;
    int counter = 1; int max_counts = n*n*n;
    int rho_min = 0; double rho_max = 1.0e1;
    
    double tolerance = 1.0e-8;
    
    double h = (rho_max - rho_min)/((double) n);
    
    double** a;
    a = (double**) matrix(n, n, sizeof(double));
    double** v;
    v = (double**) matrix(n, n, sizeof(double));
    double* rho = new double[n];
    
    // Fill matrices and rho:
    initialize(n, h, rho, a, v);
    
    
    double t = 0, s = 0, c = 0, tau = 0;
    
    int k = n-1, l = n-2;
    double amax = a[k][l];
    
    double akk = a[k][k], all = a[l][l], akl = amax;
    double aik = 0, ail = 0;
    
    
    while(fabs(amax) > tolerance){
        // update max value:
        if(counter > 1){
            amax = 0;
            find_max_off_diag(n, a, &amax, &k, &l);
        }
        
        tau = (a[l][l] - a[k][k])/(2*amax);
        if(tau>0){
            t = 1/(tau+sqrt(1+tau*tau));
        }
        else{
            t = -1/(-tau+sqrt(1+tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = t*c;
        
        // new matrix elements and vectors:
        for(int i=0; i<n; i++){
            if(i!=k && i!=l){
                aik = a[i][k]; ail = a[i][l];
                a[i][k] = aik*c - ail*s;
                a[k][i] = a[i][k];
                a[i][l] = ail*c + aik*s;
                a[l][i] = a[i][l];
            }
        }
        a[k][k] = akk*c*c - 2*akl*c*s + all*s*s;
        a[l][l] = all*c*c + 2*akl*c*s + akk*s*s;
        a[k][l] = (akk-all)*c*s + akl*(c*c - s*s);
        
        counter ++;
        
    }
    
    // print eigenvalues:
    cout << "The eigenvalues are:" << endl;
    for(int i=0; i<n; i++){
        cout << a[i][i] << endl;
    }
    cout << "number of counts:" << counter;
    
    
    delete [] rho;
    free_matrix((void **) a);
    free_matrix((void **) v);
    
    
    return 0;
}



void initialize(int n, double h, double* rho, double** a, double** v){
    // initialize rho:
    rho[0] = 0;
    for(int i=0; i<n; i++){
        rho[i] = rho[i-1] + h;
    }
    
    // initialize matrix:
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                a[i][j] = 2/(h*h) + rho[i]*rho[i];
                v[i][j] = 1;
            }
            else if(i==j+1 or i==j-1){
                a[i][j] = -1/(h*h);
            }
            else{
                a[i][j] = 0;
                v[i][j] = 0;
            }
        }
    }
}


void find_max_off_diag(int n, double** a, double* amax, int* k, int* l){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i!=j && fabs(a[i][j])>=fabs(*amax)){
                *amax = a[i][j];
                *k = i;
                *l = j;
            }
        }
    }
}



















