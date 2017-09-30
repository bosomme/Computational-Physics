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



void initialize(int, double, double*, double**, double**, int, double);
void find_max_off_diag(int, double**, double*, int*, int*);
void Jacobi(int, double, double**, double**);
void min_eigenvalue(int, double**, int, double*);


void unit_tests();


int main(int argc, const char *argv[])
{
    cout.precision(5);
    ///*
    int n = 300;
    double rho_max = 4.0;
    double h = rho_max/(n);
    
    double tolerance = 1.0e-8;
    
    
    double **a;
    a = (double**) matrix(n, n, sizeof(double));
    double **v;
    v = (double**) matrix(n, n, sizeof(double));
    double *rho = new double[n];
    
    int interaction = 1; //0 equals no interaction, 1 equals interaction
    double w_r = 0.01;
    
    // Fill matrices and rho:
    initialize(n, h, rho, a, v, interaction, w_r);
    
    Jacobi(n, tolerance, a, v);
    
    
    int N = 3;
    double *Eigenvalues = new double[N];
    min_eigenvalue(n, a, N, Eigenvalues);
    
    // print the N minimum eigenvalues:
    cout << "The " << N << " smallest eigenvalues are:" << endl;
    for(int i=0; i<N; i++){
        cout << Eigenvalues[i] << endl;
    }
    //cout << "number of counts: " << counter << endl;
    
    
    free_matrix((void **) a);
    free_matrix((void **) v);
    delete [] rho;
    delete [] Eigenvalues;
    //*/
    
    //unit_tests();       //Running unit tests, including conservation of orthonormality and max nondiagonal value of matrix.
    
    return 0;
}




// Function to initialize matrices and rho
// Takes interaction argument to decide if what potensial to use
void initialize(int n, double h, double* rho, double** a, double** v, int interaction, double w_r){
    // initialize rho:
    double hh = h*h;
    rho[0] = h;
    for(int i=1; i<n; i++){
        rho[i] = rho[i-1] + h;
    }
    
    // initialize matrix:
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                if(interaction==0){
                    a[i][j] = 2/hh + rho[i]*rho[i];
                }
                else if(interaction==1){
                    a[i][j] = 2/hh + w_r*rho[i]*rho[i] + 1/rho[i];
                }
                v[i][j] = 1;
            }
            else if(i==j+1 or i==j-1){
                a[i][j] = -1/hh;
            }
            else{
                a[i][j] = 0;
                v[i][j] = 0;
            }
        }
    }
}



// Function that finds the maximum value of the off-diagonal elements of a matrix
// Return the value as amax, and indices of the element as k and l
void find_max_off_diag(int n, double** a, double* amax, int* k, int* l){
    *amax = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i!=j && abs(a[i][j])>=abs(*amax)){
                *amax = a[i][j];
                *k = i;
                *l = j;
            }
        }
    }
}


// Function to implement the Jacobi algorithm, do the rotations, and find eigenvalues and eigenvectors
// Returns eigenvalues as diagonal of matrix a, and eigenvectors as rows in matrix v
void Jacobi(int n, double tolerance, double** a, double** v){
    int counter = 1;
    
    double t = 0, s = 0, c = 0, tau = 0;
    
    int k, l;
    double amax;
    find_max_off_diag(n, a, &amax, &k, &l);
    
    double akk, all;
    double aik = 0, ail = 0;
    
    double vik = 0, vil = 0;
    
    
    while(abs(amax) > tolerance){
        // update max value:
        find_max_off_diag(n, a, &amax, &k, &l);
        akk = a[k][k]; all = a[l][l];
        
        tau = (all - akk)/(2*amax);
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
            
            // Eigenvectors:
            vik = v[i][k]; vil = v[i][l];
            v[i][k] = vik*c - vil*s;
            v[i][l] = vil*c + vik*s;
            
        }
        
        a[k][k] = akk*c*c - 2*amax*c*s + all*s*s;
        a[l][l] = all*c*c + 2*amax*c*s + akk*s*s;
        a[k][l] = 0.0;
        a[l][k] = 0.0;
        
        
        counter ++;
    }
}


// Function to find the N minimum eigenvalues
// Returns the eigenvalues in array Ev
void min_eigenvalue(int n, double** a, int N, double* Ev){
    //N is how many eigenvalues we would like
    double max_value = 0;
    int g = 0;
    
    for(int i=0; i<N; i++){
        Ev[i] = a[i][i];
    }
    
    for(int i=0; i<N; i++){
        if(Ev[i] > max_value){
            max_value = Ev[i];
            g = i;
        }
    }
    for(int i=N; i<n; i++){
        if(a[i][i] <= max_value){
            Ev[g] = a[i][i];
            max_value = 0;
            for(int i=0; i<N; i++){
                if(Ev[i] > max_value){
                    max_value = Ev[i];
                    g = i;
                }
            }
        }
    }
}



/* 
    Testing
*/
int test_orthogonality();
int test_max_value_of_matrix();
int test_min_value();

// Function to implement all unit tests
void unit_tests(){
    if(test_orthogonality()==0 && test_max_value_of_matrix()==0 && test_min_value()==0){
        cout << "All good!" << endl;
    }
    else if(test_orthogonality() == 1){
        cout << "Problem with Jacobi" << endl;
        cout << "Function returns eigenvectors which are not orthogonal" << endl;
    }
    else if(test_max_value_of_matrix() == 2){
        cout << "Problem with find_max_off_diag" << endl;
        cout << "Function not returning max nondiagonal value of matrix" << endl;
    }
    else if(test_min_value() == 3){
        cout << "Problem with min_eigenvalue" << endl;
        cout << "Function not finding smallest eigenvalues" << endl;
    }
}

// Test to check if Jacobi-function returns orthonormal eigenvectors
// (That the rotations conserves the orthonormality)
int test_orthogonality(){
    int n = 5;
    double rho_max = 4.0;
    double h = rho_max/(n);
    
    double tolerance = 1.0e-8;
    
    int interaction = 0; double w_r = 0.0;      //No interaction
    
    double **a;
    a = (double**) matrix(n, n, sizeof(double));
    double **v;
    v = (double**) matrix(n, n, sizeof(double));
    double *rho = new double[n];
    
    initialize(n, h, rho, a, v, interaction, w_r);
    
    Jacobi(n, tolerance, a, v);
    
    
    double delta_ieqj = 0;double delta_inoteqj = 0;
    double test_tolerance = 1.0e-2;
    
    for(int i=0; i<n; i++){
        delta_ieqj += v[0][i]*v[0][i];
        delta_inoteqj += v[0][i]*v[1][i];
    }
    
    
    free_matrix((void **) a);
    free_matrix((void **) v);
    delete [] rho;
    
    
    if(abs(1-delta_ieqj) > test_tolerance){
        return 1;
    }
    else if(abs(delta_inoteqj) > test_tolerance){
        return 1;
    }
    else{
        return 0;
    }
}

// Test to check if the function find_max_off_diag returns the correct off-diagonal maximum
int test_max_value_of_matrix(){
    double **T_M;
    T_M = (double**) matrix(5, 5, sizeof(double));
    for(int i=0; i<5; i++){
        for(int j=0; j<5; j++){
            T_M[i][j] = i + 2*j;
        }
    }
    // Max-value of entire matrix is 20 on T_M[4][4]
    // But off-diagonal max is 19 on T_M[4][3]

    double max_value; int k, l;
    find_max_off_diag(5, T_M, &max_value, &k, &l);
    
    free_matrix((void **) T_M);
    
    if(k == 3 && l == 4 && max_value == 11.0){
        return 0;
    }
    else{
        return 2;
    }
}

// Test to check if the function min_eigenvalue returns the N minimum values in a diagonal matrix
int test_min_value(){
    double diag_values[5] = {5,3,7,6,2};
    double *values = new double[5];
    double **Matr = (double**) matrix(5, 5, sizeof(double));
    
    for(int i=0; i<5; i++){
        Matr[i][i] = diag_values[i];
    }
    
    min_eigenvalue(5, Matr, 3, values);
    if(values[0]==5 && values[1]==3 && values[2]==2){
        return 0;
    }
    else{
        return 3;
    }
}



















