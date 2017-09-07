//
//  main.cpp
//  LU_decomposition
//
//  Created by Børge Olsen-Hagen Sømme on 07.09.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    mat A = randu<mat>(5,5);
    mat B = randu<mat>(5,5);
    
    cout << A*B << endl;
    
    return 0;
}
