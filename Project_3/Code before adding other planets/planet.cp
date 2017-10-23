//
//  planet.cpp
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 22.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#include "planet.h"

planet::planet(){
    mass_of_sun = 1.0;
    mass_of_planet = 0.000003;
}

planet::planet(double mass){
    mass_of_sun = 1.0;
    mass_of_planet = mass;
}

double planet::Kinetic_Energy(double v_x, double v_y){
    double vv = v_x*v_x + v_y*v_y;
    return 0.5*mass_of_planet*vv;
}

double planet::Potential_Energy(double G, double r){
    return -G*mass_of_planet*mass_of_sun/r;
}

