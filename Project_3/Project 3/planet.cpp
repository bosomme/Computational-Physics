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
    x = 0.0;
    y = 0.0;
    v_x = 0.0;
    v_y = 0.0;
}

planet::planet(double mass, double pos_x, double pos_y, double vx, double vy){
    mass_of_sun = 1.0;
    mass_of_planet = mass;
    x = pos_x;
    y = pos_y;
    v_x = vx*365;
    v_y = vy*365;
}

double planet::distance(planet other){
    double X = x-other.x; double Y = y-other.y;
    return sqrt(X*X + Y*Y);
}

void planet::Gravitational_Force(planet &other, double G, double &Fx, double &Fy){
    double r = this->distance(other);
    double F = -G*mass_of_planet*other.mass_of_planet/(r*r);
    Fx += F*(x-other.x)/r; Fy += F*(y-other.y)/r;
}


double planet::Kinetic_Energy(double v_x, double v_y){
    double vv = v_x*v_x + v_y*v_y;
    return 0.5*mass_of_planet*vv;
}

double planet::Potential_Energy(double G, double r){
    return -G*mass_of_planet*mass_of_sun/r;
}


