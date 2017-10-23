//
//  planet.h
//  Project 3
//
//  Created by Børge Olsen-Hagen Sømme on 22.10.2017.
//  Copyright © 2017 Børge Olsen-Hagen Sømme. All rights reserved.
//

#ifndef planet_h
#define planet_h

#define _USE_MATH_DEFINES

#include <cmath>

class planet
{
public:
    double mass_of_sun;
    double mass_of_planet;
    
    // Initializers
    planet();
    planet(double mass);
    
    // Functions
    double Kinetic_Energy(double velocity_x, double velocity_y);
    double Potential_Energy(double G, double r);
    
};

#endif /* planet_h */
