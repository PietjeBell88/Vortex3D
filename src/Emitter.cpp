// Copyright (c) 2009, Pietje Bell <pietjebell@ana-chan.com>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Pietje Bell Group nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


///////////
// Headers
#include "Emitter.h"

#include "ParticleArray.h"
#include "Particle.h"


///////////////
// Constructor
Emitter::Emitter(const double &p_density, const double &p_diameter, const double &p_velocity, 
                 const string &dimensions, const double &radius, const double &p_rate, const int &reset_particles) 
{    
    // Save the variables.
    this->p_density = p_density;
    this->p_diameter = p_diameter;
    this->p_velocity = p_velocity;
    this->reset_particles = reset_particles;
    this->p_rate = p_rate;

    // Calculate the delimiter, grid size (p_grid) and the stepsize (dx, dy and dz).
    double x1,x2,y1,y2,z1,z2;
    int X,Y,Z;
    sscanf(dimensions.c_str(),"[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",
        &x1,&X,&x2,&y1,&Y,&y2,&z1,&Z,&z2); //e.g. [-4:30:4,0:1:0,4:1:4]"
    delimiter = x1*radius, x2*radius,
        y1*radius, y2*radius,
        z1*radius, z2*radius;

    // Calculate the grid size.
    p_grid = X,Y,Z;
    p_N = product(p_grid);

    // Calculate the stepsizes. Avoid dividing by zero.
    if (p_grid(0) <= 1) {
        dx = 0;
    } 
    else {
        dx = (delimiter(0,1)-delimiter(0,0))/(p_grid(0)-1);
    }
    if (p_grid(1) <= 1) {
        dy = 0;
    } 
    else {
        dy = (delimiter(1,1)-delimiter(1,0))/(p_grid(1)-1);
    }
    if (p_grid(2) <= 1) {
        dz = 0;
    } 
    else {
        dz = (delimiter(2,1)-delimiter(2,0))/(p_grid(2)-1);
    }
}


/////////
// Reset
double Emitter::reset(const int &p, const double &relative_time, ParticleArray *particles) 
{
    // This is a default reset() function which probably will not have to be overrided.
    Particle removed = particles->remove(p); 

    if (reset_particles != 0) {
        particles->add(startPos(p), startVel(p), relative_time);
    }

    return (relative_time - removed.spawnTime());
}

