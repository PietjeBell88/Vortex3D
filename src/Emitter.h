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

#pragma once

///////////
// Headers
#include <string.h>

#include <blitz/tinymat.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

#include "Typedefs.h"


////////////////////////
// Forward Declarations
class ParticleArray;
class Particle;


/////////////
// Namespace
using blitz::TinyVector;
using blitz::TinyMatrix;
using std::string;


///////////
// Emmiter
class Emitter
{
protected:
    // Local variables.
    TinyVector<int, 3> p_grid; // X x Y x Z grid of particles
    TinyMatrix<double, 3, 2> delimiter; // The edges of the box.
    double dx, dy, dz; // The stepsizes.

    // For the ones below, please see --help.
    double p_density;
    double p_diameter;
    double p_velocity;

    double p_rate;
    int reset_particles;
    int p_N;

    // Functions only used local in init(), update() and reset().
    virtual Vector3d startPos( int p ) = 0;

    virtual Vector3d startVel( int p ) = 0;

public:
    Emitter( double p_density, double p_diameter, double p_velocity, 
             const string &dimensions, double radius, double p_rate,
             int reset_particles );

    virtual ~Emitter();

    virtual void init( ParticleArray *particles ) = 0;

    virtual void update( double relative_time, ParticleArray *particles ) = 0;

    // reset():
    // In:     - (1) The number of the particle to be reset.
    //         - (2) The current relative time (as fraction of going around time T_l).
    // In/Out: - (3) The particlearray of which the particle is to be removed.
    // Out:    - (return) The relative time (fraction of going around time T_l) the particle spent in the box.
    virtual double reset( int p, double relative_time, ParticleArray *particles );
};
