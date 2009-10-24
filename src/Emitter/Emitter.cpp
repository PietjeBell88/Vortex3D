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

#include "..\ParticleArray.h"
#include "..\Particle.h"


///////////////
// Constructor
Emitter::Emitter( const Vortex3dParam &param )
{
    // Save the variables.
    this->p_density = param.p_density;
    this->p_diameter = param.p_diameter;
    this->p_velocity = param.p_velocity;

    this->grid = param.emitter_grid;
    this->delimiter = param.emitter_delimiter;
    this->dx = param.emitter_dx;
    this->dx = param.emitter_dy;
    this->dx = param.emitter_dz;
    this->p_N = param.p_N;

    this->p_rate = param.p_rate;

    this->reset_particles = param.reset_particles;
}


//////////////
// Destructor
Emitter::~Emitter() {}


/////////
// Reset
double Emitter::reset( int p, double relative_time, ParticleArray *particles )
{
    // This is a default reset() function which probably will not have to be overrided.
    Particle removed = particles->remove( p );

    if ( reset_particles != 0 )
        particles->add( startPos( p ), startVel( p ), relative_time );

    return (relative_time - removed.spawnTime());
}
