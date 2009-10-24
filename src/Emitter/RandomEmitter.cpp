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
#include "RandomEmitter.h"

#include "..\..\external\MTRand.h"

#include "..\ParticleArray.h"
#include "..\Particle.h"


///////////////
// Constructor
RandomEmitter::RandomEmitter( const Vortex3dParam &param ) :
    Emitter( param )
{
    last_emit_time = 0;
}

//////////////
// Destructor
RandomEmitter::~RandomEmitter() {}


//////////////////////////////////////////
// Particle Property Generators (private)
Vector3d RandomEmitter::startPos( int p )
{
    /*
     * This function directly calculates the starting position from the
     * particle number; this comes in handy when resetting particles, if
     * they leave the cube, to the position they started;
     */
    MTRand myran;

    double x = delimiter(0, 0) + (delimiter(0, 1) - delimiter(0, 0)) * myran.rand53();
    double y = delimiter(1, 0) + (delimiter(1, 1) - delimiter(1, 0)) * myran.rand53();
    double z = delimiter(2, 0) + (delimiter(2, 1) - delimiter(2, 0)) * myran.rand53();

    return Vector3d( x, y, z );
}

Vector3d RandomEmitter::startVel( int p )
{
    return Vector3d( 0, 0, p_velocity );
}


///////////////////////////////////////////
// Init, Update (public), Reset is default
void RandomEmitter::init( ParticleArray *particles ) {}

void RandomEmitter::update( double relative_time, ParticleArray *particles )
{
    int to_emit = static_cast<int> ( floor( (relative_time - last_emit_time) * p_rate ) );

    if ( to_emit >= 1 )
        last_emit_time = relative_time;

    while ( (particles->getMaxLength() - particles->getLength()) >= 1
            && to_emit >= 1 )
    {
        particles->add( startPos( 0 ), startVel( 0 ), relative_time );
        to_emit--;
    }
}
