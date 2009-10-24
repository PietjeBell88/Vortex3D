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
#include "Output.h"

#include "../Typedefs.h"
#include "../Vortex/Vortex.h"
#include "../ParticleArray.h"
#include "../Particle.h"


//////////////
// Constructor
Output::Output( const Vortex3dParam &param, Vortex *the_vortex, FILE * f ) 
{
    this->outputtype = param.outputtype;

    this->grid = param.roi_grid;
    this->delimiter = param.roi_delimiter;
    this->dx = param.roi_dx;
    this->dy = param.roi_dy;
    this->dz = param.roi_dz;

    this->the_vortex = the_vortex;

    this->f = f;

    this->param = param;
}


//////////////
// Destructor
Output::~Output() {}


void Output::getConcentration( const ParticleArray &particles, ScalarField *concentration )
{
    *concentration = 0;

    // Readability. Hopefully optimized away.
    const int & length = particles.getLength();

    // Loop over every particle and increase the count by one for the box it's in.
    for ( int p = 0; p < length; p++ )
    {
        const Vector3d & pos = particles.getParticle(p).getPos();

        int i = static_cast<int> (floor((pos(0) - delimiter(0, 0)) / dx));
        int j = static_cast<int> (floor((pos(1) - delimiter(1, 0)) / dy));
        int k = static_cast<int> (floor((pos(2) - delimiter(2, 0)) / dz));

        ++(*concentration)(i, j, k);
    }

    if ( sum( *concentration ) != 0 )
        *concentration = (*concentration) / sum( *concentration );
}

//
void Output::printFileHeader() { 
}

//
void Output::printFileFooter() { 
}


//
void Output::writeToFile( double time, const ParticleArray &particles )
{
    static bool first_call = true;

    switch ( outputtype ) {
        case 1:
            writeTrajectories( first_call, time, particles );
            break;
        case 2:
            writeConcentration( first_call, time, particles );
            break;
        case 3:
            writeVelocityField( first_call, time );
            break;
        default:
            std::cout << "ERROR: Unknown outputtype.";
    } 

    first_call = false;
}

