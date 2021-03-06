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
#include <fstream>
#include "..\Typedefs.h"
#include "..\Main.h" 


////////////////////////
// Forward Declarations
class Vortex;
class ParticleArray;


///////////////
// Output
class Output
{
protected:
    /////////////
    // Variables
    int outputtype;        // Particle positions, concentration, velocity field etc..  (See help)

    TGrid grid;            // X x Y x Z grid of resolution
    TDelimiter delimiter;  // Offsets
    double dx, dy, dz;     // deltas
    double timestep;       // Output interval in relative time.

    Vortex *the_vortex;    // The Vortex (used when emitting velocity field)

    FILE *f;

    Vortex3dParam param;   // To output the cmdline options to file when wanted.


    /////////////
    // Functions
    void getConcentration( const ParticleArray &particles, ScalarField *concentration );

    virtual void printFileHeader();
    virtual void printFileFooter();

    virtual void writeConcentration( bool first_call, double time, const ParticleArray &particles ) = 0;
    virtual void writeTrajectories( bool first_call, double time, const ParticleArray &particles ) = 0;
    virtual void writeVelocityField( bool first_call, double time ) = 0;

public:
    Output( const Vortex3dParam &param, Vortex *the_vortex );

    virtual ~Output();

    void writeToFile( double time, const ParticleArray &particles );
};
