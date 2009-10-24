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
#include "TecplotOutput.h"

#include "../ParticleArray.h" 
#include "../Particle.h" 
#include "../Vortex/Vortex.h" 


TecplotOutput::TecplotOutput( FILE * f, const Vortex3dParam &param, Vortex *the_vortex ) : 
                          Output(f, param, the_vortex) 
{
}
 

//////////////
// Destructor
TecplotOutput::~TecplotOutput() {}
    

inline void TecplotOutput::writeTrajectories( bool first_call, double time, const ParticleArray &particles )
{
    if ( first_call )
        fprintf( f, "#T     N     X     Y     Z     Speed\n" );

    for ( int i = 0; i < particles.getLength(); i++ )
    {
        // Readability
        const int & num = particles.getParticle( i ).getNum();
        const Vector3d & pos = particles.getParticle( i ).getPos();
        const double & speed = particles.getParticle( i ).speed();

        // Write output to file
        fprintf( f, "%e     %d     %e     %e     %e     %e\n",
                 time, num, pos(0), pos(1), pos(2), speed);
    }
    fprintf( f, "\n" );
}

// Concentration
inline void TecplotOutput::writeConcentration( bool first_call, double time, const ParticleArray &particles )
{
    if ( first_call )
    {
        fprintf( f,
            "TITLE=\"Example: Multi-Zone XY Line Plot Wwith Variable Sharing\"\n \
            VARIABLES=\"X\" \"Y\" \"Z\" \"C\"\n \
            ZONE T = \"%f seconds\", I = %d, J = %d, K = %d, DATAPACKING = POINT\n",
            time, grid(0), grid(1), grid(2) );
        first_call = false;

        ScalarField concentration( grid(0), grid(1), grid(2) );
        getConcentration( particles, &concentration );

        for ( int i = 0; i < grid(0); i++ )
        {
            double x = delimiter(0, 0) + i * dx;
            for ( int j = 0; j < grid(1); j++ )
            {
                double y = delimiter(1, 0) + j * dy;
                for ( int k = 0; k < grid(2); k++ )
                {
                    double z = delimiter(2, 0) + k * dz;
                    fprintf(f, "%e     %e     %e     %e     %e\n",
                        time, x, y, z, concentration(i, j, k));
                }
            }
        }
        fprintf( f, "\n" );
    }

    else
    {
        fprintf( f,
            "ZONE T = \"%f seconds\", I = %d, J = %d, K = %d, DATAPACKING = POINT\n \
            VARSHARELIST=([1-3]=1)\n",
            time, grid(0), grid(1), grid(2) );

        ScalarField concentration( grid(0), grid(1), grid(2) );
        getConcentration( particles, &concentration );

        for ( int i = 0; i < grid(0); i++ )
        {
            for ( int j = 0; j < grid(1); j++ )
            {
                for ( int k = 0; k < grid(2); k++ )
                {
                    fprintf( f, "%e\n", concentration(i, j, k) );
                }
            }
        }
        fprintf( f, "\n" );
    }
}

// Velocity Field
inline void TecplotOutput::writeVelocityField( bool first_call, double time )
{
    fprintf( f,
             "TITLE=\"Simple Data File\"\nVARIABLES=\"X\" \"Y\" \"Z\" \"UX\" \"UY\" \"UZ\"\nZONE I = %d, J = %d, K = %d, DATAPACKING = POINT\n",
             grid(0), grid(1), grid(2) );

    const VectorField &v = the_vortex->getVectorField();

    for ( int i = 0; i < grid(0); i++ )
    {
        double x = delimiter(0, 0) + i * dx;
        for ( int j = 0; j < grid(1); j++ )
        {
            double y = delimiter(1, 0) + j * dy;
            for ( int k = 0; k < grid(2); k++ )
            {
                double z = delimiter(2, 0) + k * dz;
                const Vector3d &vel = v(i, j, k);
                fprintf( f, "%e     %e     %e     %e     %e     %e\n",
                    x, y, z, vel(0), vel(1), vel(2) );
            }
        }
    }
    fprintf( f, "\n" );
}

