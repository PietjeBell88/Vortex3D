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
#include <blitz/tinymat.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

#include "Typedefs.h"

using std::string;


////////////////////////
// Forward Declarations
class Emitter;
class Vortex;
class ParticleArray;


///////////////////
// Global Settings
struct Vortex3dParam {
    // Cmdline parameters
    int vortextype, reset_particles, maxparticles, grav, emittertype, outputtype, outputformat, outputintervalmethod;
    double radius, velocity, angle, p_density, p_diameter, fl_mu, fl_density,
           duration, p_velocity, dtscale, p_rate, outputinterval;
    bool interpolate, rotategrav;
    string parameters, roi, dimensions, datafile;

    // Calculated parameters
    int max_t;
    double fl_nu;
    double systemtime;
    Vector3d gravity;
    double dt;
    double beta; // This ratio is used in the equation of motion
    double tau_a;
    double terminal_velocity;
    double t_ref;

    // Emitter
    TGrid emitter_grid;
    TDelimiter emitter_delimiter;
    double emitter_dx, emitter_dy, emitter_dz;
    int p_N;

    // ROI
    TGrid roi_grid;
    TDelimiter roi_delimiter;
    double roi_dx, roi_dy, roi_dz;
};

/////////////
// Functions


void getConcentration( const ParticleArray &particles, const Vortex3dParam &param,
                       ScalarField *concentration );

void moveParticles( Vortex *the_vortex, Emitter *the_emitter,
                    ParticleArray *particles, const Vortex3dParam &param );
void checkParticles( Vortex *the_vortex, Emitter *the_emitter,
                     ParticleArray *particles, double relative_time,
                     double *average_fall_time, double *particles_out );

inline void writeProgress( int perc );

void show_help();

void parse( int argc, char* argv[], Vortex3dParam *param );

void readGridDelimiterDelta( const string &fstring, const double &radius, TGrid *grid, 
                             TDelimiter *delimiter, double *dx, double *dy, double *dz );

void printParam( const Vortex3dParam &param );

