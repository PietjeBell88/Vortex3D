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
#include "ParticleArray.h"
#include "Vortex.h"
#include "Emitter.h"

using blitz::TinyVector;
using blitz::TinyMatrix;
using std::string;


///////////////////
// Global Settings
struct Settings {
    double systemtime;
    double gravity;
    double dt;
    int reset_particles;
    double beta;                                    // This ratio is used in the equation of motion
    double tau_a;
    string datafile;
    int outputtype;
    // Used for concentration
    TinyVector<int,3> grid;    // X x Y x Z grid of resolution
    TinyMatrix<double,3,2> delimiter; //offsets
    bool tecplot;
    double dx, dy, dz;
};

/////////////
// Functions
void ReadROI(const string &roi, const double &radius, Settings *options);

void GetConcentration(ParticleArray *particles, const Settings &options, ScalarField *concentration);

void MoveParticles(Vortex *theVortex, Emitter *theEmitter, ParticleArray *particles, const Settings &options);
void CheckParticles(Vortex *theVortex, Emitter *theEmitter, ParticleArray *particles, const double &relative_time, double *average_fall_time, double *particles_out);

inline void WriteProgress(int perc);

inline void WriteToFile(const double time, const Settings &options, FILE * f, ParticleArray *particles, Vortex *thevortex);

