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


/////////////
// Namespace
using blitz::TinyVector;
using blitz::TinyMatrix;
using std::string;


////////////////
// Class Vortex
class Vortex
{
protected:
    /////////////
    // Variables
    double angle;       // Angle of the vortex with the z-axis, clockwise when looking along the y-axis in a rh coordinate system.
    double radius;      // Typical radius of the vortex
    double velocity;    // Velocity at radius

    double fl_mu;       // Fluid Dynamic Viscocity (Pa s).
    double fl_density;  // Fluid Density (kg/m3).
    double fl_nu;       // Fluid Kinematic Viscosity (m2/s).

    bool interpolate;   // true = Interpolate a precalculated grid of the fluid velocity at the position of a particle, false = evaluate fluid velocity functions at particle position


    ///////////////////////////////////////////////
    // Variables (Only Used when interpolate=true)
    TinyVector<int, 3> grid;             // X x Y x Z grid of particles
    TinyMatrix<double, 3, 2> delimiter;  //offsets

    double dx, dy, dz;                   // Stepsizes of the grid.

    VectorField v;
    VectorField accelfluid;


    /////////////
    // Functions

    // Transformation Matrices
    TinyMatrix<double, 3, 3> Cil2Cart( const double &phi );
    TinyMatrix<double, 3, 3> Rotate_x( const double &angle );

    // Interpolation
    Vector3d Interpolate3DCube( const VectorField &v, const Vector3d &pos );

    // Vortex Velocity and Du/Dt Getters
    virtual Vector3d velocityCylinder( const double &r, const double &phi, const double &z ) = 0;
    virtual Vector3d dudtCylinder( const double &r, const double &phi, const double &z ) = 0;

    virtual Vector3d velocityCarthesian( const Vector3d &pos );
    virtual Vector3d dudtCarthesian( const Vector3d &pos );

    virtual Vector3d velocityAngle( const Vector3d &pos );
    virtual Vector3d dudtAngle( const Vector3d &pos );

    // Vortex Grid Initialization (when interpolating)
    void SetupVortexGrid( VectorField *v, VectorField *accelfluid );

public:
    Vortex( const double &radius, const double &velocity, const double &angle,
            const double &fl_mu, const double &fl_density,
            const bool &interpolate, const string &roi );

    virtual ~Vortex();

    void initInterpolate();
    bool outsideBox( const Vector3d &pos );

    // And just 3 "important" public functions. The rest is private :>
    Vector3d getDuDtAt( const Vector3d &pos );
    Vector3d getVelocityAt( const Vector3d &pos );
    VectorField getVectorField();
};
