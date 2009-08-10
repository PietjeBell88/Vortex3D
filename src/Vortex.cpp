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
#include "Vortex.h"


///////////////
// Constructor
Vortex::Vortex( const double &radius, const double &velocity, const double &angle,
                const double &fl_mu, const double &fl_density,
                const bool &interpolate, const string &roi )
{
    this->radius = radius;
    this->velocity = velocity;
    this->angle = angle;
    this->fl_mu = fl_mu;
    this->fl_density = fl_density;
    this->fl_nu = fl_mu / fl_density;
    this->interpolate = interpolate;

    double x1, x2, y1, y2, z1, z2;
    int X, Y, Z;

    sscanf( roi.c_str(), "[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",
            &x1, &X, &x2, &y1, &Y, &y2, &z1, &Z, &z2 ); //e.g. [-4:30:4,0:1:0,4:1:4]"

    delimiter = x1 * radius, x2 * radius,
                y1 * radius, y2 * radius,
                z1 * radius, z2 * radius;
    grid = X, Y, Z;

    // Calculate the stepsizes. Avoid dividing by zero.
    if ( grid(0) <= 1 )
        dx = 0;
    else
        dx = (delimiter(0, 1) - delimiter(0, 0)) / (grid(0) - 1);

    if ( grid(1) <= 1 )
        dy = 0;
    else
        dy = (delimiter(1, 1) - delimiter(1, 0)) / (grid(1) - 1);

    if ( grid(2) <= 1 )
        dz = 0;
    else
        dz = (delimiter(2, 1) - delimiter(2, 0)) / (grid(2) - 1);
}


//////////////
// Destructor
Vortex::~Vortex() {}


////////////////////////////
// Initialize Interpolation
// Can't do interpolation initialization in the constructor, results in pure function call.
// As MSDN says: "Do not call overridable methods in constructors".
void Vortex::initInterpolate()
{
    v.resize( grid(0), grid(1), grid(2) );
    accelfluid.resize( grid(0), grid(1), grid(2) );
    SetupVortexGrid( &v, &accelfluid );
}


////////////////////
// Public Functions
bool Vortex::outsideBox( const Vector3d &pos )
{
    if ( pos(0) <= delimiter(0, 0) || pos(0) > delimiter(0, 1) ||
         pos(1) <= delimiter(1, 0) || pos(1) > delimiter(1, 1) ||
         pos(2) <= delimiter(2, 0) || pos(2) > delimiter(2, 1) )
    {
        return true;
    }
    return false;
}

Vector3d Vortex::getDuDtAt( const Vector3d &pos )
{
    if ( interpolate == true )
        return Interpolate3DCube( accelfluid, pos );
    else
        return dudtAngle( pos );
}

Vector3d Vortex::getVelocityAt( const Vector3d &pos )
{
    if ( interpolate == true )
        return Interpolate3DCube( v, pos );
    else
        return velocityAngle( pos );
}

VectorField Vortex::getVectorField()
{
    return v;
}


///////////////////////////////////////////////////////////////////
// Transformation matrix from [v_r, v_phi, v_z] to [v_x, v_y, v_z]
TinyMatrix<double, 3, 3> Vortex::Cil2Cart( const double &phi )
{
    TinyMatrix<double, 3, 3> M;

    M = cos(phi), -sin(phi), 0,
        sin(phi), cos(phi),  0,
        0,        0,         1;

    return M;
}


//////////////////////////////////////////////////////////////////////////////////////
// Rotation matrix for rotation around the x axis ("folding the y-axis to the z-axis")
TinyMatrix<double, 3, 3> Vortex::Rotate_x( const double &angle )
{
    TinyMatrix<double, 3, 3> M;

    M = 1, 0,          0,
        0, cos(angle), -sin(angle),
        0, sin(angle), cos(angle);

    return M;
}


///////////////////
// Interpolation
Vector3d Vortex::Interpolate3DCube( const VectorField &v, const Vector3d &pos )
{
    /*
     This function interpolates the velocity in matrix U to the position of a particle at P.
     It does this by interpolating linearly by assigning weights to each corner of the surrounding box.
     This function requires uniform gridspacing.
     */

    /*
     Start by finding the lower index values of the box surrounding the particle.
     Of course, this step requires that the size of the index does not exceed the integer gridRange.
     */

    int i = static_cast<int> (floor((pos(0) - delimiter(0, 0)) / dx));
    int j = static_cast<int> (floor((pos(1) - delimiter(1, 0)) / dy));
    int k = static_cast<int> (floor((pos(2) - delimiter(2, 0)) / dz));

    /*
     * Calculate the weighting factors for each corner of the box
     * Note: 0 <= x < 1 (and same for y and z)
     * Because the position of the particle can be negative, make it positive first
     * and then do a modulo. Doing modulo on a negavite value can be confusing and inconsistent.
     */

    double x = fmod(pos(0) - delimiter(0, 0), dx) / dx;
    double y = fmod(pos(1) - delimiter(1, 0), dy) / dy;
    double z = fmod(pos(2) - delimiter(2, 0), dz) / dz;

    //Do a weighted addition of all the corners of the cube surrounding the particle.
    //Please note that v(i,j,k) is indeed a Vector3d
    return v(i, j, k) * (1 - x) * (1 - y) * (1 - z) + v(i + 1, j, k) * x * (1
            - y) * (1 - z) + v(i, j + 1, k) * (1 - x) * y * (1 - z) + v(i, j, k
            + 1) * (1 - x) * (1 - y) * z + v(i + 1, j, k + 1) * x * (1 - y) * z
            + v(i, j + 1, k + 1) * (1 - x) * y * z + v(i + 1, j + 1, k) * x * y
            * (1 - z) + v(i + 1, j + 1, k + 1) * x * y * z;
}


////////////////////////////////////////
// Vortex Velocity and Du/Dt Getters
Vector3d Vortex::velocityCarthesian( const Vector3d &pos )
{
    const double &x = pos(0);
    const double &y = pos(1);
    const double &z = pos(2);

    double r = sqrt( x * x + y * y );
    double phi = atan2( y, x );

    return product( Cil2Cart( phi ), velocityCylinder( r, phi, z ) );
}

Vector3d Vortex::dudtCarthesian( const Vector3d &pos )
{
    const double &x = pos(0);
    const double &y = pos(1);
    const double &z = pos(2);

    double r = sqrt( x * x + y * y );
    double phi = atan2( y, x );

    return product( Cil2Cart( phi ), dudtCylinder( r, phi, z ) );
}

Vector3d Vortex::velocityAngle( const Vector3d &pos )
{
    return product( Rotate_x( angle ), velocityCarthesian( product( Rotate_x( -angle ), pos ) ) );
}

Vector3d Vortex::dudtAngle( const Vector3d &pos )
{
    return product( Rotate_x( angle ), dudtCarthesian( product( Rotate_x( -angle ), pos ) ) );
}


//////////////////////////
// Initialize Vortex Grid
void Vortex::SetupVortexGrid(VectorField *v, VectorField *accelfluid)
{
#pragma omp parallel for
    for ( int i = 0; i < grid(0); i++ )
    {
        // What x-coordinate are we at?
        double x = delimiter(0, 0) + i * dx;
        for ( int j = 0; j < grid(1); j++ )
        {
            // What y-coordinate are we at?
            double y = delimiter(1, 0) + j * dy;
            for ( int k = 0; k < grid(2); k++ )
            {
                // What z-coordinate are we at?
                double z = delimiter(2, 0) + k * dz;

                // Calculate and set the velocities for each direction in the VectorField
                // Of course you can't have a negative index, so save each at their index+1
                (*v)(i, j, k) = velocityAngle( Vector3d( x, y, z ) );
                (*accelfluid)(i, j, k) = dudtAngle( Vector3d( x, y, z ) );
            }
        }
    }
}
