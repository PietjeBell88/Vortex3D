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
#include "BurgersVortex.h"


/////////////
// Namespace
using std::cout;


///////////////
// Constructor
BurgersVortex::BurgersVortex( const string &parameters, double radius, 
                              double velocity, double angle, double fl_mu, 
                              double fl_density, bool interpolate, const string &roi ) :
    Vortex( radius, velocity, angle, fl_mu, fl_density, interpolate, roi )
{
    //Evaluated LambartW function for alpha
    alpha = 2.5128624172523393539654752332184326538328336634026 * fl_nu / (radius * radius);
    kappa = (-1) * velocity * radius / (-1 + exp( -0.5 * alpha * radius * radius / fl_nu ));

    cout << "alpha: " << alpha << "\n";
    cout << "kappa: " << kappa << "\n";

    if ( parameters == "" )
    {
        stretching_r = 1;
        stretching_z = 1;
    }
    else
    {
        sscanf( parameters.c_str(), "%lf,%lf", &stretching_r, &stretching_z ); //"1,1"
    }
}


//////////////
// Destructor
BurgersVortex::~BurgersVortex() {}


/////////////////////////////
// Vortex Property Functions
Vector3d BurgersVortex::velocityCylinder( double r, double phi, double z )
{
    /* check for r == 0, otherwise there's a "divide by zero" and the function
     *  return NaN (inf * 0) */
    return Vector3d( -alpha * stretching_r * r,
                     (r == 0) ? 0 : (kappa / r * (1 - exp(-alpha * r * r / (2 * fl_nu)))),
                     2 * stretching_z * alpha * z );
}

Vector3d BurgersVortex::dudtCylinder( double r, double phi, double z )
{
    /* check for r == 0, otherwise there's a "divide by zero" and the function
     *  return NaN (inf * 0) */
    return Vector3d( stretching_r * stretching_r * alpha * alpha * r,
                     (r == 0) ? 0 : -((kappa * alpha * exp(-1 / 2 * alpha * r * r
                             / fl_nu)) / fl_nu - (kappa * (1 - exp(-1 / 2 * alpha * r
                             * r / fl_nu))) / (r * r)) * stretching_r * alpha * r,
                     4 * stretching_z * stretching_z * alpha * alpha * z);
}
