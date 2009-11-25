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

// TODOLIST
// TODO: move the equation of motion etc to a seperate class
// TODO: binary output
// TODO: write a better parameter class
// TODO: error checks on parameters
// TODO: Outputting data in HDF5 or whatever
// TODO: different ways of specifying what the vortex should look like

///////////
// Headers
#include "Main.h"

#include <stdio.h>

#include "..\external\getopt_pp.h"

#include "Vortex\BurgersVortex.h"

#include "Output\Output.h"
#include "Output\ByteOutput.h"
#include "Output\MatlabOutput.h"
#include "Output\TextOutput.h"
#include "Output\TecplotOutput.h"

#include "Emitter\Emitter.h"
#include "Emitter\GridEmitter.h"
#include "Emitter\GridOnceEmitter.h"
#include "Emitter\RandomEmitter.h"

#include "ParticleArray.h"
#include "Particle.h"


/////////////
// Namespace
using std::string;
using std::cout;


////////////
// Sutaato!
void moveParticles( Vortex *the_vortex, Emitter *the_emitter,
                    ParticleArray *particles, const Vortex3dParam &param )
{
    // Readability reference variables (Hopefully these will be optimized away).
    const Vector3d & gravity = param.gravity;
    const double & dt = param.dt;
    const double & beta = param.beta;
    const double & tau_a = param.tau_a;

    /* 
     * Include Drag, Gravity, Stresses and Added Mass in the equation of motion.
     * See Formula 11 in M.F. Cargnelutti and Portela's "Influence of the resuspension on 
     * the particle sedimentation in wall-bounded turbulent flows")
     */ 
    
#pragma omp parallel for
    for ( int p = 0; p < particles->getLength(); p++ )
    {
        // Get the particle.
        Particle particle = particles->getParticle( p );

        // Readability for rhs (hopefully optimized away).
        const Vector3d & p_pos = particle.getPos();
        const Vector3d & p_vel = particle.getVel();

        // Calculate the fluid velocity and accelertion at the position of the particle.
        const Vector3d & v_vel = the_vortex->getVelocityAt( p_pos ); 
        const Vector3d & Du_Dt = the_vortex->getDuDtAt( p_pos );

        // Particle equation of motion.
        const Vector3d dv = (1 / tau_a * (v_vel - p_vel) + (beta - 1) / (beta + 0.5)
                            * gravity + 1.5 / (beta + 0.5) * Du_Dt) * dt;

        // Give the particle his new velocity and position.
        particle.setVel( p_vel + dv ); 
        particle.setPos( p_pos + p_vel * dt ); 

        // Write the particle back to the array (maybe, and hopefully, it will be written 
        // directly to the array instead of making a temporary object like "particle" is.
        particles->setParticle( p, particle );
    }
}


// checkParticles():
// In:     - (1) The Vortex: To know how big the box is.
//         - (2) Emitter: to reset particles.
//         - (3) Particle Array, access to particle positions
//         - (4) Relative Time, needed for average fall time.
// In/Out: - (3) The particlearray of which the particle is to be removed.
// Out:    - (return) The relative time (fraction of going around time T_l) the particle spent in the box.
void checkParticles( Vortex *the_vortex, Emitter *the_emitter,
                     ParticleArray *particles, double relative_time,
                     double *average_fall_time, double *particles_out )
{
#pragma omp critical
    {
        // For each particle, check if it is outside or at the edge of the box; 
        // the edge counts as problematic too because then the interpolation fails.
        int p = 0;

        // Loop until there are no more particles left
        while ( p < particles->getLength() )
        {
            const Vector3d & p_pos = particles->getParticle( p ).getPos(); // Reference variable for readability

            // Check if the particle is at the edge, or outside of the box
            if ( the_vortex->outsideBox( p_pos ) )
            {
                double life_time = the_emitter->reset( p, relative_time, particles );
                
                (*average_fall_time) = ((*particles_out) * (*average_fall_time)
                        + life_time) / ((*particles_out) + 1);

                ++(*particles_out);
            }
            // The last particle which was moved to p might also be out-of-range and thus
            // also has to be checked. Only when nothing was moved can you safely increase p.
            else
            {
                p++;
            }
        }
    }
}

inline void writeProgress( int perc )
{
    printf( "[%d%%]\r", perc );
    fflush( stdout);
}

int main( int argc, char* argv[] )
{
    Vortex3dParam param;

    ///////////////////////
    // Parse the parameters
    parse( argc, argv, &param );
    printParam( param );

    ////////////////////
    // Making the Vortex
    Vortex *the_vortex;

    switch ( param.vortextype ) {
        case 1:
            the_vortex = new BurgersVortex( param );
            break;
        default:
            cout << "Unknown Vortex type";
            exit(1);
    }

    // See the comment at initInterpolate() as to why this is here.
    if ( param.interpolate == true )
        the_vortex->initInterpolate();

    /////////////////////
    // Making the Emitter
    Emitter *the_emitter;
    
    switch ( param.emittertype ) {
        case 1:
            the_emitter = new GridOnceEmitter( param );
            break;
        case 2:
            the_emitter = new GridEmitter( param );
            break;
        case 3:
            the_emitter = new RandomEmitter( param );
            break;
        default:
            cout << "Unknown Emitter type";
            exit( 1 );
    } 

    /////////////////////
    // Making the Output
    Output * outputter;

    switch ( param.outputformat ) {
        case 1:
            outputter = new ByteOutput( param, the_vortex );
            break;
        case 2:
            outputter = new MatlabOutput( param, the_vortex );
            break;
        case 3:
            outputter = new TextOutput( param, the_vortex );
            break;
        case 4:
            outputter = new TecplotOutput( param, the_vortex );
            break;
        default:
            cout << "Unknown outputformat.";
    }

    //////////////
    // Allocating memory for the array that holds the particles, and initializing (possibly emitting the first particles);
    ParticleArray particles( param.maxparticles );
    the_emitter->init( &particles );


    //////////////////
    // READY, SET, GO!

    // Needed to calculate the average fall time.
    double average_fall_time = 0;
    double particles_out = 0;

    //
    int t = 0;

    outputter->writeToFile( 0.0, particles );

    // If outputtype = 3, we're done so we can output and stop.
    if ( param.outputtype == 3 )
        exit( 0 );
 
    double time_next_output = param.outputinterval;

    // Maximum number iterations (when some particles never leave the box, or when reset_particles = 1)
    for ( t = 1; t <= param.max_t; t++ )
    {
        double relative_time = t * param.duration / param.max_t; // Relative time in fraction of T_l (going around time); goes from 0->duration
        double time = (t * param.dt) / param.max_t; // Absolute time in seconds.

        writeProgress( (t * 100) / param.max_t ); // Display progress on stdout, e.g. [45%]

        // Move the particles
        moveParticles( the_vortex, the_emitter, &particles, param );
        checkParticles( the_vortex, the_emitter, &particles, relative_time,
                &average_fall_time, &particles_out ); // Check the particles to see if any of them are outside the box. Also updates

        // Write to file
        if ( relative_time >= time_next_output )
        {
            outputter->writeToFile( time, particles );
            time_next_output += param.outputinterval;
        }

        the_emitter->update( relative_time, &particles );

        // Break when there are no more particles to plot.
        if ( particles.getLength() == 0 )
            break;
    }


    delete the_vortex;
    delete the_emitter;
    delete outputter;
    return 0;
}

////////
// Help 
void show_help()
{
    printf( 
            "\nGeneral Options:\n"
            "  -h, --help                                         Produce help message.\n"
            "      --reset_particles <int> (=0)                   0 = Remove particle when it leaves the box.\n"
            "                                                     1 = Reset particle to start position.\n"
            "      --datafile <string> (=test.txt)                The path to the output file.\n"
            "      --outputtype <int> (=1)                        1: Particle index, position and absolute velocity.\n"
            "                                                     2: Relative concentration.\n"
            "                                                     3: Vortex VectorField.\n"
            "      --outputformat <int> (=1)                      1: Byte\n"
            "                                                     2: Matlab\n"
            "                                                     3: Text\n"
            "                                                     4: Tecplot\n"
            "      --outputinterval <double> (=1.0)               Write every <n cycles | relative_time>.\n"
            "      --outputintervalmode <int> (=1)                Change behavior of outputinterval.\n"
            "                                                       1: every n cycles\n"
            "                                                       2: every relative_time\n"
            "      --interpolate                                  Use interpolation instead of direct evaluation.\n"
            "      --duration <double> (=1.0)                     Duration of computation as fraction T_l.\n"
            "      --maxparticles <int> (=1000)                   Maximum particles, no new particles will be emitted if the number of particles exceeds this parameter.\n"
            "      --gravity <double> (=1.0)                      Fraction of 9.81m/s^2.\n"
            "      --dtscale <double> (=0.5)                      dt = dtscale * systemtime.\n"
            "      --rotategrav                                   Rotate not the vortex but the gravity.\n"
            "\n"
            "Output Options\n"
            "\n"
            "      Vortex Properties:\n"
            "      --vortextype <int> (=1)                        1: Burgers Vortex, (stretching_r,stretching_z)\n"
            "      --radius <double> (=0.1)                       The typical radius of the vortex.\n"
            "      --velocity <double> (=0.001)                   Rotational velocity at the radius.\n"
            "      --parameters <string> (=)                      Other values the specified vortex might need, seperated by comma (see vortextype).\n"
            "      --angle <double> (=0.0)                        Angle of the vortex with the z-axis.\n"
            "      --roi <string> (=[-5:101:5,-5:101:5,-5:101:5]) Region of interest, middle value is not used when not using interpolation.\n"
            "\n"
            "      Particle/Emitter Properties:\n"
            "      --emittertype <int> (=1)                       1: Emit the grid specified with --dimensions once.\n"
            "                                                     2: Emit continuously from the --dimensions grid.\n"
            "                                                     3: Emit continuously and random from the --dimensions grid.\n"
            "      --dimensions <string> (=[-4:30:4,0:1:0,4:1:4]) X x Y x Z grid of particles\n"
            "                                                     Format: [start_x:steps_x:end_x,\n"
            "                                                              start_y:steps_y:end_y,\n"
            "                                                              start_z:steps_z:end_z].\n"
            "      --p_rate <double> (=100.0)                     Amount of particles emitted per T_l.\n"
            "      --p_density <double> (=1000.0)                 Density of the particles (kg/m3).\n"
            "      --p_diameter <double> (=5.0E-5)                Diameter of the particles (m).\n"
            "      --fl_mu <double> (=1.8e-005)                   Fluid Dynamic Viscocity (Pa s).\n"
            "      --fl_density <double> (=1.0)                   Fluid Density (kg/m3).\n"
            "      --p_velocity <double> (=0.0)                   Particle's initial velocity (z-direction only) in terms of the terminal velocity.\n"
          );
}


/////////////////////
// Parse commandline 
void parse( int argc, char* argv[], Vortex3dParam *param ) {
    
    using GetOpt::GetOpt_pp;
    using GetOpt::Option;
    using GetOpt::OptionPresent;

    GetOpt_pp ops(argc, argv);

    if ( ops >> OptionPresent( 'h', "help" ) )
    {
        show_help();
        exit( 1 );
    }

    // Parse the cmdline
    ops >> Option( 'a', "vortextype", param->vortextype, 1 )
        >> Option( 'a', "radius", param->radius, 0.1 )
        >> Option( 'a', "velocity", param->velocity, 0.001 )
        >> Option( 'a', "parameters", param->parameters, "" )
        >> Option( 'a', "angle", param->angle, 0.0 )
        >> Option( 'a', "roi", param->roi, "[-5:101:5,-5:101:5,-5:101:5]" )
        >> Option( 'a', "emittertype", param->emittertype, 1 )
        >> Option( 'a', "dimensions", param->dimensions, "[-4:30:4,0:1:0,4:1:4]" )
        >> Option( 'a', "p_rate", param->p_rate, 100.0 )
        >> Option( 'a', "p_density", param->p_density, 1000.0 )
        >> Option( 'a', "p_diameter", param->p_diameter, 5E-5 )
        >> Option( 'a', "fl_mu", param->fl_mu, 18E-6 )
        >> Option( 'a', "fl_density", param->fl_density, 1.0 )
        >> Option( 'a', "p_velocity", param->p_velocity, 0.0 )
        >> Option( 'a', "reset_particles", param->reset_particles, 0 )
        >> Option( 'a', "datafile", param->datafile, "test.data" )
        >> Option( 'a', "outputtype", param->outputtype, 1 )
        >> Option( 'a', "outputformat", param->outputformat, 1 )
        >> Option( 'a', "outputinterval", param->outputinterval, 1.0 )
        >> Option( 'a', "outputintervalmethod", param->outputintervalmethod, 1 )
        >> OptionPresent( 'a', "interpolate", param->interpolate )
        >> Option( 'a', "duration", param->duration, 1.0 )
        >> Option( 'a', "maxparticles", param->maxparticles, 1000 )
        >> Option( 'a', "gravity", param->grav, 1 )
        >> OptionPresent( 'a', "rotategrav", param->rotategrav )
        >> Option( 'a', "dtscale", param->dtscale, 0.5 )
     ;
    //////////////////////////////////////
    // Do some checking on the parameters
    
    // Gravity (depends on whether you rotate the gravity or the vortex)
    if (!param->rotategrav)
        param->gravity = Vector3d( 0, 0, -9.81 * param->grav );
    else
        param->gravity = Vector3d( 0, sin(param->angle), -cos(param->angle) );
    
    // Read grid+delimiter+deltas for ROI.
    readGridDelimiterDelta( param->roi, param->radius, &param->roi_grid, &param->roi_delimiter, 
                            &param->roi_dx, &param->roi_dy, &param->roi_dz );
    // Read grid+delimiter+deltas for Emitter.
    readGridDelimiterDelta( param->dimensions, param->radius, &param->emitter_grid, &param->emitter_delimiter, 
                            &param->emitter_dx, &param->emitter_dy, &param->emitter_dz );

    // Calculate the amount of particles in one grid.
    param->p_N = product( param->emitter_grid );

    // If the outputtype is 3, you want the vortex velocity field. Therefore, interpolate should be 1
    if ( param->outputtype == 3 )
        param->interpolate = true;

    // In case of the GridEmitters (1 and 2) maxparticles should be at least the size of one grid.
    if ( param->emittertype == 1 || param->emittertype == 2 )
    {
        if ( param->maxparticles < param->p_N )
            param->maxparticles = param->p_N;
    }

    /////////////////////////////////////////////////////////
    // Calculate some extra parameters from the cmdline ones
    param->fl_nu = param->fl_mu / param->fl_density;
    param->beta = param->p_density / param->fl_density; 

    param->systemtime = param->p_density * param->p_diameter * param->p_diameter / (18 * param->fl_nu);
    param->tau_a = (param->beta + 0.5) / param->beta * param->systemtime;

    param->dt = param->dtscale * param->systemtime;
    
    param->v_ref = 2 * PI * param->radius / param->velocity;

    param->terminal_velocity = (1 - (1 / param->beta)) * param->systemtime * -9.81;
    param->p_velocity = param->p_velocity * param->terminal_velocity;

    param->max_t = (int) ( param->duration * param->v_ref / param->dt );

    /////////////////////////////////////////////////////////////
    // And some more checking?

    // Can't output more data than we have:
    if ( param->outputintervalmethod == 1 )
    {
        if ( param->outputinterval < 1 )
            printf( "Error: outputinterval should be >= 1 on intervalmode 1." );

        param->outputinterval = param->outputinterval * param->dt / param->v_ref;
    }
    else if ( param->outputintervalmethod == 2 )
    {
        if ( param->outputinterval <= 0 )
            param->outputinterval = 0;
        
        // If the absolute time between outputs is smaller than dt (the iteration) 
        // exit with error.
        if ( param->outputinterval * param->v_ref / param->dt < 1 ) 
        {
            printf( "Error: outputcycle is smaller than 1: %g.\n"     
                    "       Please lower dtscale to accomodate.\n", param->outputinterval * param->v_ref / param->dt );
            exit(1);
        }
    }
    else
    {
        printf( "Unknown outputintervalmethod: %d\n", param->outputintervalmethod );
        exit(1);
    }
}

///////////////////////////
// Parse formatted strings
void readGridDelimiterDelta( const string &fstring, const double &radius, TGrid *grid, 
                             TDelimiter *delimiter, double *dx, double *dy, double *dz )
{
    double x1, x2, y1, y2, z1, z2;
    int X, Y, Z;

    // Read the values into the variables.
    sscanf( fstring.c_str(), "[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]", &x1, &X, &x2,
            &y1, &Y, &y2, &z1, &Z, &z2 ); //e.g. [-4:30:4,0:1:0,4:1:4]"

    *delimiter = x1 * radius, x2 * radius,
                 y1 * radius, y2 * radius,
                 z1 * radius, z2 * radius;
    *grid = X, Y, Z;

    // Readability:
    const TDelimiter & _delimiter = *delimiter;
    const TGrid & _grid = *grid;

    // Set dx/dy/dz, and do some checks on them.
    *dx = ( _delimiter(0, 1) - _delimiter(0, 0) ) / ( _grid(0) - 1 );
    *dy = ( _delimiter(1, 1) - _delimiter(1, 0) ) / ( _grid(1) - 1 );
    *dz = ( _delimiter(2, 1) - _delimiter(2, 0) ) / ( _grid(2) - 1 );

    if ( _grid(0) <= 1 )
        *dx = 0;

    if ( _grid(1) <= 1 )
        *dy = 0;

    if ( _grid(2) <= 1 )
        *dz = 0;
}

void printParam( const Vortex3dParam &param )
{
    printf( "System Time (tau_p): %.5g\n", param.systemtime );
    printf( "System Time (tau_a): %.5g\n", param.tau_a );
    printf( "Maximum time:        %.5g\n", param.max_t * param.dt );
    printf( "Timestep size (dt):  %.5g\n", param.dt );
    printf( "Amount of timesteps: %d\n",   param.max_t );
}