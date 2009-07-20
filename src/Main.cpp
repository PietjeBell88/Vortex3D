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
// error checks on parameters
// Outputting data in HDF5 or whatever
// different ways of specifying what the vortex should look like 

#pragma once
#define _CRT_SECURE_NO_WARNINGS


///////////
// Headers
#include <stdio.h>
#include <fstream>
#include <boost/program_options.hpp>

#include "Main.h"

#include "BurgersVortex.h"

#include "GridEmitter.h"
#include "GridOnceEmitter.h"
#include "RandomEmitter.h"

#include "Particle.h"


/////////////
// Namespace
using namespace std;
using namespace blitz;
namespace po = boost::program_options;



/////////////
void ReadROI(const string &roi, const double &radius, Settings *options) {
    double x1,x2,y1,y2,z1,z2;
    int X,Y,Z;

    sscanf(roi.c_str(),"[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",&x1,&X,&x2,&y1,&Y,&y2,&z1,&Z,&z2); //"[-4:30:4,0:1:0,4:1:4]", also see Emitter
    options->delimiter = x1*radius, x2*radius,
        y1*radius, y2*radius,
        z1*radius, z2*radius;
    options->grid = X,Y,Z;

    TinyVector<int,3> & grid = options->grid;    // readability
    TinyMatrix<double,3,2> & delimiter = options->delimiter;        //readability

    options->dx = (delimiter(0,1)-delimiter(0,0))/(grid(0)-1);
    options->dy = (delimiter(1,1)-delimiter(1,0))/(grid(1)-1);
    options->dz = (delimiter(2,1)-delimiter(2,0))/(grid(2)-1);

    if (grid(0) <= 1) {
        options->dx = 0;
    }
    if (grid(1) <= 1) {
        options->dy = 0;
    }
    if (grid(2) <= 1) {
        options->dz = 0;
    }
}

// I donno where this function belongs.. a lot of code comes straight from Vortex.cpp, it even uses its ROI
void GetConcentration(ParticleArray *particles, const Settings &options, ScalarField *concentration) {        // *TODO* const ParticleArray
    *concentration = 0;
    const TinyMatrix<double,3,2> &delimiter = options.delimiter;            //readability
    const TinyVector<int,3> & grid = options.grid;                        //readability
    const double & dx = options.dx;                                        //readability
    const double & dy = options.dy;                                        //readability
    const double & dz = options.dz;                                        //readability
    const int & length = particles->GetLength();

    for (int p = 0; p < length; p++) {
        const Vector3d & pos = particles->GetParticle(p).GetPos();

        int i = static_cast<int> ( floor((pos(0) - delimiter(0,0))/dx) );
        int j = static_cast<int> ( floor((pos(1) - delimiter(1,0))/dy) );
        int k = static_cast<int> ( floor((pos(2) - delimiter(2,0))/dz) );

        ++(*concentration)(i,j,k);
    }
    if (sum(*concentration) != 0) { 
        *concentration = (*concentration)/sum(*concentration); 
    }
}

void MoveParticles(Vortex *theVortex, Emitter *theEmitter, ParticleArray *particles, const Settings &options) {

    // What variables are used in this function, and make some nice reference variables for them to improve readability.
    const double & systemtime = options.systemtime;    
    const double & gravity = options.gravity;    
    const double & dt = options.dt;
    const int & reset_particles = options.reset_particles;
    const double & beta = options.beta;                                            
    const double & tau_a = options.tau_a;

    /*
    * This function calculates the new particle velocity and position based
    * on each particles current position and velocity and the vortex's velocity vector field.
    * It return the amount of remaining particles, so that the main loop knows when to stop,
    * that is, when amount of particles in the box equals 0.
    */

    /*
    * Move the particles, this can be multithreaded, if it will, Vector3d u and dv have to go inside the loop.
    * Only the case where the particles are being reset can be included in this (possibly multithreaded) loop.
    * Include Drag, Gravity, Stresses and Added Mass in the equation of motion.
    * See Formula 11 in M.F. Cargnelutti and Portela's "Influence of the resuspension on the particle sedimentation in wall-bounded turbulent flows")
    */
#pragma omp parallel for
    for (int p = 0; p < particles->GetLength(); p++) {

        Vector3d & p_pos = particles->GetParticle(p).GetPos();                                                // Reference variable for readability
        Vector3d & p_vel = particles->GetParticle(p).GetVel();                                                // Reference variable for readability

        Vector3d u, dv, v_vel, Du_Dt;    

        v_vel = theVortex->GetVelocity( p_pos );                                                // Calculate the fluid velocity at the position of the particle by either interpolation or evaluation
        Du_Dt = theVortex->GetDuDt( p_pos );                

        dv = (1/tau_a * (v_vel - p_vel) + (beta - 1)/(beta + 0.5) * Vector3d(0,0,gravity) + 1.5/(beta + 0.5) * Du_Dt) * dt;    // equation of motion * dt    
        //cout << (beta - 1)/(beta + 0.5) * Vector3d(0,0,gravity) << "\n";
        //cout << dv << "\n\n";

        p_vel = p_vel + dv;                                                            // Calculate the new velocity of this particle
        p_pos = p_pos + p_vel * dt;                                            // Calculate the new position of this particle by adding the new velocity multiplied with dt
    }
}


// CheckParticles():
// In:     - (1) The Vortex: To know how big the box is.
//         - (2) Emitter: to reset particles.
//         - (3) Particle Array, access to particle positions
//         - (4) Relative Time, needed for average fall time.
// In/Out: - (3) The particlearray of which the particle is to be removed.
// Out:    - (return) The relative time (fraction of going around time T_l) the particle spent in the box. 
void CheckParticles(Vortex *theVortex, Emitter *theEmitter, ParticleArray *particles, const double &relative_time, double *average_fall_time, double *particles_out) 
{
#pragma omp critical
    {

        // For each particle, check if it is outside or at the edge of the box; the
        // edge counts as problematic too because then the interpolation fails.
        int p = 0;
        while (p < particles->GetLength()) {                                                                    // loop until there are no more particles left

            const Vector3d & p_pos = particles->GetParticle(p).GetPos();                                                // Reference variable for readability
            const Vector3d & p_vel = particles->GetParticle(p).GetVel();                                                // Reference variable for readability

            if (theVortex->OutsideBox( p_pos )) {                                               // Check if the particle is at the edge, or outside of the box
                double life_time = theEmitter->Reset(p, relative_time, particles);
                (*average_fall_time) = ( (*particles_out) * (*average_fall_time) + life_time ) / ( (*particles_out) + 1);
                ++(*particles_out);
            }
            else {                                                                                // The last particle which was moved to p might also be out-of-range and thus
                p++;                                                                            // also has to be checked. Only when nothing was moved can you safely increase p.
            }
        }
    }
}

inline void WriteProgress(int perc) {
    printf("[%d%%]\r", perc );            
    fflush(stdout);
}

inline void WriteToFile(const double time, const Settings &options, FILE * f, ParticleArray *particles, Vortex *thevortex) {    // *TODO* const ParticleArray
    static bool first_call = true;

    // Print the particle positions
    if (options.outputtype == 1 && !options.tecplot) {
        if (first_call) {
            fprintf(f,"#T     N     X     Y     Z     Speed\n");
            first_call = false;
        }
        for (int i = 0; i < particles->GetLength(); i++) {
            // Readability
            const int & num = particles->GetParticle(i).GetNum();
            const Vector3d & pos = particles->GetParticle(i).GetPos();
            const double & speed = particles->GetParticle(i).GetSpeed();

            // Write output to file
            fprintf(f,"%e     %d     %e     %e     %e     %e\n", time, num, pos(0), pos(1), pos(2), speed);
        }
        fprintf(f,"\n");
    }

    if (options.outputtype == 1 && options.tecplot) {
        if (first_call) {
            fprintf(f,"#T     N     X     Y     Z     Speed\n");
            first_call = false;
        }
        for (int i = 0; i < particles->GetLength(); i++) {
            // Readability
            const int & num = particles->GetParticle(i).GetNum();
            const Vector3d & pos = particles->GetParticle(i).GetPos();
            const double & speed = particles->GetParticle(i).GetSpeed();

            // Write output to file
            fprintf(f,"%e     %d     %e     %e     %e     %e\n", time, num, pos(0), pos(1), pos(2), speed);
        }
        fprintf(f,"\n");
    }


    // Print the concentration
    else if (options.outputtype == 2 && !options.tecplot) {
        if (first_call) {
            fprintf(f,"#T     X     Y     Z     C\n");
            first_call = false;
        }
        ScalarField concentration(options.grid(0), options.grid(1), options.grid(2));
        GetConcentration(particles, options, &concentration);
        for (int i = 0; i < options.grid(0); i++) {
            double x = options.delimiter(0,0) + i * options.dx;
            for (int j = 0; j < options.grid(1); j++) {
                double y = options.delimiter(1,0) + j * options.dy;
                for (int k = 0; k < options.grid(2); k++) {
                    double z = options.delimiter(2,0) + k * options.dz;
                    fprintf(f,"%e     %e     %e     %e     %e\n", time, x, y, z, concentration(i,j,k));
                }
            }
        }
        fprintf(f,"\n");
    }

    else if (options.outputtype == 2 && options.tecplot) {
        if (first_call) {
            fprintf(f,"TITLE=\"Example: Multi-Zone XY Line Plot Wwith Variable Sharing\"\n \
                      VARIABLES=\"X\" \"Y\" \"Z\" \"C\"\n \
                      ZONE T = \"%f seconds\", I = %d, J = %d, K = %d, DATAPACKING = POINT\n", 
                      time, options.grid(0), options.grid(1), options.grid(2));
            first_call = false;

            ScalarField concentration(options.grid(0), options.grid(1), options.grid(2));
            GetConcentration(particles, options, &concentration);
            
            for (int i = 0; i < options.grid(0); i++) {
                double x = options.delimiter(0,0) + i * options.dx;
                for (int j = 0; j < options.grid(1); j++) {
                    double y = options.delimiter(1,0) + j * options.dy;
                    for (int k = 0; k < options.grid(2); k++) {
                        double z = options.delimiter(2,0) + k * options.dz;
                        fprintf(f,"%e     %e     %e     %e     %e\n", time, x, y, z, concentration(i,j,k));
                    }
                }
            }
            fprintf(f,"\n");
        }

        else {
            fprintf(f,"ZONE T = \"%f seconds\", I = %d, J = %d, K = %d, DATAPACKING = POINT\n \
                      VARSHARELIST=([1-3]=1)\n", 
                      time, options.grid(0), options.grid(1), options.grid(2));

            ScalarField concentration(options.grid(0), options.grid(1), options.grid(2));
            GetConcentration(particles, options, &concentration);

            for (int i = 0; i < options.grid(0); i++) {
                for (int j = 0; j < options.grid(1); j++) {
                    for (int k = 0; k < options.grid(2); k++) {
                        fprintf(f,"%e\n", concentration(i,j,k));
                    }
                }
            }
            fprintf(f,"\n");
        }
    }


    // Print the Vortex Field
    // Mayavi2 version
    else if (options.outputtype == 3 && !options.tecplot) {
        if (first_call) {
            fprintf(f,"#T     X     Y     Z     UX     UY     UZ\n");
            first_call = false;
        }
        const VectorField &v = thevortex->GetVectorField();

        for (int i = 0; i < options.grid(0); i++) {
            double x = options.delimiter(0,0) + i * options.dx;
            for (int j = 0; j < options.grid(1); j++) {
                double y = options.delimiter(1,0) + j * options.dy;
                for (int k = 0; k < options.grid(2); k++) {
                    double z = options.delimiter(2,0) + k * options.dz;
                    const Vector3d &vel = v(i,j,k);
                    fprintf(f,"%e     %e     %e     %e     %e     %e     %e\n", time, x, y, z, vel(0), vel(1), vel(2));
                }
            }
        }
        fprintf(f,"\n");
    }

    // TECPLOT Version
    else if (options.outputtype == 3 && options.tecplot) {
        if (first_call) {
            fprintf(f,"TITLE=\"Simple Data File\"\nVARIABLES=\"X\" \"Y\" \"Z\" \"UX\" \"UY\" \"UZ\"\nZONE I = %d, J = %d, K = %d, DATAPACKING = POINT\n", options.grid(0), options.grid(1), options.grid(2));
            first_call = false;
        }
        const VectorField &v = thevortex->GetVectorField();

        for (int i = 0; i < options.grid(0); i++) {
            double x = options.delimiter(0,0) + i * options.dx;
            for (int j = 0; j < options.grid(1); j++) {
                double y = options.delimiter(1,0) + j * options.dy;
                for (int k = 0; k < options.grid(2); k++) {
                    double z = options.delimiter(2,0) + k * options.dz;
                    const Vector3d &vel = v(i,j,k);
                    fprintf(f,"%e     %e     %e     %e     %e     %e\n", x, y, z, vel(0), vel(1), vel(2));
                }
            }
        }
        fprintf(f,"\n");
    }


    else {
        cout << "Invalid output method.";
        exit(1);
    }
}


int main(int ac, char* av[])
{
    int vortextype, reset_particles, maxparticles, grav, emittertype, outputtype;
    double radius, velocity, angle, p_density, p_diameter, fl_mu, fl_density, duration, p_velocity, outputinterval, dtscale, p_rate;
    bool interpolate;
    string parameters, roi, dimensions, datafile;


    // Declare the supported options.
    po::options_description vortex("Vortex Properties");
    vortex.add_options()
        ("vortextype", po::value<int>(&vortextype)->default_value(1), "1: Burgers Vortex, (stretching_r,stretching_z)")
        ("radius", po::value<double>(&radius)->default_value(0.1), "The typical radius of the vortex.")
        ("velocity", po::value<double>(&velocity)->default_value(0.001), "Rotational velocity at the radius.")    
        ("parameters", po::value<string>(&parameters)->default_value(""), "Other values the specified vortex might need, seperated by comma (see vortextype).")
        ("angle", po::value<double>(&angle)->default_value(0.0), "Angle of the vortex with the z-axis, clockwise when looking along the x-axis in a rh coordinate system.")
        ("roi", po::value<string>(&roi)->default_value("[-5:101:5,-5:101:5,-5:101:5]"), "Region of interest, middle value is not used when not using interpolation.")
        ;

    po::options_description emitter("Particle/Emitter Properties");
    emitter.add_options()
        ("emittertype", po::value<int>(&emittertype)->default_value(1), "1: Emit the grid specified with --dimensions once.\n2: Emit continuously from the --dimensions grid.\n3: Emit continuously and random from the --dimensions grid.")
        ("dimensions", po::value<string>(&dimensions)->default_value("[-4:30:4,0:1:0,4:1:4]"), "X x Y x Z grid of particles\nFormat: [start_x:steps_x:end_x, start_y:steps_y:end_y, start_z:steps_z:end_z]. Where start_x etc are given as a fraction of the radius of the vortex.")    
        ("p_rate", po::value<double>(&p_rate)->default_value(100), "Amount of particles emitted per T_l")
        ("p_density", po::value<double>(&p_density)->default_value(1000.0), "Density of the particles (kg/m3).")
        ("p_diameter", po::value<double>(&p_diameter)->default_value(5E-5), "Diameter of the particles (m).")
        ("fl_mu", po::value<double>(&fl_mu)->default_value(18E-6), "Fluid Dynamic Viscocity (Pa s).")
        ("fl_density", po::value<double>(&fl_density)->default_value(1.0), "Fluid Density (kg/m3).")
        ("p_velocity", po::value<double>(&p_velocity)->default_value(0), "Particle's initial velocity (z-direction only) in terms of the terminal velocity.")
        ;

    po::options_description general("General Options");
    general.add_options()
        ("help", "Produce help message.")
        ("reset_particles", po::value<int>(&reset_particles)->default_value(0), "0 = Remove particle when it leaves the box.\n1 = Reset particle to start position.")
        ("datafile", po::value<string>(&datafile)->default_value("test.txt"), "The path to the output file.")
        ("outputtype", po::value<int>(&outputtype)->default_value(1), "1: Particle index, position and absolute velocity.\n2: Relative concentration.\n3: Vortex VectorField")
        ("outputinterval", po::value<double>(&outputinterval)->default_value(0), "Fraction of T_l, e.g. 0.1 = emit every 0.1th T_l (0 = output every cycle")
        ("interpolate", po::value<bool>(&interpolate)->default_value(false), "true = Interpolate velocity field.\nfalse = Evaluate for every particle.")
        ("duration", po::value<double>(&duration)->default_value(1.0), "Duration of computation as fraction T_l")
        ("maxparticles", po::value<int>(&maxparticles)->default_value(1000), "Maximum particles, no new particles will be emitted if the number of particles exceeds this parameter.")
        ("gravity", po::value<int>(&grav)->default_value(1), "0 = No gravity (0)\n1 = Gravity (9.81).")
        ("dtscale", po::value<double>(&dtscale)->default_value(0.5), "dt = dtscale * systemtime.")
        ("tecplot", "Print data format for tecplot")
        ;

    po::options_description cmdline_options;
    cmdline_options.add(vortex).add(emitter).add(general);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).run(), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << general << "\n";
        cout << vortex << "\n";
        cout << emitter << "\n";
        return 1;
    }

    //gedoe
    Settings options;
    options.tecplot = vm.count("tecplot");
    options.systemtime = p_density * p_diameter * p_diameter /(18*fl_mu);    
    options.gravity = -9.81 * grav;
    options.dt = dtscale*options.systemtime;
    options.reset_particles = reset_particles;
    options.beta = p_density/fl_density;                                            // This ratio is used in the equation of motion
    options.tau_a = (options.beta + 0.5)/options.beta *options.systemtime;
    options.datafile = datafile;
    options.outputtype = outputtype;

    p_velocity = (1-(1/options.beta)) * p_velocity * options.systemtime * -9.81;
    ReadROI(roi, radius, &options);                                                    // Write some settings to options

    // Parameter Checking:
    // If the outputtype is 3, you want the vortex velocity field. Therefore, interpolate should be 1 
    if ( outputtype == 3) {
        interpolate = true;
    }
    // In case of the GridEmitters (1 and 2) maxparticles should be at least the size of one grid.
    if ( emittertype == 1 || emittertype == 2) {
        double x1,x2,y1,y2,z1,z2;
        int X,Y,Z;
        sscanf(dimensions.c_str(),"[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",           
            &x1,&X,&x2,&y1,&Y,&y2,&z1,&Z,&z2); //e.g. [-4:30:4,0:1:0,4:1:4]"
        int p_N = X*Y*Z;
        if (maxparticles < p_N) { 
            maxparticles = p_N;
        }
    }

    FILE * stat = fopen( datafile.c_str(), "w" );
    fprintf(stat, "#Settings: vortextype %d, radius %e, velocity %e, parameters \"%s\", angle %e, roi \"%s\", emittertype %d", vortextype, radius, velocity, parameters.c_str(), angle, roi.c_str(), emittertype);
    fprintf(stat, ", dimensions \"%s\", p_rate %e, p_density %e, p_diameter %e, fl_mu %e, fl_density %e, p_velocity %e", dimensions.c_str(), p_rate, p_density, p_diameter, fl_mu, fl_density, p_velocity);
    fprintf(stat, ", reset_particles %d, datafile \"%s\", outputtype %d, outputinterval %e, interpolate %d, duration %e, maxparticles %d, gravity %d\n", reset_particles, datafile.c_str(), outputtype, outputinterval, interpolate, duration, maxparticles, grav);
    fprintf(stat, "#Other: Systemtime %e, TerminalVelocity %e\n", options.systemtime, (1-(1/options.beta)) * options.systemtime * -9.81);
    fclose(stat);



    // Making the Vortex
    Vortex * theVortex;

    if ( vortextype == 1 ) {
        theVortex = new BurgersVortex(parameters, radius, velocity, angle, fl_mu, fl_density, interpolate, roi); 
    }
    else { 
        cout << "Unknown Vortex type";
        exit(1);        
    }

    // See the comment at InitInterpolate() as to why this is here.
    if (interpolate == true) {
        theVortex->InitInterpolate();
    }

    // Making the Emitter
    //GridOnceEmitter tempEmitter(p_density, p_diameter, p_velocity, dimensions, radius, p_rate, reset_particles);
    Emitter * theEmitter;

    if ( emittertype == 1 ) {
        theEmitter = new GridOnceEmitter(p_density, p_diameter, p_velocity, dimensions, radius, p_rate, reset_particles);
    }
    else if( emittertype == 2 ) {
        theEmitter = new GridEmitter(p_density, p_diameter, p_velocity, dimensions, radius, p_rate, reset_particles);
    }
    else if ( emittertype == 3 ) {
        theEmitter = new RandomEmitter(p_density, p_diameter, p_velocity, dimensions, radius, p_rate, reset_particles);
    }
    else {
        cout << "Unknown Emitter type";
        exit(1);        
    }


    // Allocating memory for the array that holds the particles, and initializing (possibly emitting the first particles);
    ParticleArray particles(maxparticles);
    theEmitter->Init(&particles);



    // READY, SET, GO!
    int max_t = duration * 2* 3.14159265359289783238 * radius / velocity / options.dt;
    int interval = outputinterval; 

    // Needed to calculate the average fall time. 
    double average_fall_time = 0;
    double particles_out = 0;

    // 
    int t = 0;

    FILE * f = fopen( datafile.c_str(), "a" );                                            // C style fprintf's instead of fstream and stuff, because i read somewhere that fprintf is faster
    WriteToFile(0, options, f, &particles, theVortex);

    // If outputtype = 3, we're done so we can output and stop.
    if (outputtype == 3) {
        fclose(f);
        exit(0);
    }


    for (t = 1; t <= max_t; t++) {                                                        // Maximum number iterations (when some particles never leave the box, or when reset_particles = 1)
        double relative_time = t * duration / max_t;                                    // Relative time in fraction of T_l (going around time); goes from 0->duration 
        double time = (t*options.dt)/max_t;                                                // Absolute time in seconds.

        WriteProgress(  (t*100)/max_t );                                                // Display progress on stdout, e.g. [45%]

        // move the particles
        MoveParticles(theVortex, theEmitter, &particles, options);        // Move the particles, and keep track of how many particles remain inside the box.
        CheckParticles(theVortex, theEmitter, &particles, relative_time, &average_fall_time, &particles_out);        // Check the particles to see if any of them are outside the box. Also updates 

        // write to file
        if (relative_time > interval) {                                                    // 
            WriteToFile(time, options, f, &particles, theVortex);
            interval += outputinterval;
        }

        theEmitter->Update(relative_time, &particles);

        if( particles.GetLength() == 0 ) { break; }                                        // Break when there are no more particles to plot

    }

    fprintf(f, "\n# Amount of particles that left the box: %e\n# Mean amount of time the particles spent in the box:%e",particles_out,average_fall_time);
    fclose(f);    
    delete theVortex;
    delete theEmitter;
    return 0;
}

