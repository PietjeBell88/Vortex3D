#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Main.h"
#include "Emitter.h"
#include "Particle.cpp"
#include <blitz/array.h>

using namespace std;
using namespace blitz;

class Emitter {
protected:
    TinyVector<int,3> p_grid;                        // X x Y x Z grid of particles
    TinyMatrix<double,3,2> delimiter;                // offsets
    double dx,dy,dz;                                // delta's

    double p_density;                                // Density of the particles (kg/m3)
    double p_diameter;                                // Diameter of the particles (m).
    double p_velocity;
    double p_rate;
    int reset_particles;
    int p_N;                                        // Total amount of particles.

public:
    Emitter(const double &p_density, const double &p_diameter, const double &p_velocity, 
        const string &dimensions, const double &radius, const double &p_rate, const int &reset_particles) 
    {    
        this->p_density = p_density;
        this->p_diameter = p_diameter;
        this->p_velocity = p_velocity;
        this->reset_particles = reset_particles;
        this->p_rate = p_rate;

        double x1,x2,y1,y2,z1,z2;
        int X,Y,Z;
        sscanf(dimensions.c_str(),"[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",
            &x1,&X,&x2,&y1,&Y,&y2,&z1,&Z,&z2); //e.g. [-4:30:4,0:1:0,4:1:4]"
        delimiter = x1*radius, x2*radius,
            y1*radius, y2*radius,
            z1*radius, z2*radius;

        p_grid = X,Y,Z;
        p_N = product(p_grid);

        dx = (delimiter(0,1)-delimiter(0,0))/(p_grid(0)-1);
        dy = (delimiter(1,1)-delimiter(1,0))/(p_grid(1)-1);
        dz = (delimiter(2,1)-delimiter(2,0))/(p_grid(2)-1);

        if (p_grid(0) <= 1) {
            dx = 0;
        }
        if (p_grid(1) <= 1) {
            dy = 0;
        }
        if (p_grid(2) <= 1) {
            dz = 0;
        }
    }

    virtual void init(ParticleArray *particles) = 0;

    virtual void update(const double &relative_time, ParticleArray *particles) = 0;

    /*
     * reset():
     * In:     - (1) The number of the particle to be reset.
     *         - (2) The current relative time (as fraction of going around time T_l).
     * In/Out: - (3) The particlearray of which the particle is to be removed.
     * Out:    - (return) The relative time (fraction of going around time T_l) the particle spent in the box. 
     */
    virtual double reset(const int &p, const double &relative_time, ParticleArray *particles) 
    {
        Particle removed = particles->remove(p);  // The particle that be removed;
        
        if (reset_particles != 0) {
            particles->add(getStartPos(p), getStartVel(p), relative_time);
        }

        return (relative_time - removed.getSpawnTime());
    }

    virtual vector3d getStartPos(const int &p) = 0;

    virtual vector3d getStartVel(const int &p) = 0;
};