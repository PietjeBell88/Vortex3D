#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Main.h"
#include "Emitter.cpp"
#include "Emitter.h"
#include "Particle.cpp"
#include <blitz/array.h>
#include "..\external\MTRand.h" 

using namespace std;
using namespace blitz;

class RandomEmitter : public Emitter {
protected:
    double last_emit_time;                            // Relative time at which the last particles were emitted
public:
    RandomEmitter(const double &p_density, const double &p_diameter, const double &p_velocity, const string &dimensions, const double &radius, const double &p_rate, const int &reset_particles) : Emitter(p_density, p_diameter, p_velocity, dimensions, radius, p_rate, reset_particles) {    // "[double:int:double,double:int:double,double:int:double] <-- relative to R_v
        last_emit_time = 0;
    }

    virtual void init(ParticleArray *particles) {
    }

    virtual void update(const double &relative_time, ParticleArray *particles) {
        int to_emit = static_cast<int> ( floor((relative_time - last_emit_time)*p_rate ));
        if (to_emit >= 1) {
            last_emit_time = relative_time;
        }
        while (( particles->getMaxLength() - particles->getLength()) >= 1 && to_emit >= 1) {        // 
            particles->add( getStartPos(0), getStartVel(0), relative_time );
            to_emit--;
        }
    }

    vector3d getStartPos(const int &p) {
        /* 
        * This function directly calculates the starting position from the
        * particle number; this comes in handy when resetting particles, if
        * they leave the cube, to the position they started;
        */
        MTRand myran;
        double x = delimiter(0,0) + (delimiter(0,1) - delimiter(0,0)) * myran.rand53(); // myran.doub() is in range 0 to 1
        double y = delimiter(1,0) + (delimiter(1,1) - delimiter(1,0)) * myran.rand53();
        double z = delimiter(2,0) + (delimiter(2,1) - delimiter(2,0)) * myran.rand53();

        return vector3d(x,y,z);
    }

    vector3d getStartVel(const int &p) {
        return vector3d(0,0,p_velocity);
    }
};