#pragma once
#include "Particle.cpp"
#include <blitz/array.h>

class ParticleArray {
private:
    int length;
    int nextIndex;
    Array<Particle,1> particles;

public:

    ParticleArray(const int &initiallength) {
        particles.resize(initiallength);
        nextIndex = 0;
        length = 0;
    }

    Particle &getParticle(const int &p) {
        return particles(p);
    }

    const int getLength() {
        return length;
    }

    const int getMaxLength() {
        return particles.size();
    }


    void add(const vector3d &pos, const vector3d &vel, const double &relative_time) {
        particles(length) = Particle(nextIndex, pos, vel, relative_time);
        length++;
        nextIndex++;
    }


    Particle remove(const int &p) {
        Particle temp = particles(p);
        particles(p) = particles(length-1);                                                    // [ 1 2 3 4 5 ] at length 5, with particle nr 2 (index 1) outside of the box becomes
        length--;                                                                            // [ 1 5 3 4 5 ] with length 4;
        return temp;
    }
};