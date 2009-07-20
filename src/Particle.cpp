#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Particle.h"

class Particle {
private:
    vector3d pos;
    vector3d vel;
    int num;
    double spawntime;

public:    
    Particle() {
        pos, vel, num = 0;
    }

    Particle(const int &num, const vector3d &pos, const vector3d &vel, const double &spawntime) {
        this->num = num;
        this->pos = pos;
        this->vel = vel;
        this->spawntime = spawntime;
    }

    //getters and setters
    vector3d &getPos() {
        return pos;
    }
    vector3d &getVel() {
        return vel;
    }
    int &getNum() {
        return num;
    }


    double getSpawnTime() {
        return spawntime;
    }


    // Other functions
    const double getSpeed() {
        return sqrt(dot(vel,vel));
    }
};