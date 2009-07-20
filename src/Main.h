#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <fstream>
#include <blitz/array.h>
#include <blitz/tinymat.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "Particle.cpp"
#include "ParticleArray.cpp"



using namespace std;
using namespace blitz;

typedef TinyVector<double,3> vector3d;
typedef Array<vector3d,3> vectorField;
typedef Array<double,3> scalarField;

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
