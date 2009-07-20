#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "BurgersVortex.h"
#include "Main.h"
#include "Vortex.cpp"
#include <math.h>

class BurgersVortex : public Vortex {
protected:
    // BurgersVortex Specific Properties
    double kappa;                // Determines Strength
    double alpha;                // "stretching"
    double stretching_r;            // 
    double stretching_z;            // 


    virtual vector3d velocityCylinder(const double &r, const double &phi, const double &z) {
        /* check for r == 0, otherwise there's a "divide by zero" and the function
        *  return NaN (inf * 0) */
        return vector3d(    -alpha * stretching_r * r,
            (r==0) ? 0 : ( kappa / r * (1 - exp(-alpha * r*r/(2*fl_nu))) ),
            2*stretching_z*alpha*z    );
    }

    virtual vector3d dudtCylinder(const double &r, const double &phi, const double &z) {
        /* check for r == 0, otherwise there's a "divide by zero" and the function
        *  return NaN (inf * 0) */
        return vector3d(    stretching_r*stretching_r*alpha*alpha*r,
            (r==0) ? 0 : -( (kappa*alpha*exp(-1/2*alpha*r*r/fl_nu))/fl_nu - (kappa*(1-exp(-1/2*alpha*r*r/fl_nu)))/(r*r) )*stretching_r*alpha*r,
            4*stretching_z*stretching_z*alpha*alpha*z    );

    }

public:
    ////////////////////////////////////////
    //            Constructors              //
    ////////////////////////////////////////
    BurgersVortex(const string &parameters, const double &radius, const double &velocity, const double &angle, const double &fl_mu, const double &fl_density, const bool &interpolate, const string &roi) : Vortex(radius, velocity, angle, fl_mu, fl_density, interpolate, roi) {
        alpha = 2.5128624172523393539654752332184326538328336634026 * fl_nu/(radius*radius); //Evaluated LambartW function o_O
        kappa = (-1) * velocity*radius/(-1 + exp(-0.5*alpha*radius*radius/fl_nu));
        cout << "alpha: " << alpha << "\n";
        cout << "kappa: " << kappa << "\n";
        if (parameters == "") {
            this->stretching_r = 1;
            this->stretching_z = 1;
        } else {
            sscanf(parameters.c_str(),"%lf,%lf",&stretching_r,&stretching_z); //"1,1"
        }
    }
};