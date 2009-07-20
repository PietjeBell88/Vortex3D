#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Vortex.h"
#include "Main.h"

class Vortex {
protected:
    // Vortex Properties
    double angle;                // Angle of the vortex with the z-axis, clockwise when looking along the y-axis in a rh coordinate system.
    double radius;                // Typical radius of the vortex
    double velocity;            // Velocity at radius

    double fl_mu;                // Fluid Dynamic Viscocity (Pa s).
    double fl_density;            // Fluid Density (kg/m3).
    double fl_nu;                // Fluid Kinematic Viscosity (m2/s).

    // Calculation properties
    bool interpolate;            // true = Interpolate a precalculated grid of the fluid velocity at the position of a particle, false = evaluate fluid velocity functions at particle position

    // Values below are only used when interpolate = true
    TinyVector<int,3> grid;    // X x Y x Z grid of particles
    TinyMatrix<double,3,2> delimiter; //offsets

    double dx, dy, dz;

    vectorField v;                                                            
    vectorField accelfluid;                                                    

    //**************************************//
    //                Functions                //
    //**************************************//
    //

    //    virtual double find_settle_r() = 0;

    // Not to be overwritten
    /*
    inline double v_tot(const double &r) {
    // Calculates the total planar velocity at a certain r in the z=0 plane
    return sqrt(pow2(velocity_phi) + pow2(velocity_r));
    }*/

    // return a matrix that transforms [v_r, v_phi, v_z] to [v_x, v_y, v_z]                            
    TinyMatrix<double, 3, 3> cil2cart(const double &phi) {
        TinyMatrix<double, 3, 3> M;
        M =    cos(phi), -sin(phi), 0,
            sin(phi), cos(phi), 0,
            0, 0, 1;
        return M;
    }

    // Rotation matrix for rotation around the x axis ("folding the y-axis to the z-axis")
    TinyMatrix<double, 3, 3> rotate_x(const double &angle) {
        TinyMatrix<double, 3, 3> M;
        M =    1, 0, 0,
            0, cos(angle), -sin(angle),
            0, sin(angle), cos(angle);
        return M;
    }


    ///////////////////
    // Interpolation //
    ///////////////////
    vector3d interpolate3DCube(const vectorField &v, const vector3d &pos) {
        /*
        This function interpolates the velocity in matrix U to the position of a particle at P.
        It does this by interpolating linearly by assigning weights to each corner of the surrounding box.
        This function requires uniform gridspacing.
        */

        /*
        Start by finding the lower index values of the box surrounding the particle.
        Of course, this step requires that the size of the index does not exceed the integer gridRange.
        */

        int i = static_cast<int> ( floor((pos(0)-delimiter(0,0))/dx) );
        int j = static_cast<int> ( floor((pos(1)-delimiter(1,0))/dy) );
        int k = static_cast<int> ( floor((pos(2)-delimiter(2,0))/dz) );

        /*
        * Calculate the weighting factors for each corner of the box
        * Note: 0 <= x < 1 (and same for y and z)
        * Because the position of the particle can be negative, make it positive first
        * and then do a modulo. Doing modulo on a negavite value can be confusing and inconsistent.
        */

        double x = fmod(pos(0) - delimiter(0,0), dx) / dx; 
        double y = fmod(pos(1) - delimiter(1,0), dy) / dy;
        double z = fmod(pos(2) - delimiter(2,0), dz) / dz;

        //Do a weighted addition of all the corners of the cube surrounding the particle.
        //Please note that v(i,j,k) is indeed a vector3d
        return \
            v(i,j,k) * (1 - x) * (1 - y) * (1 - z) + \
            v(i+1,j,k) * x * (1 - y) * (1 - z) + \
            v(i,j+1,k) * (1 - x) * y * (1 - z) + \
            v(i,j,k+1) * (1 - x) * (1 - y) * z + \
            v(i+1,j,k+1) * x * (1 - y) * z + \
            v(i,j+1,k+1) * (1 - x) * y * z + \
            v(i+1,j+1,k) * x * y * (1 - z) + \
            v(i+1,j+1,k+1) * x * y * z;
    }

    //////////////////////////////////////////
    //        Vortex Velocity and Du/Dt        //
    //////////////////////////////////////////

    virtual vector3d velocityCylinder(const double &r, const double &phi, const double &z) = 0;
    virtual vector3d dudtCylinder(const double &r, const double &phi, const double &z) = 0;

    virtual vector3d velocityCarthesian(const vector3d &pos) {
        const double &x = pos(0);
        const double &y = pos(1);
        const double &z = pos(2);

        double r = sqrt(x*x+y*y);
        double phi = atan2(y,x);

        return product(cil2cart(phi),velocityCylinder(r,phi,z));
    }
    virtual vector3d dudtCarthesian(const vector3d &pos) {
        const double &x = pos(0);
        const double &y = pos(1);
        const double &z = pos(2);

        double r = sqrt(x*x+y*y);
        double phi = atan2(y,x);

        return product(cil2cart(phi),dudtCylinder(r,phi,z));
    }

    virtual vector3d velocityAngle(const vector3d &pos) {
        return product( rotate_x(angle), velocityCarthesian( product( rotate_x(-angle), pos ) ) );
    }
    virtual vector3d dudtAngle(const vector3d &pos) {
        return product( rotate_x(angle), dudtCarthesian( product( rotate_x(-angle), pos ) ) );
    }




    //////////////////////////////
    //        Other Stuff            //
    //////////////////////////////

    /*// This one doesnt have to be implemented, but can be called and used
    virtual double findr_v(double best_r, double stepsize) {
    // finds the first maximum of the velocity in the the z=0 plane
    // by starting at zero, and taking increasingly large steps until 
    // the maximum velocity has been reached.
    double best_v = 0;    
    double minimumstepsize = 1E-100;

    while (true) {
    if (v_tot(best_r + stepsize) > best_v) {
    best_r = abs(best_r + stepsize);
    best_v = v_tot(best_r);
    stepsize = stepsize*2;
    }
    else if (v_tot(best_r - stepsize) > best_v) {
    best_r = abs(best_r - stepsize);            
    best_v = v_tot(best_r);
    stepsize = stepsize*2;
    }
    else if    (stepsize > minimumstepsize) {
    stepsize = stepsize/8;
    }
    else {
    break;
    }
    }
    return abs(best_r);
    }*/

    ////////////////////////////////////////////////////
    //            Initialize Vortex Grid                  //
    ////////////////////////////////////////////////////
    void SetupVortexGrid(vectorField *v, vectorField *accelfluid) {

#pragma omp parallel for
        for (int i = 0; i < grid(0); i++) {
            // What x-coordinate are we at?
            double x = delimiter(0,0) + i * dx;
            for (int j = 0; j < grid(1); j++) {
                // What y-coordinate are we at?
                double y = delimiter(1,0) + j * dy;    
                for (int k = 0; k < grid(2); k++) {
                    // What z-coordinate are we at?
                    double z = delimiter(2,0) + k * dz;    

                    // Calculate and set the velocities for each direction in the vectorfield
                    // Of course you can't have a negative index, so save each at their index+1
                    (*v)(i,j,k) = velocityAngle(vector3d(x,y,z));
                    (*accelfluid)(i,j,k) = dudtAngle(vector3d(x,y,z));
                }
            }
        }
    }

public:
    ////////////////////////////////////////////////////
    //                    Constructor                  //
    ////////////////////////////////////////////////////
    Vortex(const double &radius, const double &velocity, const double &angle, const double &fl_mu, const double &fl_density, const bool &interpolate, const string &roi) {
        this->radius = radius;
        this->velocity = velocity;        
        this->angle = angle;
        this->fl_mu = fl_mu;
        this->fl_density = fl_density;
        this->fl_nu = fl_mu/fl_density;    
        this->interpolate = interpolate;

        //wieoewieoeh
        double x1,x2,y1,y2,z1,z2;
        int X,Y,Z;

        sscanf(roi.c_str(),"[%lf:%d:%lf,%lf:%d:%lf,%lf:%d:%lf]",&x1,&X,&x2,&y1,&Y,&y2,&z1,&Z,&z2); //"[-4:30:4,0:1:0,4:1:4]", also see Emitter
        delimiter = x1*radius, x2*radius,
            y1*radius, y2*radius,
            z1*radius, z2*radius;
        grid = X,Y,Z;

        dx = (delimiter(0,1)-delimiter(0,0))/(grid(0)-1);
        dy = (delimiter(1,1)-delimiter(1,0))/(grid(1)-1);
        dz = (delimiter(2,1)-delimiter(2,0))/(grid(2)-1);

        if (grid(0) <= 1) {
            dx = 0;
        }
        if (grid(1) <= 1) {
            dy = 0;
        }
        if (grid(2) <= 1) {
            dz = 0;
        }
    }

    // Can't do interpolation initialization in the constructor, results in pure function call.
    // As MSDN says: "Do not call overridable methods in constructors".
    void initInterpolate() {
        v.resize(grid(0),grid(1),grid(2));                   // The velocity vector field of the vortex.
        accelfluid.resize(grid(0),grid(1),grid(2));          // The acceleration of the surrounding (vortex)fluid
        SetupVortexGrid(&v, &accelfluid);                    // Initialize vortex grid; Note that the function setup only changes v (so doesnt change options).
    }

    bool outsideBox(const vector3d &pos) { 
        if (pos(0) <= delimiter(0,0) || pos(0) > delimiter(0,1) || pos(1) <= delimiter(1,0) || pos(1) > delimiter(1,1) || pos(2) <= delimiter(2,0) || pos(2) > delimiter(2,1)) {
            return true;
        }
        return false;
    }

    // And just 2 important public functions.. the rest is.. private :>

    vector3d getDuDt(const vector3d &pos) {
        if (interpolate == true) {
            return interpolate3DCube(accelfluid, pos);
        } else {
            return dudtAngle(pos);
        }
    }

    vector3d getVelocity(const vector3d &pos) {
        if (interpolate == true) {
            return interpolate3DCube(v, pos);
        } else {
            return velocityAngle(pos);
        }
    }

    vectorField getVectorField() { 
        return v;
    }

};