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
#include "Particle.h"


////////////////
// Constructors
Particle::Particle() {
    pos = 0;
    vel = 0;
    num = 0;
}

Particle::Particle( int num, const Vector3d &pos, const Vector3d &vel, double spawntime )
{
    this->num = num;
    this->pos = pos;
    this->startPos = pos;
    this->vel = vel;
    this->spawntime = spawntime;
}


////////////////////////////////////////
// Getters and Setters
const Vector3d &Particle::getPos() const
{
    return pos;
}

const Vector3d &Particle::getVel() const
{
    return vel;
}

int Particle::getNum() const
{
    return num;
}

void Particle::setPos( const Vector3d &pos )
{
    this->pos = pos;
}

void Particle::setVel( const Vector3d &vel )
{
    this->vel = vel;
}

void Particle::setNum( int num )
{
    this->num = num;
}

////////////////////////////////////////
// Other stuff
double Particle::spawnTime() const
{
    return spawntime;
}

double Particle::speed() const
{
    return sqrt( dot( vel, vel ) );
}
