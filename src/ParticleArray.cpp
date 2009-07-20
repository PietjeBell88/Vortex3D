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
#include "ParticleArray.h"


////////////////////////////////////////
// Constructor
ParticleArray::ParticleArray(const int &initiallength) {
    particles.resize(initiallength);
    nextIndex = 0;
    length = 0;
}


////////////////////////////////////////
// Getters (and Setters by ref)
Particle &ParticleArray::GetParticle(const int &p) {
    return particles(p);
}

const int ParticleArray::GetLength() {
    return length;
}

const int ParticleArray::GetMaxLength() {
    return particles.size();
}


////////////////////////////////////////
// Add a particle.
void ParticleArray::Add(const Vector3d &pos, const Vector3d &vel, const double &relative_time) {
    particles(length) = Particle(nextIndex, pos, vel, relative_time);
    length++;
    nextIndex++;
}


////////////////////////////////////////
// Remove a particle.
Particle ParticleArray::Remove(const int &p) {
    Particle temp = particles(p);
    particles(p) = particles(length-1);                                                    // [ 1 2 3 4 5 ] at length 5, with particle nr 2 (index 1) outside of the box becomes
    length--;                                                                            // [ 1 5 3 4 5 ] with length 4;
    return temp;
}
