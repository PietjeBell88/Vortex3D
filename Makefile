all: default
# Sources
SRCS  = \
	./external/getopt_pp.cpp ./src/Typedefs.h \
	./src/Particle.cpp ./src/ParticleArray.cpp \
	./src/Vortex/Vortex.cpp ./src/Vortex/BurgersVortex.cpp \
	./src/Output/Output.cpp ./src/Output/ByteOutput.cpp ./src/Output/MatlabOutput.cpp ./src/Output/TextOutput.cpp ./src/Output/TecplotOutput.cpp \
	./src/Emitter/Emitter.cpp ./src/Emitter/GridEmitter.cpp	./src/Emitter/GridOnceEmitter.cpp ./src/Emitter/RandomEmitter.cpp \
	./src/Main.cpp \

CXXFLAGS = -O2 -I ../Include/blitz-0.9 -DNDEBUG

default:
	g++ $(CXXFLAGS) $(SRCS) -o vortex3d

