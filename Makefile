CXX=g++
FC=gfortran

CXXFLAGS=
FFLAGS=-std=legacy -fno-align-commons -fno-underscoring -finit-local-zero
LDFLAGS=
LDLIBS=-lstdc++

CROPS=Crop\ Source
SOILS=Soil\ Source

all: maizsim

crop:
	$(CXX) -c $(CXXFLAGS) $(CROPS)/*.cpp

soil:
	$(FC) -c $(FFLAGS) $(SOILS)/*.for $(SOILS)/*.FOR

maizsim: crop soil
	$(FC) $(LDFLAGS) $(LDLIBS) *.o -o maizsim

clean:
	rm -f *.o
