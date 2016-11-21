# makefile for gray2vec

CC=gcc
CXX=g++

ARCHFLAGS = -O3 -mtune=native -march=native

CFLAGS = $(ARCHFLAGS)
CFLAGS_GDAL  = `gdal-config --cflags`
LDFLAGS_CIMG =
LDFLAGS_GDAL = `gdal-config --libs`

CXXFLAGS = $(CFLAGS)

# ---------------------------------------

all: gray2vec

install: all
	cp gray2vec /usr/local/bin/

clean:
	rm -f *.o
	rm -f gray2vec


gray2vec: gray2vec.o Gray2Vec_Grid.o gdal_polygonize_mod.o gdalrasterpolygonenumerator.o
	$(CXX) $(LDFLAGS_CIMG) $(LDFLAGS_GDAL) gray2vec.o Gray2Vec_Grid.o gdal_polygonize_mod.o gdalrasterpolygonenumerator.o -o gray2vec -L.


gray2vec.o: gray2vec.cpp Gray2Vec_Grid.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS_GDAL) -o gray2vec.o gray2vec.cpp

Gray2Vec_Grid.o: Gray2Vec_Grid.cpp Gray2Vec_Grid.h gdal_polygonize_mod.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS_GDAL) -o Gray2Vec_Grid.o Gray2Vec_Grid.cpp

gdal_polygonize_mod.o: gdal_polygonize_mod.cpp gdal_polygonize_mod.h
	$(CXX) -c $(CXXFLAGS) $(CFLAGS_GDAL) -o gdal_polygonize_mod.o gdal_polygonize_mod.cpp

# this is actually only needed for GDAL 1.x
gdalrasterpolygonenumerator.o: gdalrasterpolygonenumerator.cpp
	$(CXX) -c $(CXXFLAGS) $(CFLAGS_GDAL) -o gdalrasterpolygonenumerator.o gdalrasterpolygonenumerator.cpp

