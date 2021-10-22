// Custom output example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up the number of blocks that the container is divided into.
const int n_x=6,n_y=6,n_z=6;

int main(int argc, char **argv) {

	if(argc < 8){
		printf("Too few input arguments\n");
		exit(0);
	}


	// File containing sites
	char* sites_file = argv[1];

	// Box-geometry
	double x_min = atof(argv[2]), x_max = atof(argv[3]);
	double y_min = atof(argv[4]), y_max = atof(argv[5]);
	double z_min = atof(argv[6]), z_max = atof(argv[7]);

	printf("---Calculating neighbours---\n");

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,20);

	// Import the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	con.import(sites_file);

	// Do a custom output routine to store a variety of face-based
	// statistics. Store the particle ID and position, the number of faces
	// the total face area, the order of each face, the areas of each face,
	// the vertices making up each face, and the neighboring particle (or
	// wall) corresponding to each face.
	con.print_custom("%i %n","neighbours.txt");
}
