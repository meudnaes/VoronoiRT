// Modified example code

#include "voro++.hh"
using namespace voro;

// Set up the number of blocks that the container is divided into.
// For the uniform sampled grid, it seems like this number should be big,
// larger than 10, or maybe larger than 20
const int n_x=20, n_y=20, n_z=20;

int main(int argc, char **argv) {

	if(argc < 9){
		printf("Too few input arguments\n");
		exit(0);
	}


	// File containing sites
	char* sites_file = argv[1];

	// File to write neighbours
	char* neighbours_file = argv[2];

	// Box-geometry
	const double x_min = atof(argv[3]), x_max = atof(argv[4]);
	const double y_min = atof(argv[5]), y_max = atof(argv[6]);
	const double z_min = atof(argv[7]), z_max = atof(argv[8]);

	printf("---Calculating neighbours---\n");

	// Create a container with the geometry given above, and make it
	// periodic in x and y, and non-periodic in z. Allocate space for
	// eight particles within each computational block.
	container con(x_min, x_max, y_min, y_max, z_min, z_max,
								n_x, n_y, n_z,
								true, true, false,
								8);

	// Import the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	con.import(sites_file);

	// Do a custom output routine to store a variety of face-based
	// statistics. Store the particle ID and position, the number of faces
	// the total face area, the order of each face, the areas of each face,
	// the vertices making up each face, and the neighboring particle (or
	// wall) corresponding to each face.
	con.print_custom("%i %n",neighbours_file);
}
