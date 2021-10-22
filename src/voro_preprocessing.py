#! /usr/bin/env python

import sys
import numpy as np

############## making voro input ###############

infilepath = "../../voro/rt_preprocessing/output_sites_template.cc"
outfilepath = "../../voro/rt_preprocessing/output_sites_out.cc"

z_min = sys.argv[1]
z_max = sys.argv[2]
x_min = sys.argv[3]
x_max = sys.argv[4]
y_min = sys.argv[5]
y_max = sys.argv[6]

infile = open(infilepath, 'r')
outfile = open(outfilepath, 'w')

for i, line in enumerate(infile):
    if line == "// Box-geometry\n":
        outfile.write("const double z_min={}, z_max={};\n".format(z_min, z_max))
        outfile.write("const double x_min={}, x_max={};\n".format(x_min, x_max))
        outfile.write("const double y_min={}, y_max={};\n".format(y_min, y_max))

    outfile.write(line)

infile.close()
outfile.close()
