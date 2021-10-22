# Voro++ makefile
#
# Author : Chris H. Rycroft (Harvard University / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011
# Edited by Elias Udnæs October 21st 2021

# Load the common configuration file
include config_rt.mk

# List of executables
EXECUTABLES= output_sites

# Makefile rules
all: $(EXECUTABLES)

output_sites: output_sites.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o output_sites output_sites.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
