#!/usr/bin/bash

sites_file=$1

x_min=$2
x_max=$3
y_min=$4
y_max=$5
z_min=$6
z_max=$7

# use Path to rt preprocessing folder in voro++ library
EXEC=../rt_preprocessing/output_sites

$EXEC $sites_file $x_min $x_max $y_min $y_max $z_min $z_max
