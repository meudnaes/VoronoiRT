#!/usr/bin/bash

sites_file=$1
neighbours_file=$2

x_min=$3
x_max=$4
y_min=$5
y_max=$6
z_min=$7
z_max=$8

# use Path to rt preprocessing folder in voro++ library
EXEC=../rt_preprocessing/output_sites

$EXEC $sites_file $neighbours_file $x_min $x_max $y_min $y_max $z_min $z_max
