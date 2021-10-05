//Applications of Barnes Hut Algorithm

//Include Libraries
#include <stdio.h>
#include <stdlib.h>

//Functions are declared here:

/*1. satellites() :
Runs the satellite system created in txt file as input
and builds the data structures with respect to 3D calculations in program
calculates Newton's Force of a system at a point*/

void satellites(char *filename);

/*2. charges() :
Runs the charge system created in txt file as input
and builds the data structures with respect to 2D calculations in program
calculates Coulumbs's Force of a system at a point*/

void charges(char *filename);

/*3. convertlat() & convertlon() 
Convert Latitude & Longitude to Km*/

//Latitude is y(vertical direction)
float convertlat(float y);

//Longitude is x(horizontal direction)
float convertlon(float x);