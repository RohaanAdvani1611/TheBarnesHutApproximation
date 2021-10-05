// 3 Dimensional Barnes Hut Algorithm :

// Include Libraries
#include <stdio.h>
#include <stdlib.h>

// ADT octree - child pointers are w.r.t 8 octants
typedef struct octree{
    struct octree* child[8];
    int leaves; //no: of leaf-nodes
    
    //COM is positional coordinate of system 
    float com_x;
    float com_y;
    float com_z;
    
    //Mass is mass of Node in System - key value
    float mass;
    
    //Bounds decide the octant the node is added to
    float bound_bot_x;
    float bound_bot_y;
    float bound_bot_z;
    
    float bound_mid_x;
    float bound_mid_y;
    float bound_mid_z;
    
    float bound_top_x;
    float bound_top_y;
    float bound_top_z;
} octree;

//The following are the functions :

/* 1. NewNode Octree:
Allocates memory
Initialize values - bounds are defined by 2 input arguments
Defines Location, mass and size of node's bounding cuboid in 3D space 
Return value is node itself (NULL if memory full)*/

octree* BarnesHut3D_newnode(float x1, float y1, float z1, float x2, float y2, float z2);

/* 2. Destroy Octree:
Frees node and respective child-nodes recursively*/

void BarnesHut3D_destroy(octree *node);

/* 3. Insert Node into Octree:
Takes positional & mass values and adds node to the Octree*/

int BarnesHut3D_add(octree *node, float x, float y, float z, float m);

/* 4. Insert Node into Subtree:
Node is added to appropriate subtree octant w.r.t. bounds*/

int BarnesHut3D_subtree(octree *node, float x, float y, float z, float m);

/* 5. Treecalc using BarnesHut Algorithm:
Calculates overall mass & COM positions of System*/

void BarnesHut3D_treecalc(octree *node);

/* 6. Calculate the force on an item:
Calculates the force on an item with the specified position and mass by the system.  
The Treecalc MUST be done for this function to execute properly.
Parameter x, y, z are the x, y, z positions of the item to be acted upon.
Parameter mass The mass of the item to be acted upon.*/

int BarnesHut3D_force(octree *node, float x, float y, float z, float mass);

/* 7. Display System:
Displays the system we created*/

void BarnesHut3D_traverse(octree *node);