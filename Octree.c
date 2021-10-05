// 3 Dimensional Barnes Hut Algorithm :

// Include Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Octree.h"

// 1. NewNode Octree:
octree* BarnesHut3D_newnode(float x1, float y1, float z1, float x2, float y2, float z2){
    // Allocate node memory
    octree* node = malloc(sizeof(octree));
    if (!node)
        return NULL;
    
    // Bounds are set bottom is lesser of 1 & 2 top is the greater
    node->bound_bot_x = (x1 < x2 ? x1 : x2);
    node->bound_bot_y = (y1 < y2 ? y1 : y2);
    node->bound_bot_z = (z1 < z2 ? z1 : z2);
    node->bound_top_x = (x1 < x2 ? x2 : x1);
    node->bound_top_y = (y1 < y2 ? y2 : y1);
    node->bound_top_z = (z1 < z2 ? z2 : z1);
    node->bound_mid_x = (x1 + x2) / 2;
    node->bound_mid_y = (y1 + y2) / 2;
    node->bound_mid_z = (z1 + z2) / 2;
    
    //Set deafault positional and mass data
    node->com_x = 0;
    node->com_y = 0;
    node->com_z = 0;
    node->mass = 0;
    
    //All octants are default NULL
    for (int i = 0; i < 8; i++)
        node->child[i] = NULL;
    
    node->leaves = 0;
    return node;
}

// 2. Destroy Octree:
void BarnesHut3D_destroy(octree *node) {
    //Base case
    if (!node) 
        return;
    // Use of recursion
    for (int i = 0; i < 8; i++)
        BarnesHut3D_destroy(node->child[i]);
    free(node);
}

// 3. Insert Node into Octree:
// Use of recursion node is added to system using this function 
// Then in system to subtree with subtree function creating new node in appropriate octant
// Then the function is add function is called again as return value to add node in bounds created
int BarnesHut3D_add(octree *node, float x, float y, float z, float m){
    if (!node) 
        return 0;
    // If octree is Empty, data is added to it 
    if (node->leaves == 0) {
        node->com_x = x;
        node->com_y = y;
        node->com_z = z;
        node->mass = m;
    }
    // handle a node that already contains data 
    else {
        /* If single node exists in system, take its position/mass and copy it in the 
        appropriate subtree, no longer making this a leaf */
        if (node->leaves == 1) {
            BarnesHut3D_subtree(node, node->com_x, node->com_y, node->com_z, node->mass);
        }
        /* Since this node is occupied, recursively add the data to the 
         appropriate child node */
        BarnesHut3D_subtree(node, x, y, z, m);
    }
    /* A data point was inserted into this node, therefore the leaf node count 
     must be incremented */
    (node->leaves)++;
    return (node->leaves);
}

// 4. Insert Node into Subtree:
int BarnesHut3D_subtree(octree *node, float x, float y, float z, float m) {
    if (!node) 
        return 0;
    // sub variable stores the octant number
    int sub = 0;
    float min_x, min_y, min_z;
    float max_x, max_y, max_z;
    if (x > node->bound_mid_x) {
        sub += 1;
        min_x = node->bound_mid_x;
        max_x = node->bound_top_x;
    } 
    else {
        min_x = node->bound_bot_x;
        max_x = node->bound_mid_x;
    }
    if (y > node->bound_mid_y) {
        sub += 2;
        min_y = node->bound_mid_y;
        max_y = node->bound_top_y;
    } 
    else {
        min_y = node->bound_bot_y;
        max_y = node->bound_mid_y;
    }
    if (z > node->bound_mid_z) {
        sub += 4;
        min_z = node->bound_mid_z;
        max_z = node->bound_top_z;
    } 
    else {
        min_z = node->bound_bot_z;
        max_z = node->bound_mid_z;
    }
    // If node does'nt exist it is built in appropriate octant with calculated bounds
    if (!node->child[sub])
        node->child[sub] = BarnesHut3D_newnode(min_x, min_y, min_z, max_x, max_y, max_z);
    // (the next line naturally checks for a successful malloc)
    return BarnesHut3D_add(node->child[sub], x, y, z, m);
}

// 5. Treecalc using BarnesHut Algorithm:
void BarnesHut3D_treecalc(octree *node){
    //Base case
    if (!node) 
        return;
    /* If the node is a leaf, then its mass and COM are simply the mass and COM 
    of its data, which were stored when the data was added */
    if (node->leaves == 1)
        return;
    // The mass and COM can be determined from the mass and COM of the node's children 
    node->mass = 0;
    node->com_x = 0;
    node->com_y = 0;
    node->com_z = 0;

    for (int i = 0; i < 8; i++) {
        //If any child node doesn't exist continue
        if (!node->child[i])
            continue;
        
        //Recursively calculate for subtrees
        BarnesHut3D_treecalc(node->child[i]);

        float child_mass = node->child[i]->mass;
        node->mass += child_mass;
        //COM of pt = sumation(child point mass * COM of child point) / pt mass
        node->com_x += child_mass * (node->child[i]->com_x);
        node->com_y += child_mass * (node->child[i]->com_y);
        node->com_z += child_mass * (node->child[i]->com_z);
    }
    //The final COM of BarnesHut is calculated below
    node->com_x /= node->mass;
    node->com_y /= node->mass;
    node->com_z /= node->mass;
}

// 6. Calculate the force on an item:
int BarnesHut3D_force(octree *node, float x, float y, float z, float mass){
    if (!node)
        return 0;
    // variables to store x, y, z directional force on a body 
    //fR is resultant force
    float fx, fy, fz, fR;
    fx = fy = fz = fR = 0;
    /* Calculate the radius between the system's COM and a point */
    float rad = sqrtf(powf(x - node->com_x, 2) + powf(y - node->com_y, 2) + powf(z - node->com_z, 2));
    printf("\n6.Radius of calculaton : %f", rad);
    /* With a radius of 0, the point is in itself.  
    As general physics explodes into burning death at this point,
    we'll return a zero instead */
    if (rad == 0) 
        return 0;//fx, fy, fz = 0
    float r3 = powf(rad, 3);
    //F = G(m1)(m2)(r2-r1)/(r12 ^ 3) - Newton's Force Formula
    fx = 6.672e-11F * mass * node->mass * (x - node->com_x) / r3;
    fy = 6.672e-11F * mass * node->mass * (y - node->com_y) / r3;
    fz = 6.672e-11F * mass * node->mass * (z - node->com_z) / r3;
    printf("\n7.Force on body in x direction : %f", fx);
    printf("\n8.Force on body in y direction : %f", fy);
    printf("\n9.Force on body in z direction : %f", fz);
    
    // Magnitude(fR) = Root(fx^2 + fy^2 + fz^2)
    fR = sqrtf(powf(fx, 2) + powf(fy, 2) + powf(fz, 2));
    printf("\n10.Resultant Force on body  : %f", fR);
    return 1;
}

// 7. Display System:
void BarnesHut3D_traverse(octree *node){
    //Base case
    if (!node)
        return;
    for(int i = 0; i < 8; i++){
        // octant with no child node error handled
        if(!node->child[i])
            continue;
        printf("Latitude:%f", node->child[i]->com_x);
        printf("  Longitude:%f", node->child[i]->com_y);
        printf("  Altitude:%f", node->child[i]->com_z);
        printf("  Mass:%f\n", node->child[i]->mass);
        
        // Use of recursion 
        BarnesHut3D_traverse(node->child[i]);
    }
    return;
}