// 2 Dimensional Barnes Hut Algorithm :

// Include Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Quadtree.h"

// 1. NewNode Quadtree:
quadtree* BarnesHut2D_newnode(float x1, float y1, float x2, float y2){
    //Allocate node memory
    quadtree* node = malloc(sizeof(quadtree));
    if (!node)
        return NULL;
    
    //Bounds are set bottom is lesser of 1 & 2 top is the greater
    node->bound_bot_x = (x1 < x2 ? x1 : x2);
    node->bound_bot_y = (y1 < y2 ? y1 : y2);
    node->bound_top_x = (x1 < x2 ? x2 : x1);
    node->bound_top_y = (y1 < y2 ? y2 : y1);
    node->bound_mid_x = (x1 + x2) / 2;
    node->bound_mid_y = (y1 + y2) / 2;
    
    //Set deafault positional and charge data
    node->com_x = 0;
    node->com_y = 0;
    node->charge = 0;
    
    //All quarants are default NULL
    for (int i = 0; i < 4; i++)
        node->child[i] = NULL;

    node->leaves = 0;
    return node;
}

// 2. Destroy Quadtree:
void BarnesHut2D_destroy(quadtree *node){
    //Base Case
    if (!node) 
        return;
    //Use of recursion
    for (int i = 0; i < 4; i++)
        BarnesHut2D_destroy(node->child[i]);
    free(node);
}

// 3. Insert Node into Quadtree:
// Use of recursion node is added to system using this function 
// Then in system to subtree with subtree function creating new node in appropriate quadrant
// Then the function is add function is called again as return value to add node in bounds created
int BarnesHut2D_add(quadtree *node, float x, float y, float c){
    if (!node) 
        return 0;
    // If quadtree is empty, data is added to it 
    if (node->leaves == 0) {
        node->com_x = x;
        node->com_y = y;
        node->charge = c;
    }
    // handle a node that already contains data 
    else {
        /* If single node exists in system, take its position/charge and copy it in the 
        appropriate subtree, no longer making this a leaf */
        if (node->leaves == 1) {
            BarnesHut2D_subtree(node, node->com_x, node->com_y, node->charge);
        }
        /* Since this node is occupied, recursively add the data to the 
         appropriate child node */
        BarnesHut2D_subtree(node, x, y, c);
    }
    /* A data point was inserted into this node, therefore the element count 
     must be incremented */
    (node->leaves)++;
    return (node->leaves);
}

// 4. Insert Node into Subtree:
int BarnesHut2D_subtree(quadtree *node, float x, float y, float c){
    if (!node) 
        return 0;
    // sub variable stores the quadrant number
    int sub = 0;
    float min_x, min_y;
    float max_x, max_y;
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
    // If node does'nt exist it is built in appropriate quadrant with calculated bounds
    if (!node->child[sub])
        node->child[sub] = BarnesHut2D_newnode(min_x, min_y, max_x, max_y);
    // (the next line naturally checks for a successful malloc)
    return BarnesHut2D_add(node->child[sub], x, y, c);
}

// 5. Treecalc using BarnesHut Algorithm:
void BarnesHut2D_treecalc(quadtree *node){
    //Base case
    if (!node) 
        return;
    /* If the node is a leaf, then its charge and COM are simply the charge and COM 
    of its data, which were stored when the data was added */
    if (node->leaves == 1)
        return;
    // The mass and COM can be determined from the mass and COM of the node's children 
    node->charge = 0;
    node->com_x = 0;
    node->com_y = 0;

    for (int i = 0; i < 4; i++) {
        //If any child node doesn't exist continue
        if (!node->child[i])
            continue;
        //Recursively calculate for subtrees
        BarnesHut2D_treecalc(node->child[i]);

        float child_charge = node->child[i]->charge;
        node->charge += child_charge;
        //COM of pt = sumation(child point charge * COM of child point) / pt charge
        node->com_x += child_charge * (node->child[i]->com_x);
        node->com_y += child_charge * (node->child[i]->com_y);
    }
    //The final COM of BarnesHut is calculated below
    node->com_x /= node->charge;
    node->com_y /= node->charge;
}

// 6. Calculate the force on an item:
int BarnesHut2D_force(quadtree *node, float x, float y, float charge){
    if (!node)
        return 0;
    // variables to store x, y, z directional force on a body 
    //fR is resultant force
    float fx, fy, fR;
    fx = fy = fR = 0;
    /* Calculate the radius between the system's COM and a point */
    float rad = sqrtf(powf(x - node->com_x, 2) + powf(y - node->com_y, 2));
    printf("\n5.Radius of calculaton : %f", rad);
    /* With a radius of 0, the point is in itself.  
    As general physics explodes into burning death at this point,
    we'll return a zero instead */
    if (rad == 0) 
        return 0;//fx, fy, fz = 0
    float r3 = powf(rad, 3);
    //printf("\n6.%f", r3);
    //F = K(q1)(q2)(r2-r1)/(r12 ^ 3) - Coulumbs's Force Formula
    fx = 9e9F * charge * node->charge * (x - node->com_x) / r3;
    fy = 9e9F * charge * node->charge * (y - node->com_y) / r3;
    printf("\n6.Force on body in x direction : %f", fx);
    printf("\n7.Force on body in y direction : %f", fy);
    
    // Magnitude(fR) = Root(fx^2 + fy^2)
    fR = sqrtf(powf(fx, 2) + powf(fy, 2));
    printf("\n8.Resultant Force on body  : %f", fR);
    return 1;
}

// 7. Display System:
void BarnesHut2D_traverse(quadtree *node){
    //Base case
    if (!node)
        return;
    for(int i = 0; i < 4; i++){
        // Quadrant with no child node error handled
        if(!node->child[i])
            continue;
        printf("X:%f", node->child[i]->com_x);
        printf("  Y:%f", node->child[i]->com_y);
        printf("  Charge:%f\n", node->child[i]->charge);
        
        // Use of recursion 
        BarnesHut2D_traverse(node->child[i]);
    }
    return;
}