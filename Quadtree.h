// 2 Dimensional Barnes Hut Algorithm :

// Include Libraries
#include <stdio.h>
#include <stdlib.h>

//ADT Quadtree - child pointers are w.r.t 4 quadrants
typedef struct quadtree{
    struct quadtree* child[4];
    int leaves; //no: of leaf-nodes
    
    //COM is positional coordinate of system 
    float com_x;
    float com_y;
    //Charge is charge of Node in System - key value
    float charge;
    
    //Bounds decide the quadrant the node is added to
    float bound_bot_x;
    float bound_bot_y;

    float bound_mid_x;
    float bound_mid_y;

    float bound_top_x;
    float bound_top_y;
} quadtree;

//The following are the functions :

/* 1. NewNode Quadtree:
Allocates memory
Initialize values - bounds are defined by 2 input arguments
Defines Location, charge and size of node's bounding square in 2D space 
Return value is node itself (NULL if memory full)*/

quadtree* BarnesHut2D_newnode(float x1, float y1, float x2, float y2);

/* 2. Destroy Quadtree:
Frees node and respective child-nodes recursively*/

void BarnesHut2D_destroy(quadtree *node);

/* 3. Insert Node into Quadtree:
Takes positional & charge values and adds node to the Quadtree*/

int BarnesHut2D_add(quadtree *node, float x, float y, float c);

/* 4. Insert Node into Subtree:
Node is added to appropriate subtree quadrant w.r.t. bounds*/

int BarnesHut2D_subtree(quadtree *node, float x, float y, float c);

/* 5. Treecalc using BarnesHut Algorithm:
Calculates overall charge & COM positions of System*/

void BarnesHut2D_treecalc(quadtree *node);

/* 6. Calculate the force on an item:
Calculates the force on an item with the specified position and charge by the system.  
The Treecalc MUST be done for this function to execute properly.
Parameter x, y are the x, y positions of the item to be acted upon.
Parameter charge The charge of the item to be acted upon.*/

int BarnesHut2D_force(quadtree *node, float x, float y, float charge);

/* 7. Display System:
Displays the system we created*/

void BarnesHut2D_traverse(quadtree *node);