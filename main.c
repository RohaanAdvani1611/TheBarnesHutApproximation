//DSA SY Project - SEM 4
//by - Rohaan Advani - 111903151
// Main/Try File - Runs the program

//Include Libraries and Files
#include <stdio.h>
#include "Octree.h"
#include "Quadtree.h"
#include "Applications.h"

int main(){
    int ch1, ch2;
    printf("CALCULATING N-BODY FORCES USING BARNES HUT ALGORITHM\n");
    printf("\nThe following program will be used to build a system & implement The Barnes Hut Algorithm");
    printf("\n1. Simulate for Newton's Law of Gravitation in 3D");
    printf("\n2. Simulate for Coulomb's Law of Electrostatics in 2D");
    printf("\nChoose one of the above options to proceed : ");
    scanf("%d", &ch1);
    if(ch1 == 1){
        printf("\nChoose test system 1 or 2: ");
        scanf("%d", &ch2);
        if(ch2 == 1)
            satellites("SAT1.txt");
        else if(ch2 == 2)
            satellites("SAT2.txt");
        else{
        printf("\nWrong Option Chosen");
        printf("\n--------------------------------------------------------------------");
        }
        return 0;
    }
    else if(ch1 == 2){
        printf("\nChoose test system 1 or 2: ");
        scanf("%d", &ch2);
        if(ch2 == 1)
            charges("CRG1.txt");
        else if(ch2 == 2)
            charges("CRG2.txt");
        else{
        printf("\nWrong Option Chosen");
        printf("\n--------------------------------------------------------------------");
        }
        return 0;
    }
    else{
        printf("\nWrong Option Chosen");
        printf("\n--------------------------------------------------------------------");
    }
    return 0;
}