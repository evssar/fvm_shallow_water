/*

Function for finding cell neighbours

*/

#include "sw.h"

// Check if two edges are adjacent
bool are_edges_adjacent(Eigen::Vector3d &a1, Eigen::Vector3d &a2, Eigen::Vector3d &b1, Eigen::Vector3d &b2)
{
    Eigen::Vector3d c1, c2;
    c1 = a1 - b2;
    c2 = a2 - b1;

    if (c1.norm() < 1e-6 && c2.norm() < 1e-6)
    {
        return true;
    }

    return false;
}
// Find cell neighbours
void find_cell_neighbours(sw &sw)
{
         for (int k = 0; k < sw.N_cells; k++)
        {
            sw.cells[k].neighbours[0] = -1.0;
            sw.cells[k].neighbours[1] = -1.0;
            sw.cells[k].neighbours[2] = -1.0;
            sw.cells[k].neighbours[3] = -1.0;
        }

        for (int k = 0; k < sw.N_cells; k++)
        {

            Eigen::Vector3d rk,rm,deltar,r1k,r2k,r3k,r4k,r1m,r2m,r3m,r4m;
            
            rk = sw.cells[k].r;
            
            for(int m=0;m<sw.N_cells;m++){

                if(m!=k){

                    rm=sw.cells[m].r;

                    deltar=rm-rk;

                    if(deltar.norm()<5.0*sw.cells[k].size){
                        
                        //Cell m is within 5 distances from cell k
                        
                        r1k = sw.cells[k].r1;
                        r2k = sw.cells[k].r2;
                        r3k = sw.cells[k].r3;
                        r4k = sw.cells[k].r4;

                        r1m = sw.cells[m].r1;
                        r2m = sw.cells[m].r2;
                        r3m = sw.cells[m].r3;
                        r4m = sw.cells[m].r4;

                        // Side 1
                        //////////////
                        if(are_edges_adjacent(r1k, r2k, r1m, r2m)){
                                sw.cells[k].neighbours[0] = m;
                        }
                        if(are_edges_adjacent(r1k, r2k, r2m, r3m)){
                                sw.cells[k].neighbours[0] = m;
                        }
                        if(are_edges_adjacent(r1k, r2k, r3m, r4m)){
                                sw.cells[k].neighbours[0] = m;
                        }
                        if(are_edges_adjacent(r1k, r2k, r4m, r1m)){
                                sw.cells[k].neighbours[0] = m;
                        }
                        // Side 2
                        //////////////
                        if(are_edges_adjacent(r2k, r3k, r1m, r2m)){
                                sw.cells[k].neighbours[1] = m;
                        }
                        if(are_edges_adjacent(r2k, r3k, r2m, r3m)){
                                sw.cells[k].neighbours[1] = m;
                        }
                        if(are_edges_adjacent(r2k, r3k, r3m, r4m)){
                                sw.cells[k].neighbours[1] = m;
                        }
                        if(are_edges_adjacent(r2k, r3k, r4m, r1m)){
                                sw.cells[k].neighbours[1] = m;
                        }
                        // Side 3
                        //////////////
                        if(are_edges_adjacent(r3k, r4k, r1m, r2m)){
                                sw.cells[k].neighbours[2] = m;
                        }
                        if(are_edges_adjacent(r3k, r4k, r2m, r3m)){
                                sw.cells[k].neighbours[2] = m;
                        }
                        if(are_edges_adjacent(r3k, r4k, r3m, r4m)){
                                sw.cells[k].neighbours[2] = m;
                        }
                        if(are_edges_adjacent(r3k, r4k, r4m, r1m)){
                                sw.cells[k].neighbours[2] = m;
                        }
                        //Side 4
                        ///////////////
                        if(are_edges_adjacent(r4k, r1k, r1m, r2m)){
                                sw.cells[k].neighbours[3] = m;
                        }
                        if(are_edges_adjacent(r4k, r1k, r2m, r3m)){
                                sw.cells[k].neighbours[3] = m;
                        }
                        if(are_edges_adjacent(r4k, r1k, r3m, r4m)){
                                sw.cells[k].neighbours[3] = m;
                        }
                        if(are_edges_adjacent(r4k, r1k, r4m, r1m)){
                                sw.cells[k].neighbours[3] = m;
                        }
                    }
                }
            }
        }
        

}
