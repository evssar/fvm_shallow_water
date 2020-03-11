#include <iostream>
#include <string>
#include <chrono>

#include "sw.h"

int main(int argc, char * argv [] )
{
	if(argc<2){
		std::cout<<"Usage: sw_solver <case>"<<std::endl;
		std::exit(-1);
	}

	sw sw;
	sw.case_name=argv[1];

	std::cout<<"Reading config...";
	read_config(sw);
	std::cout<<"Done!"<<std::endl;

	std::cout<<"Using";
	if(sw.high_order){
		std::cout<<" high order";
	}else{
		std::cout<<" low order";
	}
	std::cout<<" approximation in conjunction with the ";
	if(sw.riemann_solver){
		std::cout<<"Roe Riemann solver.";
	}else{
		std::cout<<"Lax-Friedriech Riemann solver.";
	}

	std::cout<<std::endl;

	std::cout<<"Reading grid...";
	read_grid(sw);
	std::cout<<"Done!"<<std::endl;
	std::cout<<"N_cells="<<sw.N_cells<<std::endl;
	
	std::cout<<"Finding cell neighbours...";
	find_cell_neighbours(sw);
	std::cout<<"Done!"<<std::endl;

	//Boundary cell
	//Note: If a side of the cell does not have a neighbour the cell is a boundary cell
	for(int i=0;i<sw.N_cells;i++){
		for(int j=0;j<sw.cells[i].neighbours.size();j++){
			if(sw.cells[i].neighbours[j]==-1){
				sw.cells[i].ghost=1;
				break;
			}
		}
	}

	//Initial conditions
	for(int i=0;i<sw.N_cells;i++){
		if(sw.cells[i].r[0]<=-0.3){
			sw.cells[i].Q_old(0)=1.5;
			sw.cells[i].Q_old(1)=0;
			sw.cells[i].Q_old(2)=0;
		}else{
			sw.cells[i].Q_old(0)=1;
			sw.cells[i].Q_old(1)=0;
			sw.cells[i].Q_old(2)=0;	
		}
	}

	std::cout<<"Performing "<<sw.N_timesteps<<" explicit steps..."<<std::endl;

	double t=0;

	//Run timeloop
	for(int n=0;n<sw.N_timesteps;n++){

		std::cout<<"Step "<<n<<", t="<<t<<std::endl;

		//Loop over each cell
		//Note: This is not efficient, in reality one would loop over every edge. This 
		//version of the code is used for demonstration purposes so efficiency is not sought.
		for(int i=0;i<sw.N_cells;i++){

			//If cell is not a boundary cell
			if(sw.cells[i].ghost==0){
				
				int neighbour;
				double norm_d;
				Eigen::Vector3d sum_Fij;
				Eigen::Vector3d Ql,Qr,Fij,d,n;

				sum_Fij(0)=0;
				sum_Fij(1)=0;
				sum_Fij(2)=0;

				//Loop over each cell side
				for(int j=0;j<4;j++){

					neighbour=sw.cells[i].neighbours[j];

					if(j==0){
						d=sw.cells[i].r2-sw.cells[i].r1;
						Ql=sw.cells[neighbour].Q_old;
						Qr=sw.cells[i].Q_old;
					}

					if(j==1){
						d=sw.cells[i].r3-sw.cells[i].r2;
						Ql=sw.cells[i].Q_old;
						Qr=sw.cells[neighbour].Q_old;
					}

					if(j==2){
						d=sw.cells[i].r4-sw.cells[i].r3;
						Ql=sw.cells[i].Q_old;
						Qr=sw.cells[neighbour].Q_old;
					}

					if(j==3){
						d=sw.cells[i].r1-sw.cells[i].r4;
						Ql=sw.cells[neighbour].Q_old;
						Qr=sw.cells[i].Q_old;
						
					}

					norm_d=d.norm();
					n=d.cross(sw.cells[i].n);
					n/=n.norm();
					
					if(sw.high_order==1){
						//Reconstruct Ql,Qr states to obtain a higher order solution.
						//Ql=Ql+phi*grad_Ql*delta_rl;
						//Qr=Qr+phi*grad_Qr*delta_rr;
					}

					Fij=flux(Ql,Qr,n,sw.riemann_solver); 
					
					sum_Fij+=Fij*norm_d;

				}
				sw.cells[i].R=sum_Fij;	
			}
		}
		
		//Advance solution in time by using Euler integration
		//Note: Runge-Kutta, Adams-Bashforth etc for higher order.
		for(int i=0;i<sw.N_cells;i++){
			sw.cells[i].Q=sw.cells[i].Q_old-(sw.dt/sw.cells[i].S)*sw.cells[i].R;
			sw.cells[i].Q_old=sw.cells[i].Q;
		}

		//Set boundary conditions
		for(int i=0;i<sw.N_cells;i++){
			if(sw.cells[i].ghost==1){
				int neighbour;
				//Along x, min y
				if(sw.cells[i].neighbours[0]==-1){
					neighbour=sw.cells[i].neighbours[2];
					sw.cells[i].Q_old=sw.cells[neighbour].Q_old;
					sw.cells[i].Q_old[2]*=-1.0;
				}
				//Along x, max y
				if(sw.cells[i].neighbours[2]==-1){
					neighbour=sw.cells[i].neighbours[0];
					sw.cells[i].Q_old=sw.cells[neighbour].Q_old;
					sw.cells[i].Q_old[2]*=-1.0;
				}
				//Along y, max x
				if(sw.cells[i].neighbours[1]==-1){
					neighbour=sw.cells[i].neighbours[3];
					sw.cells[i].Q_old=sw.cells[neighbour].Q_old;
					sw.cells[i].Q_old[1]*=-1.0;
				}
				//Along y, min x
				if(sw.cells[i].neighbours[3]==-1){
					neighbour=sw.cells[i].neighbours[1];
					sw.cells[i].Q_old=sw.cells[neighbour].Q_old;
					sw.cells[i].Q_old[1]*=-1.0;
				}
			}
		}
		t+=sw.dt;
		write_results(sw,n);
	}

	std::cout<<"Done!"<<std::endl;

	
	return 0;
}
