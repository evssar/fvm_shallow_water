

/*

Function for calculating numerical flux

*/

#include "sw.h"

// Harten's entropy fix
double entropy_fix(double lambda)
{
  double epsilon = 2.0;
  if (abs(lambda) >= epsilon){
    return abs(lambda);
  }else{
    return(lambda*lambda + epsilon*epsilon)/(2*epsilon);
  }
}

//Numerical flux
Eigen::Vector3d flux(Eigen::Vector3d &Ql,Eigen::Vector3d &Qr,Eigen::Vector3d &n,int riemann_solver){

	double hl,ul,vl;
	double hr,ur,vr;
	const double g=9.80665;
	Eigen::Vector3d F,G,delta,w,H;

	//left states
	hl=Ql(0);
	ul=Ql(1)/Ql(0);
	vl=Ql(2)/Ql(0);

	//right states
	hr=Qr(0);
	ur=Qr(1)/Qr(0);
	vr=Qr(2)/Qr(0);

	//difference
	delta=Qr-Ql;

	//x-direction

	F(0) = 0.5*(ul*hl+ur*hr);
	F(1) = 0.5*(ul*ul*hl+0.5*g*hl*hl+ur*ur*hr+0.5*g*hr*hr);
	F(2) = 0.5*(ul*vl*hl+ur*vr*hr);

	//y-direction
	G(0) = 0.5*(vl*hl+vr*hr);
  	G(1) = 0.5*(ul*vl*hl+ur*vr*hr);
  	G(2) = 0.5*(vl*vl*hl+0.5*g*hl*hl+vr*vr*hr+0.5*g*hr*hr);

	//Lax-Friedriech
	//////////////////////////////////////////////////////////////////////
	
	if(riemann_solver==0){

		double lambda_u,lambda_v;

		lambda_u = 0.5*std::abs(ul+ur)+std::sqrt(g*0.5*(hl+hr));
		w=0.5*lambda_u*delta;

		F=F-w;

		lambda_v = 0.5*std::abs(vl+vr)+std::sqrt(g*0.5*(hl+hr));
		w=0.5*lambda_v*delta;

		G=G-w;

	}

	//Roe 	
	//////////////////////////////////////////////////////////////////////

	if(riemann_solver==1){

		double h_bar, u_tilde,v_tilde,c_tilde,lambda1,lambda2,lambda3;
		Eigen::Vector3d r1,r2,r3,I1,I2,I3;

		h_bar = 0.5*(hl+hr);
		u_tilde = (sqrt(hl)*ul+sqrt(hr)*ur)/(sqrt(hl)+sqrt(hr));
		v_tilde = (sqrt(hl)*vl+sqrt(hr)*vr)/(sqrt(hl)+sqrt(hr));
		c_tilde = sqrt(g*h_bar);

		r1(0)=1;
		r1(1)=u_tilde-c_tilde;
		r1(2)=v_tilde;

		r2(0)=0;
		r2(1)=0;
		r2(2)=-1;

		r3(0)=1;
		r3(1)=u_tilde+c_tilde;
		r3(2)=v_tilde;

		I1(0)=0.5+u_tilde/(2.0*c_tilde);
		I1(1)=-1.0/(2.0*c_tilde);
		I1(2)=0;

		I2(0)=v_tilde;
		I2(1)=0;
		I2(2)=-1;

		I3(0)=0.5-u_tilde/(2.0*c_tilde);
		I3(1)=1.0/(2.0*c_tilde);
		I3(2)=0;

		lambda1 = u_tilde-c_tilde;
		lambda2 = u_tilde;
		lambda3 = u_tilde+c_tilde;

		w=entropy_fix(lambda1)*I1.dot(delta)*r1;
		w+=entropy_fix(lambda2)*I2.dot(delta)*r2;
		w+=entropy_fix(lambda3)*I3.dot(delta)*r3;
		w*=0.5;

		F=F-w;

		r1(0)=1;
		r1(1)=u_tilde;
		r1(2)=v_tilde-c_tilde;

		r2(0)=0;
		r2(1)=1;
		r2(2)=0;

		r3(0)=1;
		r3(1)=u_tilde;
		r3(2)=v_tilde+c_tilde;

		I1(0)=0.5+v_tilde/(2.0*c_tilde);
		I1(1)=0;
		I1(2)=-1.0/(2.0*c_tilde);

		I2(0)=-u_tilde;
		I2(1)=1;
		I2(2)=0;

		I3(0)=0.5-v_tilde/(2.0*c_tilde);
		I3(1)=0;
		I3(2)=1.0/(2.0*c_tilde);

		lambda1 = v_tilde-c_tilde;
		lambda2 = v_tilde;
		lambda3 = v_tilde+c_tilde;

		w=entropy_fix(lambda1)*I1.dot(delta)*r1;
		w+=entropy_fix(lambda2)*I2.dot(delta)*r2;
		w+=entropy_fix(lambda3)*I3.dot(delta)*r3;
		w*=0.5;

		G=G-w;
		
	}

	//H=F*nx+G*ny

	H(0)=F(0)*n(0)+G(0)*n(1);
	H(1)=F(1)*n(0)+G(1)*n(1);
	H(2)=F(2)*n(0)+G(2)*n(1);

	return H;
}
