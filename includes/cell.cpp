/*

Cell class

*/

#include "cell.h"


cell::cell(Eigen::Vector3d r1, Eigen::Vector3d r2, Eigen::Vector3d r3, Eigen::Vector3d r4)
{

  this->r1 = r1;
  this->r2 = r2;
  this->r3 = r3;
  this->r4 = r4;
  this->r=1.0/4.0*(r1+r2+r3+r4);
  this->n = r3-r1;
  this->n = this->n.cross(r4-r2);
  this->S = 0.5 * this->n.norm();
  this->ghost=0;

  this->n = this->n / this->n.norm();
  this->m = this->r1 - this->r;
  this->m = this->m / this->m.norm();
  this->l = this->m.cross(this->n);
  this->l = this->l / this->l.norm();

  this->neighbours[0]=-1;
  this->neighbours[1]=-1;
  this->neighbours[2]=-1;
  this->neighbours[3]=-1;

  Eigen::Vector3d d1,d2,d3,d4;
  double norm_d;

  this->size = 0;

  d1=r2-r1;
  norm_d=d1.norm();
  if(norm_d>this->size){
    this->size=norm_d;
  }
  d2=r3-r2;
  norm_d=d2.norm();
  if(norm_d>this->size){
    this->size=norm_d;
  }
  d3=r4-r3;
  norm_d=d3.norm();
  if(norm_d>this->size){
    this->size=norm_d;
  }
  d4=r1-r4;
  norm_d=d4.norm();
  if(norm_d>this->size){
    this->size=norm_d;
  }
  
}

cell::~cell(void) {}