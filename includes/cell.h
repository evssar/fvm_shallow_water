#ifndef _cell_h
#define _cell_h

#include <array>
#include <Eigen/Dense>

class cell
{

public:
  Eigen::Vector3d r1;
  Eigen::Vector3d r2;
  Eigen::Vector3d r3;
  Eigen::Vector3d r4;
  Eigen::Vector3d r;
  Eigen::Vector3d l;
  Eigen::Vector3d m;
  Eigen::Vector3d n;
  std::array<int,4> neighbours;
  double S;
  double size;
  int ghost;
  Eigen::Vector3d Q;
  Eigen::Vector3d Q_old;
  Eigen::Vector3d R;
  cell(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d,Eigen::Vector3d);
  ~cell(void);
};


#endif
