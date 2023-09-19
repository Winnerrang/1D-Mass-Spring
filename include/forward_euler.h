#include <Eigen/Dense>

#ifndef FORWARD_EULER_H
#define FORWARD_EULER_H
//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Forward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Forward Euler time integration

template<typename FORCE>
inline void forward_euler(Eigen::VectorXd& q, Eigen::VectorXd& qdot, double dt, double mass, FORCE& force) {
	// y^t+1 = y^t + dt * 1/m f(y^t)
	auto qdot_old = qdot;
	Eigen::VectorXd f;
	force(f, q, qdot);
	qdot = qdot + dt * (1 / mass) * f;
	q = q + dt * qdot_old;

}
#endif
