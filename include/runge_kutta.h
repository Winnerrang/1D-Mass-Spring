#include <Eigen/Dense>
#include "forward_euler.h"
//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {


	Eigen::VectorXd f1, f2, f3, f4;
	Eigen::VectorXd q1, q2, q3, q4, qdot1, qdot2, qdot3, qdot4;
	q1 = q;
	qdot1 = qdot;

	q2 = q1;
	qdot2 = qdot1;
	forward_euler(q2, qdot2, dt / 2, mass, force);

	q3 = q2;
	qdot3 = qdot2;
	forward_euler(q3, qdot3, dt / 2, mass, force);

	q4 = q3;
	qdot4 = qdot3;
	forward_euler(q4, qdot4, dt, mass, force);

	
	force(f1, q1, qdot1);
	force(f2, q2, qdot2);
	force(f3, q3, qdot3);
	force(f4, q4, qdot4);

	qdot += dt / 6.0 * (f1 + 2.0 * f2 + 2.0 * f3 + f4)/mass;
	q += dt / 6.0 * (qdot1 + 2.0 * qdot2 + 2.0 * qdot3 + qdot4);

}