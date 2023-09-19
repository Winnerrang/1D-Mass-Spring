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
	Eigen::VectorXd vkappa1, vkappa2, vkappa3, vkappa4, qkappa1, qkappa2, qkappa3, qkappa4;
	double temp_dt;
	q1 = q;
	qdot1 = qdot;

	Eigen::VectorXd f;
	// calcuate kappa1
	force(f, q, qdot);
	vkappa1 = (1 / mass) * f;
	qkappa1 = qdot;
	

	//calculate y^{t+1/2}
	temp_dt = dt / 2;
	qdot2 = qdot + temp_dt * vkappa1;
	q2 = q + temp_dt * qkappa1;

	//calculate kappa2
	force(f, q2, qdot2);
	vkappa2 = (1 / mass) * f;
	qkappa2 = qdot2;

	//calculate y^{t+1}
	temp_dt = dt / 2;
	qdot3 = qdot + temp_dt * vkappa2;
	q3 = q + temp_dt * qkappa2;

	//calculate kappa3
	force(f, q3, qdot3);
	vkappa3 = (1 / mass) * f;
	qkappa3 = qdot3;

	//calculate y^{t+1} (fake one)
	temp_dt = dt;
	qdot4 = qdot + temp_dt * vkappa3;
	q4 = q + temp_dt * qkappa3;

	//calculate kappa4
	force(f, q4, qdot4);
	vkappa4 = (1 / mass) * f;
	qkappa4 = qdot4;

	// calculate real y^{t+1}
	qdot = qdot + dt * (vkappa1 + 2 * vkappa2 + 2 * vkappa3 + vkappa4) / 6;
	q = q + dt * (qkappa1 + 2 * qkappa2 + 2 * qkappa3 + qkappa4) / 6;


}