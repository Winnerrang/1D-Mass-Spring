#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness) {
    dV.resize(1);

    //dU/dq = kq
    dV(0) = stiffness * q(0);
}