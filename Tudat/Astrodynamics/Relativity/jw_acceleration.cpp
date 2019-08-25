#include "Tudat/Astrodynamics/Relativity/jw_acceleration.h"

namespace tudat {

namespace jw {

    Eigen::Vector3d schwarzschild_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    ) {
        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa = state_sa.segment(0,3);
        Eigen::Vector3d vel_sa = state_sa.segment(3,3);

        double c = tudat::physical_constants::SPEED_OF_LIGHT;
        double r = pos_sa.norm();

        double mu_c2_r3 = mu_actor / (c*c*r*r*r);

        return mu_c2_r3 * (
            (4.0*mu_actor/r + vel_sa.dot(vel_sa)) * pos_sa +
            4*(pos_sa.dot(vel_sa)) * vel_sa
        );

//        Eigen::Vector3d output_vector;
//        output_vector.setZero();
//        return output_vector;
    }

}

}
