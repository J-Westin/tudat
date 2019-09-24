#include "Tudat/Astrodynamics/Relativity/jw_acceleration.h"

#include <iostream>

#include <Eigen/Geometry>

namespace tudat {

namespace jw {

    Eigen::Vector3d kinetic_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    ) {
        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa   = state_sa.segment(0,3);

        Eigen::Vector3d vel_a    = state_actor.segment(3,3);

        const double c = tudat::physical_constants::SPEED_OF_LIGHT;

        double beta_a = vel_a.norm()/c;

        double r = pos_sa.norm();

        return -(2.0*beta_a*beta_a*mu_actor/(r*r*r))*pos_sa;

    }

    Eigen::Vector3d schwarzschild_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    ) {
        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa = state_sa.segment(0,3);
        Eigen::Vector3d vel_sa = state_sa.segment(3,3);

        const double c = tudat::physical_constants::SPEED_OF_LIGHT;
        double r = pos_sa.norm();

        double mu_c2_r3 = mu_actor / (c*c*r*r*r);

        return mu_c2_r3 * (
            (4.0*mu_actor/r - vel_sa.dot(vel_sa)) * pos_sa +
            4*(pos_sa.dot(vel_sa)) * vel_sa
        );

//        Eigen::Vector3d output_vector;
//        output_vector.setZero();
//        return output_vector;
    }

    Eigen::Vector3d cb_velocity_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    ) {
        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa = state_sa.segment(0,3);
        //Eigen::Vector3d vel_sa = state_sa.segment(3,3);
        Eigen::Vector3d vel_a  = state_actor.segment(3,3);

        double r = pos_sa.norm();
        double r_dot_v = pos_sa.dot(vel_a);
        double c = tudat::physical_constants::SPEED_OF_LIGHT;
        double pre_factor = (mu_actor*r_dot_v)/(2.0*c*c*r*r*r);

        return pre_factor * ( ((3.0*r_dot_v)/(r*r)) * pos_sa - 2.0*vel_a);
    }

    Eigen::Vector3d lense_thirring_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        Eigen::Vector3d angular_momentum_actor
    ) {
//        std::cout << __FILE__ << "\n Line " << __LINE__ << std::endl;
//        throw std::runtime_error("Lense Thirring function not implemented");
//        return Eigen::Vector3d::Zero();

        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa = state_sa.segment(0,3);
        Eigen::Vector3d vel_sa = state_sa.segment(3,3);

        double r = pos_sa.norm();

        double G = tudat::physical_constants::GRAVITATIONAL_CONSTANT;
        double c = tudat::physical_constants::SPEED_OF_LIGHT;

        return (2*G/(c*c*r*r*r)) * (
            (3.0/(r*r)) * pos_sa.dot(angular_momentum_actor) * pos_sa.cross(vel_sa) +
            vel_sa.cross(angular_momentum_actor)
        );
    }

    Eigen::Vector3d wavi_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    ) {
//        std::cout << __FILE__ << "\n Line " << __LINE__ << std::endl;
//        throw std::runtime_error("Wavi function not implemented");
//        return Eigen::Vector3d::Zero();

        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector3d pos_sa   = state_sa.segment(0,3);
        Eigen::Vector3d vel_sa   = state_sa.segment(3,3);
        Eigen::Vector3d vel_a    = state_actor.segment(3,3);

        const double c = tudat::physical_constants::SPEED_OF_LIGHT;

        double r = pos_sa.norm();

        return 4*mu_actor / (c*c*r*r*r) * (
            vel_a.dot(vel_sa) * pos_sa -
            pos_sa.dot(vel_sa) * vel_a
        );
    }

    Eigen::Vector3d de_sitter_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        Eigen::Vector6d state_primary,
        double mu_primary
    ) {
//        std::cout << __FILE__ << "\n Line " << __LINE__ << std::endl;
//        throw std::runtime_error("de Sitter function not implemented");
//        return Eigen::Vector3d::Zero();
        Eigen::Vector6d state_sa = state_subject - state_actor;
        Eigen::Vector6d state_ap = state_actor - state_primary;

        Eigen::Vector3d v_sa = state_sa.segment(3,3);
        Eigen::Vector3d r_ap = state_ap.segment(0,3);
        Eigen::Vector3d v_ap = state_ap.segment(3,3);

        double c = tudat::physical_constants::SPEED_OF_LIGHT;

        //double r = r_sa.norm();
        double R = r_ap.norm();

        double pre_factor = 3.0*mu_primary/(c*c*R*R*R);

        return pre_factor * (
            (v_ap.cross(r_ap)).cross(v_sa)
        );
    }

}

}
