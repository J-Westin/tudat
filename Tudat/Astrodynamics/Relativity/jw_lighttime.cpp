#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Relativity/jw_lighttime.h"

#include <cmath>

namespace tudat {

namespace gnv {

    //! Function to calculate first order relativistic light time correction due to a gravitating point mass.
    double calculate_shapiro_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double ppn_gamma ) {

            // Calculate Euclidean geometric distances between transmitter, receiver and gravitating body.
            Eigen::Vector3d pos_tx = state_tx.segment(0,3);
            Eigen::Vector3d pos_rx = state_rx.segment(0,3);
            Eigen::Vector3d pos_cb = state_cb.segment(0,3);

            double r_rx = ( pos_rx - pos_cb ).norm( );
            double r_tx = ( pos_tx - pos_cb ).norm( );
            double R_tr = ( pos_tx - pos_rx ).norm( );

            // Calculate and return light time correction.
            return ( 1.0 + ppn_gamma ) * mu_cb * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT *
                std::log(
                    ( r_rx + r_tx + R_tr ) / ( r_rx + r_tx - R_tr )
                );

    }

    double calculate_velocity_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double ppn_gamma
    ) {
        Eigen::Vector3d pos_tx, pos_rx, pos_cb, vel_cb;
        pos_tx = state_tx.segment(0,3);
        pos_rx = state_rx.segment(0,3);
        pos_cb = state_cb.segment(0,3);
        vel_cb = state_cb.segment(3,3);

        double c = physical_constants::SPEED_OF_LIGHT;

        Eigen::Vector3d R_tr = pos_rx - pos_tx;
        double d_tr = R_tr.norm();

        Eigen::Vector3d N_tr = (1.0/d_tr)*R_tr;

        Eigen::Vector3d beta = (1.0/c)*vel_cb;
        double b = beta.norm();

        double g = 1.0/std::sqrt( 1.0 - b*b );

        Eigen::Vector3d R_ct, R_cr;
        R_ct = pos_tx + (g*g/(1.0+g)) * (beta.dot(pos_tx - pos_cb)) * beta - pos_cb;
        R_cr = pos_rx + (g*g/(1.0+g)) * (beta.dot(pos_rx - pos_cb)) * beta - pos_cb;

        double d_cr = R_cr.norm();

        double factor, frac_t1, frac_t2;

        factor = 2*mu_cb*g*(1.0 - N_tr.dot(beta))/(c*c*c);
        frac_t1 = (R_ct + g*d_tr*beta).norm() + d_cr;
        frac_t2 = g*d_tr*(1.0 - beta.dot(N_tr));

        return factor * std::log( (frac_t1 + frac_t2) / (frac_t1 - frac_t2) );

    }

    double calculate_j2_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double j2r2,
        const Eigen::Vector3d k,
        const double ppn_gamma
    ) {
        Eigen::Vector3d pos_tx, pos_rx, pos_cb;
        pos_tx = state_tx.segment(0,3);
        pos_rx = state_rx.segment(0,3);
        pos_cb = state_cb.segment(0,3);

        Eigen::Vector3d pos_ct, pos_cr, n_ct, n_cr;
        pos_ct = pos_tx - pos_cb;
        pos_cr = pos_rx - pos_cb;
        n_ct = pos_ct.normalized();
        n_cr = pos_cr.normalized();

        double r_ct, r_cr, R_tr;
        r_ct = pos_ct.norm();
        r_cr = pos_cr.norm();
        R_tr = (pos_cr - pos_ct).norm();

        double c = tudat::physical_constants::SPEED_OF_LIGHT;
        double factor =
            mu_cb * j2r2 * R_tr /( c*c*c * r_ct*r_cr * (1.0 + n_ct.dot(n_cr)) );

        double k_nt, k_nr, k_ntr;
        k_nt  = k.dot(n_ct);
        k_nr  = k.dot(n_cr);
        k_ntr = k.dot(n_ct + n_cr);

        double t1, t2, t3;

        t1 = (1.0 - k_nt*k_nt) / r_ct;
        t2 = (1.0 - k_nr*k_nr) / r_cr;
        t3 = (1.0/r_ct + 1.0/r_cr) * k_ntr*k_ntr / (1.0 + n_ct.dot(n_cr));

        return factor * (t1 + t2 - t3);
    }

    double calculate_second_order_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double ppn_gamma
    ) {
        Eigen::Vector3d pos_tx, pos_rx, pos_cb;
        pos_tx = state_tx.segment(0,3);
        pos_rx = state_rx.segment(0,3);
        pos_cb = state_cb.segment(0,3);

        Eigen::Vector3d pos_ct, pos_cr, n_ct, n_cr;
        pos_ct = pos_tx - pos_cb;
        pos_cr = pos_rx - pos_cb;
        n_ct = pos_ct.normalized();
        n_cr = pos_cr.normalized();

        double r_ct, r_cr, R_tr;
        r_ct = pos_ct.norm();
        r_cr = pos_cr.norm();
        R_tr = (pos_cr - pos_ct).norm();

        double c = tudat::physical_constants::SPEED_OF_LIGHT;
        double factor = mu_cb * mu_cb * R_tr / (c*c*c*c*c * r_ct * r_cr);
        double udot = n_ct.dot(n_cr);

        double t1, t2;
        t1 = (7.0 + 8.0*ppn_gamma) * std::acos(udot) / (4.0*std::sqrt(1.0 - udot*udot));
        t2 = (1.0 + ppn_gamma) * (1.0 + ppn_gamma) / (1.0 + udot);

        return factor * (t1 + t2);

    }


    Eigen::Matrix< double, 1, 3 > calculate_jw_lighttime_gradient(
        const double mu_cb,
        const Eigen::Vector3d& pos_tx,
        const Eigen::Vector3d& pos_rx,
        const Eigen::Vector3d& centralBodyPosition,
        const bool evaluateGradientAtReceiver,
        const double ppn_gamma ) {
            throw std::runtime_error("jw lighttime gradient not implemented");

            Eigen::Matrix<double, 1, 3> zero_vector;
            zero_vector.setZero();
            return zero_vector;
    }

} // namespace gnv

} // namespace tudat
