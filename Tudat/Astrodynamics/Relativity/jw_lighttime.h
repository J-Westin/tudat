#ifndef JW_LIGHTTIME_H
#define JW_LIGHTTIME_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include <cmath>
#include <vector>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat {

namespace gnv {

    double calculate_shapiro_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double ppn_gamma
    );

    double calculate_velocity_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double dt,
        const double ppn_gamma
    );

    double calculate_j2_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double j2r2,
        const Eigen::Vector3d k,
        const double ppn_gamma
    );

    double calculate_second_order_lighttime(
        const double mu_cb,
        const Eigen::Vector6d& state_tx,
        const Eigen::Vector6d& state_rx,
        const Eigen::Vector6d& state_cb,
        const double ppn_gamma
    );

    Eigen::Matrix< double, 1, 3 > calculate_jw_lighttime_gradient(
        const double bodyGravitationalParameter,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const Eigen::Vector3d& centralBodyPosition,
        const bool evaluateGradientAtReceiver,
        const double ppnParameterGamma = 1.0
    );



} // namespace gnv

} // namespace tudat

#endif // JW_LIGHTTIME_H
