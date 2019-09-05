#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Relativity/jw_lighttime.h"

namespace tudat {

namespace jw {

    //! Function to calculate first order relativistic light time correction due to a gravitating point mass.
    double calculate_jw_lighttime(
        const double bodyGravitationalParameter,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const Eigen::Vector6d& centralBodyState,
        const double ppnParameterGamma ) {

            // Calculate Euclidean geometric distances between transmitter, receiver and gravitating body.
            Eigen::Vector3d centralBodyPosition = centralBodyState.segment(0,3);
            double distanceToReceiver = ( receiverPosition - centralBodyPosition ).norm( );
            double distanceToTransmitter = ( transmitterPosition - centralBodyPosition ).norm( );
            double linkEuclideanDistance = ( transmitterPosition - receiverPosition ).norm( );

            // Calculate and return light time correction.
            return ( 1.0 + ppnParameterGamma ) * bodyGravitationalParameter * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * std::log(
                        ( distanceToReceiver + distanceToTransmitter + linkEuclideanDistance ) /
                        ( distanceToReceiver + distanceToTransmitter - linkEuclideanDistance ) );

    }

    Eigen::Matrix< double, 1, 3 > calculate_jw_lighttime_gradient(
        const double bodyGravitationalParameter,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const Eigen::Vector3d& centralBodyPosition,
        const bool evaluateGradientAtReceiver,
        const double ppnParameterGamma ) {
            throw std::runtime_error("jw lighttime gradient not implemented");

            Eigen::Matrix<double, 1, 3> zero_vector;
            zero_vector.setZero();
            return zero_vector;
    }

} // namespace jw

} // namespace tudat
