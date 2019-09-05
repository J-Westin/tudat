#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/jw_lighttime.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/jw_lighttime_correction.h"


namespace tudat {

namespace observation_models {

    //! Function to calculate first order relativistic light time correction due to set of gravitating point masses.
    double jw_lighttime_calculator::calculateLightTimeCorrection(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime )
    {
        // Retrieve ppn parameter gamma.
        double ppnParameterGamma = ppnParameterGammaFunction_( );

        // Initialize correction to zero.
        currentTotalLightTimeCorrection_ = 0.0;

        double evaluationTime = TUDAT_NAN;
        // Iterate over all gravitating bodies.
        for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
        {
            evaluationTime = transmissionTime + lightTimeEvaluationContribution_.at( i ) * ( receptionTime - transmissionTime );
            // Calculate correction due to current body and add to total.
            currentLighTimeCorrectionComponents_[ i ] = jw::calculate_jw_lighttime(
                        perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                        transmitterState.segment( 0, 3 ), receiverState.segment( 0, 3 ),
                        perturbingBodyStateFunctions_[ i ]( evaluationTime ),
                        ppnParameterGamma );
            currentTotalLightTimeCorrection_ += currentLighTimeCorrectionComponents_[ i ];
        }

        return currentTotalLightTimeCorrection_;
    }

    //! Function to compute the partial derivative of the light-time correction w.r.t. link end position
    Eigen::Matrix< double, 3, 1 > jw_lighttime_calculator::
    calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated )
    {
        throw std::runtime_error("jw_lighttime gradient not implemented");
        Eigen::Matrix<double, 3, 1> zero_vector;
        zero_vector.setZero();
        return zero_vector;
    }



} // namespace observation_models

} // namespace tudat
