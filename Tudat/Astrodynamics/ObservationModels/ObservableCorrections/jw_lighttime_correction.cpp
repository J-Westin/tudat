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



} // namespace observation_models

} // namespace tudat
