#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/jw_lighttime.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/jw_lighttime_correction.h"

#include <iostream>

#ifndef YEET_DEBUG
    #define YEET_DEBUG std::cout << __FILE__ << "\nLine " << __LINE__ << std::endl
#endif // YEET_DEBUG

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

        double evaluationTime; // = TUDAT_NAN;

        double current_shapiro      = 0.0;
        double current_second_order = 0.0;
        double current_velocity     = 0.0;
        double current_j2           = 0.0;

        double total_shapiro      = 0.0;
        double total_second_order = 0.0;
        double total_velocity     = 0.0;
        double total_j2           = 0.0;

        // Iterate over all gravitating bodies.
        for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ ) {
            //evaluationTime = transmissionTime + lightTimeEvaluationContribution_.at( i ) * ( receptionTime - transmissionTime );
            //evaluationTime = compute_evaluation_time(transmitterState, receiverState, transmissionTime, receptionTime);
            evaluationTime =
                compute_evaluation_time(
                    transmitterState,
                    receiverState,
                    perturbingBodyStateFunctions_[ i ](receptionTime),
                    transmissionTime,
                    receptionTime
                );

            if (shapiro_flag || velocity_flag) {
                current_shapiro = jw::calculate_shapiro_lighttime(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState,
                    receiverState,
                    perturbingBodyStateFunctions_[ i ]( evaluationTime ),
                    ppnParameterGamma
                );

                total_shapiro += current_shapiro;

                if (velocity_flag) {
                    current_velocity = jw::calculate_velocity_lighttime(
                        perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                        transmitterState,
                        receiverState,
                        perturbingBodyStateFunctions_[ i ]( evaluationTime ),
                        ppnParameterGamma
                    );

                    total_velocity += current_velocity;


                }
            }

            if (second_order_flag) {
                current_second_order = jw::calculate_second_order_lighttime(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState,
                    receiverState,
                    perturbingBodyStateFunctions_[ i ]( evaluationTime ),
                    ppnParameterGamma
                );

                total_second_order += current_second_order;
                //throw std::runtime_error("Second order lighttime correction not implemented");
            }

            if (j2_flag) {
                current_j2 = jw::calculate_j2_lighttime(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState,
                    receiverState,
                    perturbingBodyStateFunctions_[ i ]( evaluationTime ),
                    perturbing_body_j2_coefficients[i],
                    orientation_axes[i],
                    ppnParameterGamma
                );

                total_j2 += current_j2;
            }
        }



        if (velocity_flag){
            currentTotalLightTimeCorrection_ += total_velocity;
            history_shapiro.push_back(total_shapiro);
            history_velocity.push_back(total_velocity - total_shapiro);
        }
        else if (shapiro_flag) {
            currentTotalLightTimeCorrection_ += total_shapiro;
            history_shapiro.push_back(total_shapiro);
            history_velocity.push_back(0.0);
        }
        else {
            history_shapiro.push_back(0.0);
            history_velocity.push_back(0.0);
        }

        if (second_order_flag) {
            currentTotalLightTimeCorrection_ += total_second_order;
            history_second_order.push_back(total_second_order);
        }
        else {
            history_second_order.push_back(0.0);
        }

        if (j2_flag) {
            currentTotalLightTimeCorrection_ += total_j2;
            history_j2.push_back(total_j2);
        }
        else {
            history_j2.push_back(0.0);
        }

        history_time.push_back(evaluationTime);

        history_total.push_back(currentTotalLightTimeCorrection_);

        return currentTotalLightTimeCorrection_;
    }

    double jw_lighttime_calculator::compute_evaluation_time(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const Eigen::Vector6d& cb_state,
        const double transmissionTime,
        const double receptionTime
    ) {
        //return receptionTime;

        Eigen::Vector3d pos_t, pos_r, pos_cb, vel_cb;
        pos_t  = transmitterState.segment(0,3);
        pos_r  = receiverState.segment(0,3);
        pos_cb = cb_state.segment(0,3);
        vel_cb = cb_state.segment(3,3);

        double c = tudat::physical_constants::SPEED_OF_LIGHT;

        Eigen::Vector3d N_tr, beta_cb;
        N_tr = (pos_r - pos_t).normalized();
        beta_cb = (1.0/c) * vel_cb;

        Eigen::Vector3d g = N_tr - beta_cb;

        double param_2 = (g.dot(pos_r - pos_cb)) / (c*g.norm()*g.norm());

        return std::max(
            transmissionTime,
            receptionTime - std::max(
                0.0,
                param_2
            )
        );
    }


} // namespace observation_models

} // namespace tudat
