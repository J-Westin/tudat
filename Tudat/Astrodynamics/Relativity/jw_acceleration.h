#ifndef JW_ACCELERATION_H
#define JW_ACCELERATION_H

#include <memory>
#include <functional>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Basics/basicTypedefs.h"

//#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"

namespace tudat {

namespace jw {
    Eigen::Vector3d kinetic_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    Eigen::Vector3d schwarzschild_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor,
        double ppn_gamma
    );

    Eigen::Vector3d cb_velocity_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    Eigen::Vector3d lense_thirring_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        Eigen::Vector3d angular_momentum_actor
    );

    Eigen::Vector3d wavi_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    Eigen::Vector3d de_sitter_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        Eigen::Vector6d state_primary,
        double mu_primary
    );



    class jw_acceleration : public basic_astrodynamics::AccelerationModel < Eigen::Vector3d > {
    private:
        std::function< Eigen::Vector6d() > state_function_subject;
        std::function< Eigen::Vector6d() > state_function_actor;
        std::function< double() > mu_function_actor;

        Eigen::Vector3d angular_momentum_actor;

        bool newton_flag, kinetic_flag, schwarzschild_flag,
             de_sitter_flag, cb_velocity_flag, cb_acceleration_flag,
             lense_thirring_flag, wavi_flag;

        std::function< double() > ppn_gamma_function;

        //std::shared_ptr<jw_acceleration_settings> acceleration_settings;

        Eigen::Vector3d current_acceleration;
        Eigen::Vector3d current_schwarzschild_acceleration;
        Eigen::Vector3d current_cb_velocity_acceleration;
        Eigen::Vector3d current_kinetic_acceleration;
        Eigen::Vector3d current_lense_thirring_acceleration;
        Eigen::Vector3d current_wavi_acceleration;
        Eigen::Vector3d current_de_sitter_acceleration;

        std::function< Eigen::Vector6d() > state_function_primary;
        std::function< double() > mu_function_primary;

    public:

        //using namespace tudat::simulation_setup;

        jw_acceleration(
            std::function< Eigen::Vector6d() > state_function_subject,
            std::function< Eigen::Vector6d() > state_function_actor,
            std::function< double() > mu_function_actor,
            Eigen::Vector3d angular_momentum_actor,
            bool newton_flag,
            bool kinetic_flag,
            bool schwarzschild_flag,
            bool de_sitter_flag,
            bool cb_velocity_flag,
            bool cb_acceleration_flag,
            bool lense_thirring_flag,
            bool wavi_flag,
            std::function< double() > ppn_gamma_function = [ ]( ){ return 1.0; }
            //std::shared_ptr<jw_acceleration_settings> acceleration_settings
        ):
            state_function_subject(state_function_subject),
            state_function_actor(state_function_actor),
            mu_function_actor(mu_function_actor),
            angular_momentum_actor(angular_momentum_actor),
            newton_flag(newton_flag),
            kinetic_flag(kinetic_flag),
            schwarzschild_flag(schwarzschild_flag),
            de_sitter_flag(de_sitter_flag),
            cb_velocity_flag(cb_velocity_flag),
            cb_acceleration_flag(cb_acceleration_flag),
            lense_thirring_flag(lense_thirring_flag),
            wavi_flag(wavi_flag),
            ppn_gamma_function(ppn_gamma_function)
        { }


        jw_acceleration(
            std::function< Eigen::Vector6d() > state_function_subject,
            std::function< Eigen::Vector6d() > state_function_actor,
            std::function< Eigen::Vector6d() > state_function_primary_in,
            std::function< double() > mu_function_actor,
            std::function< double() > mu_function_primary_in,
            Eigen::Vector3d angular_momentum_actor,
            bool newton_flag,
            bool kinetic_flag,
            bool schwarzschild_flag,
            bool de_sitter_flag,
            bool cb_velocity_flag,
            bool cb_acceleration_flag,
            bool lense_thirring_flag,
            bool wavi_flag
            //std::shared_ptr<jw_acceleration_settings> acceleration_settings
        ):
            jw_acceleration (
                state_function_subject,
                state_function_actor,
                mu_function_actor,
                angular_momentum_actor,
                newton_flag,
                kinetic_flag,
                schwarzschild_flag,
                de_sitter_flag,
                cb_velocity_flag,
                cb_acceleration_flag,
                lense_thirring_flag,
                wavi_flag
            )
//            state_function_primary(state_function_primary),
//            mu_function_primary(mu_function_primary)
        {
            state_function_primary = state_function_primary_in;
            mu_function_primary = mu_function_primary_in;
        }


        ~jw_acceleration() { }

        Eigen::Vector3d getAcceleration () {
            return current_acceleration;
        }

        Eigen::Vector3d get_schwarzschild_acceleration() {
            return current_schwarzschild_acceleration;
        }

        Eigen::Vector3d get_cb_velocity_acceleration() {
            return current_cb_velocity_acceleration;
        }

        Eigen::Vector3d get_lense_thirring_acceleration() {
            return current_lense_thirring_acceleration;
        }

        Eigen::Vector3d get_wavi_acceleration() {
            return current_wavi_acceleration;
        }

        Eigen::Matrix<double, 8, 1> get_acceleration_components() {
            Eigen::Matrix<double, 8, 1> component_vector;
            component_vector <<
                0.0, current_kinetic_acceleration.norm(), current_schwarzschild_acceleration.norm(),
                current_de_sitter_acceleration.norm(), current_cb_velocity_acceleration.norm(), 0.0,
                current_lense_thirring_acceleration.norm(), current_wavi_acceleration.norm();

            return component_vector;
        }

        std::function< Eigen::Vector6d() > get_state_function_subject () {
            return state_function_subject;
        }

        std::function< Eigen::Vector6d() > get_state_function_actor () {
            return state_function_actor;
        }

        std::function< double() > get_mu_function_actor () {
            return mu_function_actor;
        }

        std::function< double() > get_ppn_gamma_function () {
            return ppn_gamma_function;
        }

        void updateMembers( const double currentTime = TUDAT_NAN ) {
            Eigen::Vector6d state_subject = state_function_subject();
            Eigen::Vector6d state_actor   = state_function_actor();
            Eigen::Vector6d state_primary = state_function_primary();

            double mu_actor   = mu_function_actor();
            double mu_primary = mu_function_primary();

            double ppn_gamma = ppn_gamma_function();

            current_acceleration.setZero();
            current_kinetic_acceleration.setZero();
            current_schwarzschild_acceleration.setZero();
            current_cb_velocity_acceleration.setZero();
            current_lense_thirring_acceleration.setZero();
            current_wavi_acceleration.setZero();
            current_de_sitter_acceleration.setZero();

            if (kinetic_flag) {
                current_kinetic_acceleration =
                    kinetic_acceleration (
                        state_subject,
                        state_actor,
                        mu_actor
                    );

                current_acceleration += current_kinetic_acceleration;
            }

            if (schwarzschild_flag) {
                current_schwarzschild_acceleration =
                    schwarzschild_acceleration(
                        state_subject,
                        state_actor,
                        mu_actor,
                        ppn_gamma
                    );

                current_acceleration += current_schwarzschild_acceleration;
            }

            if (cb_velocity_flag) {
                current_cb_velocity_acceleration =
                    cb_velocity_acceleration(
                        state_subject,
                        state_actor,
                        mu_actor
                    );

                current_acceleration += current_cb_velocity_acceleration;
            }

            if (de_sitter_flag) {
                current_de_sitter_acceleration =
                    de_sitter_acceleration (
                        state_subject,
                        state_actor,
                        state_primary,
                        mu_primary
                    );

                current_acceleration += current_de_sitter_acceleration;
            }

            if (lense_thirring_flag) {
                current_lense_thirring_acceleration =
                    lense_thirring_acceleration (
                        state_subject,
                        state_actor,
                        angular_momentum_actor
                    );

                current_acceleration += current_lense_thirring_acceleration;
            }

            if (wavi_flag) {
                current_wavi_acceleration =
                    wavi_acceleration (
                        state_subject,
                        state_actor,
                        mu_actor
                    );

                current_acceleration += current_wavi_acceleration;
            }

            return;
        }
    };

} // namespace jw

} // namespace tudat

#endif // JW_ACCELERATION_H
