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

    Eigen::Vector3d schwarzschild_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    Eigen::Vector3d cb_velocity_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    class jw_acceleration : public basic_astrodynamics::AccelerationModel < Eigen::Vector3d > {
    private:
        std::function< Eigen::Vector6d() > state_function_subject;
        std::function< Eigen::Vector6d() > state_function_actor;
        std::function< double() > mu_function_actor;

        bool newton_flag, kinetic_flag, schwarzschild_flag,
             de_sitter_flag, cb_velocity_flag, cb_acceleration_flag;

        //std::shared_ptr<jw_acceleration_settings> acceleration_settings;

        Eigen::Vector3d current_acceleration;
        Eigen::Vector3d current_schwarzschild_acceleration;
        Eigen::Vector3d current_cb_velocity_acceleration;

    public:

        //using namespace tudat::simulation_setup;

        jw_acceleration(
            std::function< Eigen::Vector6d() > state_function_subject,
            std::function< Eigen::Vector6d() > state_function_actor,
            std::function< double() > mu_function_actor,
            bool newton_flag,
            bool kinetic_flag,
            bool schwarzschild_flag,
            bool de_sitter_flag,
            bool cb_velocity_flag,
            bool cb_acceleration_flag
            //std::shared_ptr<jw_acceleration_settings> acceleration_settings
        ):
            state_function_subject(state_function_subject),
            state_function_actor(state_function_actor),
            mu_function_actor(mu_function_actor),
            newton_flag(newton_flag),
            kinetic_flag(kinetic_flag),
            schwarzschild_flag(schwarzschild_flag),
            de_sitter_flag(de_sitter_flag),
            cb_velocity_flag(cb_velocity_flag),
            cb_acceleration_flag(cb_acceleration_flag)
        { }

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

        Eigen::Vector6d get_acceleration_components() {
            Eigen::Vector6d component_vector;
            component_vector <<
                0.0, 0.0, current_schwarzschild_acceleration.norm(),
                0.0, current_cb_velocity_acceleration.norm(), 0.0;

            return component_vector;
        }

        void updateMembers( const double currentTime = TUDAT_NAN ) {
            Eigen::Vector6d state_subject = state_function_subject();
            Eigen::Vector6d state_actor   = state_function_actor();
            double mu_actor = mu_function_actor();

            current_acceleration.setZero();
            current_schwarzschild_acceleration.setZero();
            current_cb_velocity_acceleration.setZero();

            if (schwarzschild_flag) {
                current_schwarzschild_acceleration =
                    schwarzschild_acceleration(
                        state_subject,
                        state_actor,
                        mu_actor
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

            return;
        }
    };

} // namespace jw

} // namespace tudat

#endif // JW_ACCELERATION_H
