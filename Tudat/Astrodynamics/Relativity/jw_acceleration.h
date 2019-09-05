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

    Eigen::Vector3d velocity_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    class jw_acceleration : public basic_astrodynamics::AccelerationModel < Eigen::Vector3d > {
    private:
        std::function< Eigen::Vector6d() > state_function_subject;
        std::function< Eigen::Vector6d() > state_function_actor;
        std::function< double() > mu_function_actor;

        bool newton_flag, schwarzschild_flag;

        //std::shared_ptr<jw_acceleration_settings> acceleration_settings;

        Eigen::Vector3d current_acceleration;

    public:

        //using namespace tudat::simulation_setup;

        jw_acceleration(
            std::function< Eigen::Vector6d() > state_function_subject,
            std::function< Eigen::Vector6d() > state_function_actor,
            std::function< double() > mu_function_actor,
            bool newton_flag,
            bool schwarzschild_flag
            //std::shared_ptr<jw_acceleration_settings> acceleration_settings
        ):
            state_function_subject(state_function_subject),
            state_function_actor(state_function_actor),
            mu_function_actor(mu_function_actor),
            newton_flag(newton_flag),
            schwarzschild_flag(schwarzschild_flag)
        { }

        ~jw_acceleration() { }

        Eigen::Vector3d getAcceleration () {
            return current_acceleration;
        }

        void updateMembers( const double currentTime = TUDAT_NAN ) {
            Eigen::Vector6d state_subject = state_function_subject();
            Eigen::Vector6d state_actor   = state_function_actor();
            double mu_actor = mu_function_actor();

            current_acceleration.setZero();

            current_acceleration +=
                schwarzschild_acceleration(
                    state_subject,
                    state_actor,
                    mu_actor
                );

            current_acceleration +=
                velocity_acceleration(
                    state_subject,
                    state_actor,
                    mu_actor
                );

            return;
        }
    };

} // namespace jw

} // namespace tudat

#endif // JW_ACCELERATION_H
