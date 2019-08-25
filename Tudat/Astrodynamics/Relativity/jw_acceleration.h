#ifndef JW_ACCELERATION_H
#define JW_ACCELERATION_H

#include <functional>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat {

namespace jw {

    Eigen::Vector3d schwarzschild_acceleration(
        Eigen::Vector6d state_subject,
        Eigen::Vector6d state_actor,
        double mu_actor
    );

    class jw_acceleration : public basic_astrodynamics::AccelerationModel < Eigen::Vector3d > {
    private:
        std::function< Eigen::Vector6d() > state_function_subject;
        std::function< Eigen::Vector6d() > state_function_actor;
        std::function< double() > mu_function_actor;

        Eigen::Vector3d current_acceleration;

    public:
        jw_acceleration(
            std::function< Eigen::Vector6d() > state_function_subject,
            std::function< Eigen::Vector6d() > state_function_actor,
            std::function< double() > mu_function_actor
        ):
            state_function_subject(state_function_subject),
            state_function_actor(state_function_actor),
            mu_function_actor(mu_function_actor)
        { }

        ~jw_acceleration() { }

        Eigen::Vector3d getAcceleration () {
            return current_acceleration;
        }

        void updateMembers( const double currentTime = TUDAT_NAN ) {
            Eigen::Vector6d state_subject = state_function_subject();
            Eigen::Vector6d state_actor   = state_function_actor();
            double mu_actor = mu_function_actor();

            current_acceleration = schwarzschild_acceleration(
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
