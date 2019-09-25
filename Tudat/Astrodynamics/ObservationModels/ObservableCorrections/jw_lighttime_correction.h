#ifndef JW_LIGHTTIME_CORRECTION_H
#define JW_LIGHTTIME_CORRECTION_H

#include <cmath>
#include <vector>

#include <functional>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"

namespace tudat {

namespace observation_models {

    class jw_lighttime_calculator: public LightTimeCorrection {
    private:
        //! Set of function returning the state of the gravitating bodies as a function of time.
        std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions_;

        //! Set of functions returning the gravitational parameters of the gravitating bodies.
        std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;

        //! Names of bodies causing light-time correction.
        std::vector< std::string > perturbingBodyNames_;

        bool shapiro_flag, second_order_flag,
             velocity_flag, j2_flag;


        //! Function returning the parametric post-Newtonian parameter gamma
        /*!
         *  Function returning the parametric post-Newtonian parameter gamma, a measure for the space-time curvature due to a
         *  unit rest mass (1.0 in GR)
         */
        std::function< double( ) > ppnParameterGammaFunction_;

        //! List of values (between 0 and 1) of how far into the light-time the state of each perturbing body is to be evaluated.
        std::vector< double > lightTimeEvaluationContribution_;

        //! List of light-time correction due to each separate perturbing body, as computed by last call to
        //! calculateLightTimeCorrection.
        std::vector< double > currentLighTimeCorrectionComponents_;

        //! Total light-time correction, as computed by last call to calculateLightTimeCorrection.
        double currentTotalLightTimeCorrection_;

    public:
        //! Constructor, takes and sets gravitating body properties.
        /*!
         *  Constructor, takes and sets gravitating body properties.
         *  \param perturbingBodyStateFunctions Set of function returning the state of the gravitating bodies as a function
         *  of time.
         *  \param perturbingBodyGravitationalParameterFunctions Set of functions returning the gravitational parameters of
         *  the gravitating bodies.
         *  \param perturbingBodyNames Names of bodies causing light-time correction.
         *  \param transmittingBody Name of transmitting body
         *  \param receivingBody Name of receiving body
         *  \param ppnParameterGammaFunction Function returning the parametric post-Newtonian parameter gamma, a measure
         *  for the space-time curvature due to a unit rest mass (default 1.0; value from GR)
         */
        jw_lighttime_calculator(
            const std::vector< std::function< Eigen::Vector6d( const double ) > >& perturbingBodyStateFunctions,
            const std::vector< std::function< double( ) > >& perturbingBodyGravitationalParameterFunctions,
            const std::vector< std::string > perturbingBodyNames,
            const std::string transmittingBody,
            const std::string receivingBody,
            const bool shapiro,
            const bool second_order,
            const bool velocity,
            const bool j2,
            const std::function< double( ) >& ppnParameterGammaFunction = [ ]( ){ return 1.0; } )

             :  LightTimeCorrection( first_order_relativistic ),
                perturbingBodyStateFunctions_( perturbingBodyStateFunctions ),
                perturbingBodyGravitationalParameterFunctions_( perturbingBodyGravitationalParameterFunctions ),
                perturbingBodyNames_( perturbingBodyNames ),
                shapiro_flag(shapiro),
                second_order_flag(second_order),
                velocity_flag(velocity),
                j2_flag(j2),
                ppnParameterGammaFunction_( ppnParameterGammaFunction ) {

                    currentTotalLightTimeCorrection_ = 0.0;
                    currentLighTimeCorrectionComponents_.resize( perturbingBodyNames_.size( ) );

                    // Check if perturbing body is transmitting/receiving body, and set evaluation time settings accordingly
                    for( unsigned int i = 0; i < perturbingBodyNames.size( ); i++ )
                    {
                        if( perturbingBodyNames.at( i ) == transmittingBody )
                        {
                            lightTimeEvaluationContribution_.push_back( 0.0 );
                        }
                        else if( perturbingBodyNames.at( i ) == receivingBody )
                        {
                            lightTimeEvaluationContribution_.push_back( 1.0 );
                        }
                        else
                        {
                            lightTimeEvaluationContribution_.push_back( 0.5 );
                        }
                    }
                }

        //! Destructor
        ~jw_lighttime_calculator( ){ }

        //! Function to calculate first-order relativistic light time correction due to set of gravitating point masses.
        /*!
         *  Function to calculate first order relativistic light time correction due to set of gravitating point masses,
         *  according to Eq. (11.17) of 2010 IERS conventions. Calculation are performed by calling ca
         *  calculateFirstOrderLightTimeCorrectionFromCentralBody function for each gravitating body.
         *  \param transmitterState State of transmitter at transmission time.
         *  \param receiverState State of receiver at reception time
         *  \param transmissionTime Time of signal transmission
         *  \param receptionTime Time of signal reception
         *  \return Total light time correction due to gravitating masses defined by perturbingBodyStateFunctions_ and
         *  perturbingBodyGravitationalParameterFunctions_ member variables.
         */
        double calculateLightTimeCorrection(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime
        );

        //! Function to compute the partial derivative of the light-time correction w.r.t. observation time
        /*!
         * Function to compute the partial derivative of the light-time correction w.r.t. observation time, equal to zero in this
         * model.
         * \param transmitterState State of transmitted at transmission time
         * \param receiverState State of receiver at reception time
         * \param transmissionTime Time of signal transmission
         * \param receptionTime Time of singal reception
         * \param fixedLinkEnd Reference link end for observation
         * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
         * \return Light-time correction w.r.t. observation time
         */
        double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType fixedLinkEnd,
            const LinkEndType linkEndAtWhichPartialIsEvaluated ) {

                throw std::runtime_error("jw_lighttime time derivative not implemented");

                return 0;

        }

        //! Function to compute the partial derivative of the light-time correction w.r.t. link end position
        /*!
         * Function to compute the partial derivative of the light-time correction w.r.t. link end position
         * \param transmitterState State of transmitted at transmission time
         * \param receiverState State of receiver at reception time
         * \param transmissionTime Time of signal transmission
         * \param receptionTime Time of singal reception
         * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
         * \return Light-time correction w.r.t. link end position
         */
        Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated
        ) {
            throw std::runtime_error("jw lighttime derivative wrt time has not been implemented");

            Eigen::Matrix< double, 3, 1 > ret;
            ret.setZero();
            return ret;
        }



        std::vector< std::string > getPerturbingBodyNames( ) {
            return perturbingBodyNames_;
        }

        std::vector< std::function< double( ) > > getPerturbingBodyGravitationalParameterFunctions( ) {
            return perturbingBodyGravitationalParameterFunctions_;
        }

        double getCurrentTotalLightTimeCorrection( ) {
            return currentTotalLightTimeCorrection_;
        }

        double getCurrentLightTimeCorrectionComponent( const int bodyIndex ) {
            return currentLighTimeCorrectionComponents_.at( bodyIndex );
        }

        std::function< double( ) > getPpnParameterGammaFunction_( ) {
            return ppnParameterGammaFunction_;
        }
    };

} // namespcae observation_models

} // namespace tudat

#endif // JW_LIGHTTIME_CORRECTION_H
