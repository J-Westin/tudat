/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JW_LIGHTTIME_CORRECTION_PARTIAL_H
#define TUDAT_JW_LIGHTTIME_CORRECTION_PARTIAL_H

#include <iostream>

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/jw_lighttime_correction.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
/*!
 * Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
 * \param singleBodyLightTimeCorrection Value of light-time correction for which partial is to be evaluated.
 * \param bodyGravitationalParameter Value of gravitational parameter
 * \return Partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
 */
double compute_jw_lighttime_partial_wrt_mu(
        const double singleBodyLightTimeCorrection,
        const double bodyGravitationalParameter
    );

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
/*!
 * Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
 * \param totalLightTimeCorrection Value of light-time correction for which partial is to be evaluated.
 * \param ppnParameterGamma Value of PPN parameter gamma.
 * \return Partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
 */
double compute_jw_lighttime_partial_wrt_gamma(
        const double totalLightTimeCorrection,
        const double ppnParameterGamma
    );

//! Class for computing the partial derivatives of 1st-order relativistic light-time correction.
class jw_lighttime_correction_partial : public LightTimeCorrectionPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param correctionCalculator Object used to compute light-time correction.
     */
    jw_lighttime_correction_partial( const std::shared_ptr< observation_models::jw_lighttime_calculator > correctionCalculator )
      : LightTimeCorrectionPartial( observation_models::jw_lighttime ),
        correctionCalculator_( correctionCalculator )
    {
        //std::cout << "creating jw correction partial calculator." << std::endl;

        perturbingBodies_ = correctionCalculator_->getPerturbingBodyNames( );
        perturbingBodyGravitationalParameterFunctions_ =
                correctionCalculator_->getPerturbingBodyGravitationalParameterFunctions( );
    }

    //! Destructor.
    ~jw_lighttime_correction_partial( ){ }

    //! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
    /*!
     * Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \param bodyIndex Index of body in correctionCalculator_->currentLighTimeCorrectionComponents_ for which partial is to
     * be computed. NOTE: this index must be smaller than the number of bodies perturbing the light-time.
     * \return Partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
     */
    SingleOneWayRangePartialReturnType wrtBodyGravitationalParameter(
            const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times,
            const int bodyIndex )
    {
        if( !( static_cast< int >( perturbingBodies_.size( ) ) > bodyIndex ) )
        {
            throw std::runtime_error( "Error, bodyIndex in jw_lighttime_correction_partial not consistent with contents" );
        }

        double partialValue = compute_jw_lighttime_partial_wrt_mu(
                    correctionCalculator_->getCurrentLightTimeCorrectionComponent( bodyIndex ),
                    //correctionCalculator_->getCurrentTotalShapiro( ),
                    perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( ) );
        return std::make_pair(
                    ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ),
                    ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    //! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
    /*!
     * Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
     * \param states States of link ends for observation for which partial is to be computed
     * \param times Times of link ends for observation for which partial is to be computed
     * \return Partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
     */
    SingleOneWayRangePartialReturnType wrtPpnParameterGamma(
            const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times )
    {
//        std::cout << correctionCalculator_->getCurrentTotalLightTimeCorrection( ) << std::endl;
//        std::cout << "Last tx time: " << correctionCalculator_->last_tx_time << std::endl;
//        std::cout << "Last rx time: " << correctionCalculator_->last_rx_time << std::endl;
//        std::cout << "Last ev time: " << correctionCalculator_->last_eval_time << std::endl;
//        std::cout << "Last tx state: " << correctionCalculator_->last_tx_state << std::endl;
//        std::cout << "Last rx state: " << correctionCalculator_->last_rx_state << std::endl;
//        std::cout << "Last pert state: " << correctionCalculator_->last_pert_state << std::endl;
//
//        std::cout << "Last velocity lighttime: " << correctionCalculator_->last_vu << std::endl;
//        std::cout << "Last shapiro  lighttime: " << correctionCalculator_->last_shp << std::endl;

        double partialValue = compute_jw_lighttime_partial_wrt_gamma(
                    correctionCalculator_->getCurrentTotalLightTimeCorrection( ),
                    //correctionCalculator_->getCurrentShapiro( ),
                    correctionCalculator_->getPpnParameterGammaFunction_( )( ) );
        return std::make_pair(
                    ( Eigen::Matrix< double, 1, Eigen::Dynamic >( 1, 1 ) << partialValue ).finished( ),
                    ( times[ 0 ] + times[ 1 ] ) / 2.0 );
    }

    //! Function to get the names of bodies causing light-time correction.
    /*!
     * Function to get the names of bodies causing light-time correction.
     * \return Names of bodies causing light-time correction.
     */
    std::vector< std::string > getPerturbingBodies( )
    {
        return perturbingBodies_;
    }


protected:

    //! Object used to compute light-time correction.
    std::shared_ptr< observation_models::jw_lighttime_calculator > correctionCalculator_;

    //! Names of bodies causing light-time correction.
    std::vector< std::string > perturbingBodies_;

    //! Set of functions returning the gravitational parameters of the gravitating bodies.
    std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;

};

}

}

#endif // JW_LIGHTTIME_CORRECTION_PARTIAL_H
