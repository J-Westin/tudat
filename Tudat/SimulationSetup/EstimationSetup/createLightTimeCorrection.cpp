/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/jw_lighttime_correction.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"

#ifndef YEET_DEBUG
    #define YEET_DEBUG std::cout << __FILE__ << "\nLine: " << __LINE__ << std::endl
#endif // YEET_DEBUG

namespace tudat
{

namespace observation_models
{

//! Function to create object that computes a single (type of) correction to the light-time
std::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const std::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver )
{

    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;

    std::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    // Identify type of light time correction to be created.
    switch( correctionSettings->getCorrectionType( ) )
    {
    case first_order_relativistic:
    {
        // Check input consistency
        if( std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != nullptr )
        {
            // Retrieve list of bodies causing light time perturbation
            std::vector< std::string > perturbingBodies =
                    std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            // Retrieve mass and state functions for each perturbing body.
            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    throw std::runtime_error(
                                "Error when making 1st order relativistic light time correction, could not find body " +
                                perturbingBodies.at( i ) );
                }
                else
                {
                    // Set state function.
                    perturbingBodyStateFunctions.push_back(
                                std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                                                                         bodyMap.at( perturbingBodies[ i ] ), std::placeholders::_1 ) );

                    // Set gravitational parameter function.
                    perturbingBodyGravitationalParameterFunctions.push_back(
                                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                             bodyMap.at( perturbingBodies[ i ] )->
                                             getGravityFieldModel( ) ) );
                }
            }

            // Create light-time correction function
            lightTimeCorrection = std::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        transmitter.first, receiver.first,
                        std::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ) );

        }
        else
        {
            throw std::runtime_error(
                        "Error, correction settings type (1st order relativistic) does not coincide with data type." );
        }

        break;
    }

    case jw_lighttime: {

        std::shared_ptr<jw_lighttime_settings> jw_settings =
            std::dynamic_pointer_cast<jw_lighttime_settings>(correctionSettings);

        if (jw_settings != nullptr) {



            // Retrieve list of bodies causing light time perturbation
            std::vector< std::string > perturbingBodies =
                jw_settings->get_perturbing_bodies( );

            if (jw_settings->j2) {
                if (jw_settings->orientation_axes.size() != perturbingBodies.size()) {
                    throw std::runtime_error("Number or supplied orientation axes does not match number of bodies");
                }
                if (jw_settings->j2_coefficients.size() != perturbingBodies.size()) {
                    throw std::runtime_error("Number or supplied J2 coefficients does not match number of bodies");
                }
            }

            std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;
            std::vector< Eigen::Vector3d > unit_orientation_axes;

            // Retrieve mass and state functions for each perturbing body.
            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ ) {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 ) {
                    throw std::runtime_error(
                                "Error when making jw light time correction, could not find body " +
                                perturbingBodies.at( i ) );
                } else {
                    // Set state function.
                    perturbingBodyStateFunctions.push_back(
                        std::bind(
                            &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                            bodyMap.at( perturbingBodies[ i ] ),
                            std::placeholders::_1
                        )
                    );

                    // Set gravitational parameter function.
                    perturbingBodyGravitationalParameterFunctions.push_back(
                                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                             bodyMap.at( perturbingBodies[ i ] )->
                                             getGravityFieldModel( ) ) );

                    if (jw_settings->j2) {
                        unit_orientation_axes.push_back(jw_settings->orientation_axes[i].normalized());
                    }
                }



            }

            // Create light-time correction function
            lightTimeCorrection = std::make_shared< jw_lighttime_calculator >(
                perturbingBodyStateFunctions,
                perturbingBodyGravitationalParameterFunctions,
                jw_settings->j2_coefficients,
                unit_orientation_axes,
                perturbingBodies,
                transmitter.first,
                receiver.first,
                jw_settings->shapiro,
                jw_settings->second_order,
                jw_settings->velocity,
                jw_settings->j2,
                std::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ) );
        }

        else {
            throw std::runtime_error(
                "Error, correction settings type (jw_lighttime) does not coincide with data type." );
        }

        break;

    }

    default:
    {
        std::string errorMessage = "Error, light time correction type " +
                std::to_string( correctionSettings->getCorrectionType( ) ) + " not recognized.";
        throw std::runtime_error( errorMessage );

        break;
    }

    }
    return lightTimeCorrection;
}

}

}
