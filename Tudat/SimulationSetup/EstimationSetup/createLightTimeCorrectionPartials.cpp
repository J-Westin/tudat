/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/jw_lighttime_correction.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/jw_lighttime_correction_partial.h"

#include <iostream>

namespace tudat
{

namespace observation_partials
{

//! Function to create a partial objects from list of light time corrections.
std::vector< std::shared_ptr< LightTimeCorrectionPartial > > createLightTimeCorrectionPartials(
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrectionList )
{
    std::vector< std::shared_ptr< LightTimeCorrectionPartial > > partialList;

    std::cout << __FILE__ << std::endl
        << "lightTimeCorrectionList contains " << lightTimeCorrectionList.size() << " elements" << std::endl
        << "Entry 0 has type " << lightTimeCorrectionList.at( 0 )->getLightTimeCorrectionType( ) << std::endl;
    // Iterate over all light time corrections
    for( unsigned int i = 0; i < lightTimeCorrectionList.size( ); i++ )
    {
        // Check type of light time correction
        switch( lightTimeCorrectionList.at( i )->getLightTimeCorrectionType( ) )
        {
        case observation_models::first_order_relativistic:
        {
            //std::cout << "Current entry has type " << lightTimeCorrectionList.at( i )->getLightTimeCorrectionType( ) << std::endl;

            std::shared_ptr< observation_models::FirstOrderLightTimeCorrectionCalculator > currentCorrection =
                    std::dynamic_pointer_cast< observation_models::FirstOrderLightTimeCorrectionCalculator >(
                        lightTimeCorrectionList.at( i ) );
            if( currentCorrection == nullptr )
            {
                throw std::runtime_error( "Error when making first order light time correction partial, type id observation_models::first_order_relativistic not consistent with class type." );
            }
            else
            {
                // Create partial of first-order relativistic light-time correction
                partialList.push_back(
                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionPartial >( currentCorrection ) );
            }

            break;
        }
        case observation_models::jw_lighttime:
        {
            std::shared_ptr< observation_models::jw_lighttime_calculator > currentCorrection =
                    std::dynamic_pointer_cast< observation_models::jw_lighttime_calculator >(
                        lightTimeCorrectionList.at( i ) );
            if( currentCorrection == nullptr )
            {
                throw std::runtime_error( "Error when making first order light time correction partial, type id observation_models::jw_lighttime not consistent with class type." );
            }
            else
            {
                // Create partial of first-order relativistic light-time correction
                partialList.push_back(
                            std::make_shared< jw_lighttime_correction_partial >( currentCorrection ) );
            }

            break;
        }
        default:
            break;
        }
    }

    return partialList;
}

}

}

