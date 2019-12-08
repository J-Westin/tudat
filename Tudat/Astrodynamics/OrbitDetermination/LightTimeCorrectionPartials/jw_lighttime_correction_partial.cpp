/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/jw_lighttime_correction_partial.h"

#include <stdexcept>

namespace tudat
{

namespace observation_partials
{

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
double compute_jw_lighttime_partial_wrt_mu(
        const double singleBodyLightTimeCorrection, const double bodyGravitationalParameter )
{
    return singleBodyLightTimeCorrection / bodyGravitationalParameter;
}

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
double compute_jw_lighttime_partial_wrt_gamma(
        const double totalLightTimeCorrection, const double ppnParameterGamma )
{
    double thicc = totalLightTimeCorrection / ( ppnParameterGamma + 1.0 );
//    std::cout << "total light time correction = " << totalLightTimeCorrection << std::endl;
//    std::cout << "gamma = " << ppnParameterGamma << std::endl;
//    std::cout << "partial = " << thicc << std::endl;
//
//    throw std::runtime_error("AAAAAY STOP THIS IS SEEKENEENG");

    return thicc;
}



}

}
