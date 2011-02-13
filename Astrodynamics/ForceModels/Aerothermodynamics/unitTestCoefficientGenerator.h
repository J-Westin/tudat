/*!   \file hypersonicLocalInclinationAnalysis.cpp
 *    This file contains the unit test of the hypersonic local inclination
 *    analysis.
 *    Path              : /Astrodynamics/ForceModels\Aerothermodynamics/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : Bart Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : B.Romgens@student.tudelft.nl
 *
 *    Date created      : 4 February, 2011
 *    Last modified     : 8 February, 2011
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition,
 *        McGraw Hill, 2001.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *        edition, AIAA Education Series, 2006
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110402    D. Dirkx          First version of file.
 *      110802    D. Dirkx          Added Apollo check.
 */

#ifndef UNITTESTCOEFFICIENTGENERATOR_H
#define UNITTESTCOEFFICIENTGENERATOR_H

// Include statements.
#include "sphereSegment.h"
#include "capsule.h"
#include "hypersonicLocalInclinationAnalysis.h"
#include "writingOutputToFile.h"

namespace unit_tests
{
    bool testCoefficientGenerator( );
}

#endif // UNITTESTCOEFFICIENTGENERATOR_H

// End of file.
