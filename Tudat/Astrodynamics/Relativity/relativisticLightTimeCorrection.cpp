/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Relativity/relativisticLightTimeCorrection.h"

namespace tudat
{

namespace relativity
{

    //! Function to calculate first order relativistic light time correction due to a gravitating point mass.
    double calculateFirstOrderLightTimeCorrectionFromCentralBody( const double bodyGravitationalParameter,
                                                                  const Eigen::Vector3d& transmitterPosition,
                                                                  const Eigen::Vector3d& receiverPosition,
                                                                  const Eigen::Vector3d& centralBodyPosition,
                                                                  const double ppnParameterGamma )
    {
        // Calculate Euclidean geometric distances between transmitter, receiver and gravitating body.
        double distanceToReceiver = ( receiverPosition - centralBodyPosition ).norm( );
        double distanceToTransmitter = ( transmitterPosition - centralBodyPosition ).norm( );
        double linkEuclideanDistance = ( transmitterPosition - receiverPosition ).norm( );

        // Calculate and return light time correction.
        return ( 1.0 + ppnParameterGamma ) * bodyGravitationalParameter * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * std::log(
                    ( distanceToReceiver + distanceToTransmitter + linkEuclideanDistance ) /
                    ( distanceToReceiver + distanceToTransmitter - linkEuclideanDistance ) );

    }

    //! Function to calculate gradient of first order relativistic light time correction due to a gravitating point mass.
    Eigen::Matrix< double, 1, 3 > calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
            const double bodyGravitationalParameter,
            const Eigen::Vector3d& transmitterPosition,
            const Eigen::Vector3d& receiverPosition,
            const Eigen::Vector3d& centralBodyPosition,
            const bool evaluateGradientAtReceiver,
            const double ppnParameterGamma )
    {
        Eigen::Vector3d relativePositionVector = ( receiverPosition - transmitterPosition );
        double receiverDistance = ( receiverPosition - centralBodyPosition ).norm( );
        double transmitterDistance = ( transmitterPosition - centralBodyPosition ).norm( );
        double linkEndDistance = relativePositionVector.norm( );

        Eigen::Matrix< double, 1, 3 > gradient = ( receiverDistance + transmitterDistance ) *
                ( relativePositionVector.normalized( ) ).transpose( );
        if( evaluateGradientAtReceiver )
        {
           gradient -= relativePositionVector.norm( ) * ( receiverPosition.normalized( ) ).transpose( );
        }
        else
        {
            gradient += relativePositionVector.norm( ) * ( transmitterPosition.normalized( ) ).transpose( );\
        }

        return 2.0 * ppnParameterGamma * bodyGravitationalParameter * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * gradient /
                ( ( receiverDistance + transmitterDistance ) * ( receiverDistance + transmitterDistance ) -
                      linkEndDistance * linkEndDistance );


    }

    double calculate_velocity_lighttime(
        double mu,
        Eigen::Vector3d x_a, // Transmitter position
        Eigen::Vector3d x_b, // Receiver position
        Eigen::Vector6d s_p, // Central body state
        double gamma_ppn
    ) {
        double c = physical_constants::SPEED_OF_LIGHT;

        Eigen::Vector3d R_ab = x_b - x_a;
        double r_ab = R_ab.norm();

        Eigen::Vector3d N_ab = (1.0/r_ab)*R_ab;

        Eigen::Vector3d x_p, v_p, beta;
        x_p = s_p.segment(0, 3);
        v_p = s_p.segment(3, 3);

        beta = (1.0/c)*v_p;
        double b = beta.norm();

        double g = 1.0/std::sqrt( 1.0 - b*b );

        Eigen::Vector3d R_pa, R_pb;
        R_pa = x_a + (g*g/(1.0+g)) * (beta.dot(x_a - x_p)) * beta - x_p;
        R_pb = x_b + (g*g/(1.0+g)) * (beta.dot(x_b - x_p)) * beta - x_p;

        double r_pb = R_pb.norm();

        double factor, frac_t1, frac_t2;

        factor = 2*mu*g*(1.0 - N_ab.dot(beta))/(c*c*c);
        frac_t1 = (R_pa + g*r_ab*beta).norm() + r_pb;
        frac_t2 = g*r_ab*(1.0 - beta.dot(N_ab));

        return factor * std::log( (frac_t1 + frac_t2) / (frac_t1 - frac_t2) );
    }

}

}

