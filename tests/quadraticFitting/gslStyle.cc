#include <gsl/gsl_multifit.h>
#include <vector>

/* jmyers sep 2010

   copy-pasted from the LSST linkTracklets.cc.  Compare with Kubica's code.

 */ 



void getBestFitVelocityAndAcceleration(std::vector<double> positions, const std::vector<double> & times,
                                       double & velocity, double &acceleration, double &position0)
{
    if (positions.size() < 1) {
        velocity = 0; 
        acceleration = 0;
        position0 = 0;
        return;
    }

    // try to get these guys along the same 180-degree stretch of degrees...
    double p0 = positions.at(0);
    for (unsigned int i = 1; i < positions.size(); i++) {
        while ( positions.at(i) - p0 > 180) {
            positions.at(i) -= 360;
        }
        while ( p0 - positions.at(i) > 180) {
            positions.at(i) += 360;
        }
    }
    
    // just ignore this for now... we're only benchmarking.
    // if (positions.size() != times.size()) {
    //     throw LSST_EXCEPT(ProgrammerErrorException,
    //                       "getBestFitVelocityAndAcceleration: position and time vectors not same size!");
    // }

    /* we're using GSL for this. this is roughly adapted from the GSL
     * documentation; see
     * http://www.gnu.org/software/gsl/manual/html_node/Fitting-Examples.html
     */ 

    gsl_vector * y = gsl_vector_alloc(positions.size());
    gsl_matrix * X = gsl_matrix_alloc(positions.size(), 3);
        
    for (unsigned int i = 0; i < positions.size(); i++) {
        gsl_matrix_set(X, i, 0, 1.0);
        gsl_matrix_set(X, i, 1, times.at(i) );
        gsl_matrix_set(X, i, 2, times.at(i)*times.at(i) );
        gsl_vector_set(y, i, positions.at(i));
    }
    gsl_vector * c = gsl_vector_alloc(3); // times*c = positions 
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (positions.size(), 3);
    gsl_matrix * covariance = gsl_matrix_alloc(3,3);
    double chiSquared = 0;    
    // TBD: check return values for error, etc
    gsl_multifit_linear(X, y, c, covariance, &chiSquared, work);
    position0    = gsl_vector_get(c,0);
    velocity     = gsl_vector_get(c,1);
    acceleration = 2 * gsl_vector_get(c,2); // we need to multiply by 2 due to the 1/2 part of 1/2 at^2 

    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_matrix_free(covariance);
    gsl_multifit_linear_free(work);
}


