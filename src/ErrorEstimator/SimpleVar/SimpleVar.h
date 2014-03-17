#pragma once


#include "../ErrorEstimator.h"

#include "../../structs.h"

namespace QMC2
{

/*! \brief Calculates the simple variance of the sampled values.
 */
class SimpleVar : public ErrorEstimator {
public:
    SimpleVar(const ParParams &);
    
    double estimate_error();
    
    /*!
     * Overrides the default described in the superclass.
     * Does not store values in memory, but rather use sum variables.
     */
    void update_data(double val);

    void reset();

protected:
    
    double f; //!< sum variable used to calculate the mean
    double f2; //!< sum variable used to calulate the mean of squares.
    
    
};

}
