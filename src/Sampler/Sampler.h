#pragma once

#include "../defines.h"

#include "../structs.h"

#include <string>


namespace QMC2
{


class Walker;
class ErrorEstimator;

class Sampler {
public:

    Sampler(std::string name);

    virtual void push_values(Walker* walker)
    {
        (void) walker;
        std::cout << "This method should never be called." << std::endl;
        exit(1);
    }

    void update_mean(const double weight);

    void push_value(const double value);

    void push_mean();

    double extract_mean();

    double extract_mean_of_means();

    double extract_mean_error();

    double extract_mean_of_means_error();

    ErrorEstimator * errorEstimator;

    ErrorEstimator *getErrorEstimator() const
    {
        return errorEstimator;
    }

    enum _ErrorEstimators
    {
        SIMPLE,
        BLOCKING,
        NONE
    };

    enum _SampleState
    {
        MEAN,
        MEANOFMEANS
    };


    static int standardErrorEstimator;

    void reset();

    static void setParParams(ParParams & pp)
    {
        Sampler::m_pp = pp;
    }


protected:
    static ParParams m_pp;

    unsigned long int m_counter;
    unsigned long int mm_counter;

    double mean;
    double mean_of_means;

    std::vector<double> queued_samples;
    int sampleState;

};

}
