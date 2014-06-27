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

    virtual void push_values(const Walker* walker)
    {
        (void) walker;
        std::cout << "This method should never be called." << std::endl;
        exit(1);
    }

    void initializeErrorEstimator(const int type, int n_c);

    void push_value(const double value);

    void update_mean(const double weight);

    void push_mean();

    double result();

    double error();

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

    std::string name;


    double _extract_mean();

    double _extract_mean_of_means();

};

}
