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

    Sampler();

    virtual void push_values(Walker* walker)
    {
        (void) walker;
        std::cout << "This method should never be called." << std::endl;
        exit(1);
    }

    void queue_value(const double value)
    {
        queued_value = value;
    }

    void push_value(const double value)
    {
        queue_value(value);
        update_mean();
    }

    void update_mean(const double weight = 1);

    void push_mean();

    double extract_mean() {

#ifdef MPI_ON
        MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &m_counter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        
        double M = mean /= m_counter;

        mean = 0;
        m_counter = 0;


        return M;
    }

    double extract_mean_of_means() {

        double M = mean_of_means / mm_counter;

        mean_of_means = 0;
        mm_counter = 0;

        return M;

    }

    double extract_mean_error();

    double extract_mean_of_means_error();

    ErrorEstimator * meanErrorEstimator;
    ErrorEstimator * meanOfMeansErrorEstimator;

    enum meanType
    {
        MEAN,
        MEANOFMEANS
    };

    void setErrorEstimator(const meanType type, ErrorEstimator * errorEstimator)
    {
        switch (type) {
        case MEAN:
            meanErrorEstimator = errorEstimator;
            break;
        case MEANOFMEANS:
            meanOfMeansErrorEstimator = errorEstimator;
            break;
        default:
            std::cout << "Invalid type " << type << std::endl;
            exit(1);
            break;
        }
    }

    void reset();

    const double SKIPVALUE = 1337;


    static void setParParams(ParParams & pp)
    {
        Sampler::m_pp = &pp;
    }


protected:
    static ParParams *m_pp;

    unsigned long int m_counter;
    unsigned long int mm_counter;

    double queued_value;
    double mean;
    double mean_of_means;

    double error;

};

}
