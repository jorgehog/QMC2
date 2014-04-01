#pragma once


#include "../ErrorEstimator/SimpleVar/SimpleVar.h"

#include "../defines.h"

#include "../structs.h"

#include <string>


namespace QMC2
{


class Walker;

class Sampler {
public:

    Sampler() :
        m_counter(0),
        mm_counter(0),
        mean(0),
        mean_of_means(0),
        error(0)
    {
        setErrorEstimator(MEAN, new SimpleVar(pp));
        setErrorEstimator(MEANOFMEANS, new SimpleVar(pp));
    }

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

    void update_mean(const double weight = 1) {

        if (queued_value == SKIPVALUE) return;

        double newVal = weight*queued_value;

        meanErrorEstimator->update_data(newVal);

        mean += newVal;
        m_counter++;
    }

    void push_mean() {

        if (m_counter == 0) return;
        
        double M = extract_mean();

        meanOfMeansErrorEstimator->update_data(M);


        mean_of_means += M;
        mm_counter++;

    }

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

    double extract_mean_error() {

        double err = meanErrorEstimator->estimate_error();

        meanErrorEstimator->reset();

        return err;

    }

    double extract_mean_of_means_error() {

        double err = meanOfMeansErrorEstimator->estimate_error();

        meanOfMeansErrorEstimator->reset();

        return err;

    }

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

    void reset()
    {

        meanErrorEstimator->reset();
        meanOfMeansErrorEstimator->reset();

        mean = 0;
        m_counter = 0;

        mean_of_means = 0;
        mm_counter = 0;


    }

    const double SKIPVALUE = 1337;


    static void setParParams(const ParParams & pp)
    {
        Sampler::pp = pp;
    }


protected:
    static ParParams pp;

    unsigned long int m_counter;
    unsigned long int mm_counter;

    double queued_value;
    double mean;
    double mean_of_means;

    double error;

};

}
