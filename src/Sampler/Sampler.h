#pragma once


#include "../defines.h"

#include <string>


namespace QMC2
{


class Walker;

class Sampler {
public:

    Sampler() {
        m_counter = mm_counter = 0;
        mean = mean_of_means = 0;
    }

    void setValue(const Walker* walker) {
        queue_value(calculateValue(walker));
    }

    void queue_value(const double value) {
        queued_value = value;
    }

    void update_mean(const double weight = 1) {

        if (queued_value == SKIPVALUE) return;

        double newVal = weight*queued_value;

        mean += newVal;
        m_counter++;
    }

    void push_mean() {

        if (m_counter == 0) return;
        
        
        mean_of_means += extract_mean();
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

    const double SKIPVALUE = 1337;


protected:

    unsigned long int m_counter;
    unsigned long int mm_counter;

    double queued_value;
    double mean;
    double mean_of_means;


    virtual double calculateValue(const Walker * walker) {

        (void) walker;

        return 0.0;
    }

};

}
