#include "Sampler.h"

#include "../ErrorEstimator/Blocking/Blocking.h"
#include "../ErrorEstimator/SimpleVar/SimpleVar.h"

using namespace QMC2;

Sampler::Sampler(std::string name) :
    m_counter(0),
    mm_counter(0),
    mean(0),
    mean_of_means(0),
    name(name)
{



}

void Sampler::initializeErrorEstimator(const int type, int n_c)
{

    if (m_counter != 0 || mm_counter != 0)
    {
        throw std::runtime_error("Sampler not reset before initialized.");
    }

    sampleState = type;

    switch (standardErrorEstimator) {
    case BLOCKING:
        errorEstimator = new Blocking(n_c, m_pp, "blocking_out_" + name);

        break;
    case SIMPLE:
        errorEstimator = new SimpleVar(m_pp);

    default:
        errorEstimator = NULL;
    }
}

void Sampler::update_mean(const double weight)
{
    for (const double & sample : queued_samples)
    {

        double value = weight*sample;

        if (sampleState == MEAN)
        {
            errorEstimator->update_data(value);
        }

        mean += value;
        m_counter++;

    }

    queued_samples.clear();
}

void Sampler::push_value(const double value)
{
    queued_samples.push_back(value);
}


void Sampler::push_mean()
{
    if (m_counter == 0) return;

    double M = _extract_mean();

    if (sampleState == MEANOFMEANS)
    {
        errorEstimator->update_data(M);
    }

    mean_of_means += M;
    mm_counter++;

}

double Sampler::_extract_mean()
{

#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &m_counter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

    double M = mean /= m_counter;

    mean = 0;
    m_counter = 0;


    return M;
}

double Sampler::_extract_mean_of_means()
{

    double M = mean_of_means / mm_counter;

    mean_of_means = 0;
    mm_counter = 0;

    return M;

}

double Sampler::result()
{
    if (sampleState == MEAN)
    {
        return _extract_mean();
    }

    else if (sampleState == MEANOFMEANS)
    {
        return _extract_mean_of_means();
    }
}


double Sampler::error()
{

    double err = errorEstimator->estimate_error();

    errorEstimator->reset();

    return err;

}

void Sampler::reset()
{
    errorEstimator->reset();

    mean = 0;
    m_counter = 0;

    mean_of_means = 0;
    mm_counter = 0;

    queued_samples.clear();

}

ParParams Sampler::m_pp;

int Sampler::standardErrorEstimator = Sampler::SIMPLE;
