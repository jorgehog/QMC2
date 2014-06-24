#include "Sampler.h"

#include "../ErrorEstimator/Blocking/Blocking.h"
#include "../ErrorEstimator/SimpleVar/SimpleVar.h"

using namespace QMC2;

Sampler::Sampler(std::string name) :
    m_counter(0),
    mm_counter(0),
    mean(0),
    mean_of_means(0),
    error(0)
{
    if (standardErrorEstimator_mean == BLOCKING)
    {
        setErrorEstimator(MEAN, new Blocking(m_pp, "blocking_out_" + name));
    }
    else
    {
        setErrorEstimator(MEAN, new SimpleVar(m_pp));
    }

    if (standardErrorEstimator_mean_of_means == BLOCKING)
    {
        setErrorEstimator(MEANOFMEANS, new Blocking(m_pp, "blocking_out_" + name));
    }
    else
    {
        setErrorEstimator(MEANOFMEANS, new SimpleVar(m_pp));
    }

}

void Sampler::update_mean(const double weight)
{

    if (queued_value == SKIPVALUE) return;

    double newVal = weight*queued_value;

    meanErrorEstimator->update_data(newVal);

    mean += newVal;
    m_counter++;

}

void Sampler::push_mean()
{
    if (m_counter == 0) return;

    double M = extract_mean();

    meanOfMeansErrorEstimator->update_data(M);


    mean_of_means += M;
    mm_counter++;

}

double Sampler::extract_mean_error()
{

    double err = meanErrorEstimator->estimate_error();

    meanErrorEstimator->reset();

    return err;

}

double Sampler::extract_mean_of_means_error()
{

    double err = meanOfMeansErrorEstimator->estimate_error();

    meanOfMeansErrorEstimator->reset();

    return err;

}

void Sampler::reset()
{

    meanErrorEstimator->reset();
    meanOfMeansErrorEstimator->reset();

    mean = 0;
    m_counter = 0;

    mean_of_means = 0;
    mm_counter = 0;


}

ParParams Sampler::m_pp;

int Sampler::standardErrorEstimator_mean = Sampler::SIMPLE;
int Sampler::standardErrorEstimator_mean_of_means = Sampler::SIMPLE;
