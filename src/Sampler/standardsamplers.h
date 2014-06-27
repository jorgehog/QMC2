#pragma once

#include "Sampler.h"

#include "../Walker/Walker.h"

namespace QMC2
{

class SampleRadius : public Sampler
{
public:
    SampleRadius() :
        Sampler("radius")
    {

    }

    // Sampler interface
public:
    void push_value_from_walker(const Walker *walker)
    {
        double m = 0;

        for (const double & r : walker->abs_r)
        {
            m += r;
        }

        m  /= walker->_n_p();

        push_value(m);
    }
};

class CalculateAndSampleRadius : public Sampler
{
public:
    CalculateAndSampleRadius() :
        Sampler("radius")
    {

    }

    // Sampler interface
public:
    void push_value_from_walker(const Walker *walker)
    {
        double m = 0;
        for (const double & r2 : walker->r2)
        {
            m += sqrt(r2);
        }

        m /= walker->_n_p();

        push_value(m);
    }
};

class SampleRadiusSquared : public Sampler
{
public:
    SampleRadiusSquared() :
        Sampler("radiusSquared")
    {

    }
    // Sampler interface
public:
    void push_value_from_walker(const Walker *walker)
    {
        double m = 0;
        for (const double & r2 : walker->r2)
        {
            m += r2;
        }

        m /= walker->_n_p();

        push_value(m);
    }
};

class SampleRelativeDistance : public Sampler
{
public:
    SampleRelativeDistance() :
        Sampler("relativeDistance")
    {

    }

    // Sampler interface
public:
    void push_value_from_walker(const Walker *walker)
    {

        double m = 0;

        for (uint i = 0; i < walker->r_rel.n_cols; ++i)
        {
            for (uint j = i + 1; j < walker->r_rel.n_cols; ++j)
            {
                m += walker->r_rel(i, j);
            }
        }

        m /= walker->_n_p()*(walker->_n_p() - 1)/2;

        push_value(m);
    }
};

}
