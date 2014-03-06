#pragma once

#include <QMC2.h>
#include "HF/src/libs/system/system.h"


class HFOrbitals : public QMC2::Orbitals
{
public:
    HFOrbitals(System * system);

    // Orbitals interface
public:
    void set_qnum_indie_terms(QMC2::Walker *walker, int i) {}
    double phi(const QMC2::Walker *walker, int particle, int q_num);
    double del_phi(const QMC2::Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const QMC2::Walker *walker, int particle, int q_num);

    System* m_system;
};
