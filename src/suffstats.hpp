/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2018).
Contributors:
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

This software is a component-based library to write bayesian inference programs based on the
graphical model.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http:////www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

#pragma once

#include <cmath>
#include "interfaces.hpp"
#include "tinycompo.hpp"

/*
====================================================================================================
  ~*~ Gamma Suff Stat ~*~
==================================================================================================*/
class GammaSuffstat : public tc::Component, public LogProb, public Proxy {
    std::vector<Value<double>*> values;
    void add_value(Value<double>* p) { values.push_back(p); }

    Value<double>* k_;
    Value<double>* theta_;
    double sum{0};
    double sum_log{0};

  public:
    GammaSuffstat() {
        port("values", &GammaSuffstat::add_value);
        port("k", &GammaSuffstat::k_);
        port("theta", &GammaSuffstat::theta_);
    }

    void acquire() final {
        for (auto p : values) {
            auto value = p->get_ref();
            sum += value;
            sum_log += log(value);
        }
    }

    void release() final {
        sum = 0;
        sum_log = 0;
    }

    double get_log_prob() final {
        int N = values.size();
        double k = k_->get_ref();
        double theta = theta_->get_ref();
        return -N * log(std::tgamma(k)) - N * k * log(theta) + (k - 1) * sum_log -
               (1 / theta) * sum;
    }

    double get_log_prob_a() final {  // a = k
        int N = values.size();
        double k = k_->get_ref();
        double theta = theta_->get_ref();
        return -N * log(std::tgamma(k)) - N * k * log(theta) + (k - 1) * sum_log;
    }

    double get_log_prob_b() final {  // b = theta
        int N = values.size();
        double k = k_->get_ref();
        double theta = theta_->get_ref();
        return -N * k * log(theta) - (1 / theta) * sum;
    }
};

/*
====================================================================================================
  ~*~ Poisson Suff Stat ~*~
==================================================================================================*/
class PoissonSuffstat : public tc::Component, public LogProb, public Proxy {
    std::vector<Value<int>*> values;
    void add_value(Value<int>* p) { values.push_back(p); }

    Value<double>* lambda_;
    double sum{0};

  public:
    PoissonSuffstat() {
        port("values", &PoissonSuffstat::add_value);
        port("lambda", &PoissonSuffstat::lambda_);
    }

    void acquire() final {
        sum = 0;
        for (auto p : values) {
            auto value = p->get_ref();
            sum += value;
        }
    }

    void release() final { sum = 0; }

    double get_log_prob() final {  // a = lambda
        acquire();
        int N = values.size();
        double lambda = lambda_->get_ref();
        return -N * lambda + log(lambda) * sum;
    }
};
