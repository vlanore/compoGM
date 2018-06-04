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
struct GammaShapeScaleSSFormula {
    static double full_log_prob(double k, double theta, int N, double sum, double sum_log) {
        return -N * log(std::tgamma(k)) - N * k * log(theta) + (k - 1) * sum_log -
               (1 / theta) * sum;
    }

    static double partial_log_prob_a(double k, double, int N, double, double sum_log) {  // a = k
        return -N * log(std::tgamma(k)) + (k - 1) * sum_log;
    }

    static double partial_log_prob_b(double k, double theta, int N, double sum,
                                     double) {  // b = theta
        return -N * k * log(theta) - (1 / theta) * sum;
    }
};

struct GammaShapeRateSSFormula {
    static double full_log_prob(double alpha, double beta, int N, double sum, double sum_log) {
        return N * alpha * log(beta) - N * log(std::tgamma(alpha)) + (alpha - 1) * sum_log -
               beta * sum;
    }

    static double partial_log_prob_a(double alpha, double beta, int N, double,
                                     double sum_log) {  // a = alpha
        return N * alpha * log(beta) - N * log(std::tgamma(alpha)) + (alpha - 1) * sum_log;
    }

    static double partial_log_prob_b(double alpha, double beta, int N, double sum,
                                     double) {  // b = beta
        return N * alpha * log(beta) - beta * sum;
    }
};

template <class Formula>
class GammaSSTemplate : public tc::Component, public LogProb, public Proxy {
    std::vector<Value<double>*> values;
    void add_value(Value<double>* p) { values.push_back(p); }

    Value<double>* a_;
    Value<double>* b_;
    double sum{0};
    double sum_log{0};

  public:
    GammaSSTemplate() {
        port("values", &GammaSSTemplate::add_value);
        port("a", &GammaSSTemplate::a_);
        port("b", &GammaSSTemplate::b_);
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
        return Formula::full_log_prob(a_->get_ref(), b_->get_ref(), values.size(), sum, sum_log);
    }

    double get_log_prob_a() final {  // a = k
        return Formula::partial_log_prob_a(a_->get_ref(), b_->get_ref(), values.size(), sum,
                                           sum_log);
    }

    double get_log_prob_b() final {  // b = theta
        return Formula::partial_log_prob_b(a_->get_ref(), b_->get_ref(), values.size(), sum,
                                           sum_log);
    }
};

using GammaShapeRateSuffstat = GammaSSTemplate<GammaShapeRateSSFormula>;
using GammaShapeScaleSuffstat = GammaSSTemplate<GammaShapeScaleSSFormula>;

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

    void release() final {}

    double get_log_prob() final {  // a = lambda
        int N = values.size();
        double lambda = lambda_->get_ref();
        return -N * lambda + log(lambda) * sum;
    }
};
