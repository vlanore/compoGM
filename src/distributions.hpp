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
#include "utils.hpp"

/*
====================================================================================================
  ~*~ Exponential distribution ~*~
==================================================================================================*/
struct Exp {
    using ValueType = double;

    static double full_log_prob(double x, double lambda) { return log(lambda) - lambda * x; }

    static double partial_log_prob_x(double x, double lambda) { return -lambda * x; }

    // here a (first parameter) is lambda
    static double partial_log_prob_a(double x, double lambda) { return log(lambda) - lambda * x; }
};

/*
====================================================================================================
  ~*~ Gamma distribution ~*~
==================================================================================================*/
struct Gamma {
    using ValueType = double;

    static double full_log_prob(double x, double k, double theta) {
        return -log(std::tgamma(k)) - k * log(theta) + (k - 1) * log(x) - x / theta;
    }

    static double partial_log_prob_x(double x, double k, double theta) {
        return (k - 1) * log(x) - x / theta;
    }

    static double partial_log_prob_a(double, double k, double theta) {
        return -log(std::tgamma(k)) - k * log(theta) + (k - 1);
    }

    static double partial_log_prob_b(double x, double k, double theta) {
        return -k * log(theta) - x / theta;
    }
};

/*
====================================================================================================
  ~*~ Poisson distribution ~*~
==================================================================================================*/
struct Poisson {
    using ValueType = int;

    static double full_log_prob(int x, double lambda) {
        return x * log(lambda) - lambda - log_factorial(x);
    }

    static double partial_log_prob_x(int x, double lambda) {
        return x * log(lambda) - log_factorial(x);
    }

    static double partial_log_prob_a(int x, double lambda) { return x * log(lambda) - lambda; }
};

/*
====================================================================================================
  ~*~ Normal distribution ~*~
==================================================================================================*/
struct Normal {
    using ValueType = double;

    static double full_log_prob(double x, double mu, double sigma) {
        return -(x - mu) * (x - mu) / (2 * sigma * sigma) - 0.5 * log(2 * M_PI * sigma * sigma);
    }

    static double partial_log_prob_x(double x, double mu, double sigma) {
        return -(x - mu) * (x - mu) / (2 * sigma * sigma);
    }

    static double partial_log_prob_a(double x, double mu, double sigma) {
        return -(x - mu) * (x - mu) / (2 * sigma * sigma);
    }

    static double partial_log_prob_b(double x, double mu, double sigma) {
        return -(x - mu) * (x - mu) / (2 * sigma * sigma) - 0.5 * log(2 * M_PI * sigma * sigma);
    }
};
