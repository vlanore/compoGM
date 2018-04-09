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

#include <tinycompo.hpp>
#include "interfaces.hpp"

/*
====================================================================================================
  ~*~ UnaryNode ~*~
  An unary node is a node with exactly one parent in the graphical model. typically,  it's a node
  with a one-parameter distribution (such as Exp).
==================================================================================================*/
template <class PDS>
class UnaryNode : public Value<double>, public LogProb, public Backup, public tc::Component {
    double value{0};
    double bk_value{0};
    Value<double>* parent{nullptr};

  public:
    UnaryNode(double value) : value(value) { port("parent", &UnaryNode::parent); }
    double& get_ref() final { return value; }
    double get_log_prob() final { return PDS::full_log_prob(value, parent->get_ref()); }
    double get_log_prob_x() final { return PDS::partial_log_prob_x(value, parent->get_ref()); }
    double get_log_prob_a() final { return PDS::partial_log_prob_a(value, parent->get_ref()); }
    void backup() final { bk_value = value; }
    void restore() final { value = bk_value; }
};

/*
====================================================================================================
  ~*~ OrphanNode ~*~
  An orphan node is a node with no parent in the graphical model, typically because it has constant
  parameters.
==================================================================================================*/
template <class PDS>
class OrphanNode : public Value<double>, public LogProb, public Backup, public tc::Component {
    double value{0};
    double bk_value{0};
    std::function<double(double)> f;  // std::function is used as a way to store constructor args

  public:
    template <class... Args>
    OrphanNode(double value, Args... args)
        : value(value), f([args...](double v) { return PDS::partial_log_prob_x(v, args...); }) {}
    double& get_ref() final { return value; }
    double get_log_prob() final { return f(value); }
    void backup() final { bk_value = value; }
    void restore() final { value = bk_value; }
};
