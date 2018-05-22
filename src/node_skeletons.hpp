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
  ~*~ Constant ~*~
==================================================================================================*/
template <class ValueType>
class Constant : public Value<ValueType>, public tc::Component {
    ValueType x;  // just a buffer for computation of f(a, b, c)

  public:
    Constant(ValueType x) : x(x) { port("x", &Constant::x); }

    ValueType& get_ref() final { return x; }

    const ValueType& get_ref() const final { return x; }

    std::string debug() const final { return "Constant [" + std::to_string(x) + "]"; }
};

/*
====================================================================================================
  ~*~ DeterministicTernaryNode ~*~
==================================================================================================*/
template <class ValueType>
class DeterministicTernaryNode : public Value<ValueType>, public tc::Component {
    Value<double>*a, *b, *c;
    ValueType (*f)(double, double, double);
    mutable ValueType x;  // just a buffer for computation of f(a, b, c)

  public:
    DeterministicTernaryNode(ValueType (*f)(double, double, double)) : f(f) {
        port("a", &DeterministicTernaryNode::a);
        port("b", &DeterministicTernaryNode::b);
        port("c", &DeterministicTernaryNode::c);
    }

    ValueType& get_ref() final {
        x = f(a->get_ref(), b->get_ref(), c->get_ref());
        return x;
    }

    const ValueType& get_ref() const final {
        x = f(a->get_ref(), b->get_ref(), c->get_ref());
        return x;
    }

    std::string debug() const final {
        return "DeterministicTernaryNode [" + std::to_string(get_ref()) + "]";
    }
};

/*
====================================================================================================
  ~*~ DeterministicUnaryNode ~*~
==================================================================================================*/
template <class ValueType>
class DeterministicUnaryNode : public Value<ValueType>, public tc::Component {
    Value<double>* parent;
    ValueType (*f)(double);
    mutable ValueType x;  // just a buffer for computation of f(a, b, c)

  public:
    DeterministicUnaryNode(ValueType (*f)(double)) : f(f) {
        port("a", &DeterministicUnaryNode::parent);
    }

    ValueType& get_ref() final {
        x = f(parent->get_ref());
        return x;
    }

    const ValueType& get_ref() const final {
        x = f(parent->get_ref());
        return x;
    }

    std::string debug() const final {
        return "DeterministicUnaryNode [" + std::to_string(get_ref()) + "]";
    }
};

/*
====================================================================================================
  ~*~ BinaryNode ~*~
==================================================================================================*/
template <class PDS>
class BinaryNode : public Value<typename PDS::ValueType>,
                   public LogProb,
                   public Backup,
                   public tc::Component {
    using ValueType = typename PDS::ValueType;
    ValueType value{0};
    ValueType bk_value{0};
    Value<double>* a{nullptr};
    Value<double>* b{nullptr};

  public:
    BinaryNode(ValueType value) : value(value) {
        port("x", &BinaryNode::value);
        port("a", &BinaryNode::a);
        port("b", &BinaryNode::b);
    }
    ValueType& get_ref() final { return value; }
    const ValueType& get_ref() const final { return value; }
    double get_log_prob() final { return PDS::full_log_prob(value, a->get_ref(), b->get_ref()); }
    double get_log_prob_x() final {
        return PDS::partial_log_prob_x(value, a->get_ref(), b->get_ref());
    }
    double get_log_prob_a() final {
        return PDS::partial_log_prob_a(value, a->get_ref(), b->get_ref());
    }
    double get_log_prob_b() final {
        return PDS::partial_log_prob_b(value, a->get_ref(), b->get_ref());
    }
    void backup() final { bk_value = value; }
    void restore() final { value = bk_value; }
    std::string debug() const final { return "BinaryNode [" + std::to_string(value) + "]"; }
};

/*
====================================================================================================
  ~*~ UnaryNode ~*~
  An unary node is a node with exactly one parent in the graphical model. typically,  it's a node
  with a one-parameter distribution (such as Exp).
==================================================================================================*/
template <class PDS>
class UnaryNode : public Value<typename PDS::ValueType>,
                  public LogProb,
                  public Backup,
                  public tc::Component {
    using ValueType = typename PDS::ValueType;
    ValueType value{0};
    ValueType bk_value{0};
    Value<double>* parent{nullptr};  // FIXME, template parameter?

  public:
    UnaryNode(ValueType value) : value(value) {
        port("x", &UnaryNode::value);
        port("a", &UnaryNode::parent);
    }
    ValueType& get_ref() final { return value; }
    const ValueType& get_ref() const final { return value; }
    double get_log_prob() final { return PDS::full_log_prob(value, parent->get_ref()); }
    double get_log_prob_x() final { return PDS::partial_log_prob_x(value, parent->get_ref()); }
    double get_log_prob_a() final { return PDS::partial_log_prob_a(value, parent->get_ref()); }
    void backup() final { bk_value = value; }
    void restore() final { value = bk_value; }
    std::string debug() const override { return "UnaryNode [" + std::to_string(value) + "]"; }
};

/*
====================================================================================================
  ~*~ OrphanNode ~*~
  An orphan node is a node with no parent in the graphical model, typically because it has constant
  parameters.
==================================================================================================*/
template <class PDS>
class OrphanNode : public Value<typename PDS::ValueType>,
                   public LogProb,
                   public Backup,
                   public tc::Component {
    using ValueType = typename PDS::ValueType;
    ValueType value{0};
    ValueType bk_value{0};
    std::function<double(ValueType)> f;  // std::function is used as a way to store constructor args

  public:
    template <class... Args>
    OrphanNode(double value, Args... args)
        : value(value), f([args...](ValueType v) { return PDS::partial_log_prob_x(v, args...); }) {
        port("x", &OrphanNode::value);
    }
    ValueType& get_ref() final { return value; }
    const ValueType& get_ref() const final { return value; }
    double get_log_prob() final { return f(value); }
    void backup() final { bk_value = value; }
    void restore() final { value = bk_value; }
    std::string debug() const override { return "OrphanNode [" + std::to_string(value) + "]"; }
};
