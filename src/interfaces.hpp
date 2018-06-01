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
#include <cstdio>
#include <cstdlib>

/*
====================================================================================================
  ~*~ Go interface ~*~
  Used to start computation, e.g., to start a move scheduler.
==================================================================================================*/
struct Go {
    virtual void go() = 0;
};

/*
====================================================================================================
  ~*~ Go interface ~*~
  Used to start a move with tuning.
==================================================================================================*/
struct Move {
    virtual void move(double tuning = 1.0) = 0;
};

/*
====================================================================================================
  ~*~ LogProb interface ~*~
  Used to get a log prob from something (typically a node or suffstat). Partial versions can be
  provided but are not required.
==================================================================================================*/
struct LogProb {
    virtual double get_log_prob() = 0;                           // <- return the full log prob
    virtual double get_log_prob_x() { return get_log_prob(); };  // <- partial log probs
    virtual double get_log_prob_a() { return get_log_prob(); };  // (where a and b are generic names
    virtual double get_log_prob_b() { return get_log_prob(); };  // for parameters)
};

/*
====================================================================================================
  ~*~ Value interface ~*~
  Used to access a value by reference.
==================================================================================================*/
template <class ValueType>
struct Value {
    virtual ValueType& get_ref() = 0;
    virtual const ValueType& get_ref() const = 0;
};

/*
====================================================================================================
  ~*~ Backup interface ~*~
==================================================================================================*/
struct Backup {
    virtual void backup() = 0;
    virtual void restore() = 0;
};

/*
====================================================================================================
  ~*~ Proxy interface ~*~
  TODO used to update a proxy interface (such as a suffstat or distributed ghost). WIP.
==================================================================================================*/
struct Proxy {
    virtual void acquire() = 0;
    virtual void release() = 0;
};

/*
====================================================================================================
  ~*~ LogProbSelector ~*~
==================================================================================================*/
class LogProbSelector : public LogProb {
  public:
    enum Direction { X, A, B, Full, Invalid };

  private:
    Direction d;
    LogProb* ptr;

  public:
    LogProbSelector(Direction d = Invalid, LogProb* ptr = nullptr) : d(d), ptr(ptr) {}
    virtual ~LogProbSelector() = default;

    double get_log_prob() final {
        switch (d) {
            case X:
                return ptr->get_log_prob_x();
            case A:
                return ptr->get_log_prob_a();
            case B:
                return ptr->get_log_prob_b();
            case Full:
                return ptr->get_log_prob();
            default:
                printf("Error: trying to use uninitialized LogProbSelector!\n");
                exit(1);
        }
    }
};
