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
#include "utils.hpp"

/*
====================================================================================================
  ~*~ SimpleMHMove ~*~
  A generic Metropolis-Hastings move.
==================================================================================================*/
template <class M>
class SimpleMHMove : public Move, public tc::Component {
    using ValueType = typename M::ValueType;

    // config
    Value<ValueType>* target;
    Backup* target_backup;
    std::vector<LogProb*> log_probs;
    void add_log_prob(LogProb* ptr) { log_probs.push_back(ptr); }

    // internal stats
    int reject{0}, total{0};

  public:
    SimpleMHMove() {
        port("target", &SimpleMHMove::target);
        port("targetbackup", &SimpleMHMove::target_backup);
        port("logprob", &SimpleMHMove::add_log_prob);
    }

    void move(double tuning = 1.0) final {
        target_backup->backup();
        double log_prob_before =
            accumulate(log_probs.begin(), log_probs.end(), 0.0,
                       [](double acc, LogProb* ptr) { return acc + ptr->get_log_prob(); });
        double log_hastings = M::move(target->get_ref(), tuning);
        double log_prob_after =
            accumulate(log_probs.begin(), log_probs.end(), 0.0,
                       [](double acc, LogProb* ptr) { return acc + ptr->get_log_prob(); });
        bool accept = decide(exp(log_prob_after - log_prob_before + log_hastings));
        if (not accept) {
            target_backup->restore();
            reject++;
        }
        total++;
    }

    double accept_rate() const { return double(total - reject) / total; }
};
