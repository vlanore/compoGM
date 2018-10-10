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

#include "arrays.hpp"
#include "interfaces.hpp"
#include "introspection.hpp"

struct DirectedLogProb {
    static void _connect(tc::Assembly& a, tc::PortAddress user, tc::Address provider,
        LogProbSelector::Direction direction) {
        auto& user_ref = a.at(user.address);
        auto& provider_ref = a.at<LogProb>(provider);
        user_ref.set(user.prop, LogProbSelector(direction, &provider_ref));
    }
};

template <typename ValueType>
struct MoveToTarget : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress move, tc::Address target) {
        m.connect<tc::Use<Value<ValueType>>>(move, target);
        m.connect<tc::Use<Backup>>(tc::PortAddress("targetbackup", move.address), target);
        m.connect<DirectedLogProb>(
            tc::PortAddress("logprob", move.address), target, LogProbSelector::X);
    }
};

// move that connects two components of possibly varying dimensions (component, array or matrix for
// each) by using the correct connector
template <class Connector, class... Args>
struct AdaptiveConnect : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress user, tc::Address provider, Args... args) {
        auto is_matrix = [&m](tc::Address a) {
            if (m.is_composite(a)) {
                auto& contents = m.get_composite(a);
                if (contents.is_composite(contents.all_addresses().front().first())) {
                    return true;
                }
            }
            return false;
        };
        if (is_matrix(user.address)) {
            if (is_matrix(provider)) {  // matrix to matrix
                m.connect<ManyToMany2D<Connector>>(user, provider, args...);
            } else if (m.is_composite(provider)) {  // matrix to vector
                m.connect<ManyToMany<ManyToOne<Connector>>>(user, provider, args...);
            } else {  // matrix to component
                m.connect<ManyToOne<ManyToOne<Connector>>>(user, provider, args...);
            }
        } else if (m.is_composite(user.address)) {
            if (is_matrix(provider)) {  // vector to matrix
                m.connect<ManyToMany<OneToMany<Connector>>>(user, provider, args...);
            } else if (m.is_composite(provider)) {  // vector to vector
                m.connect<ManyToMany<Connector>>(user, provider, args...);
            } else {  // vector to component
                m.connect<ManyToOne<Connector>>(user, provider, args...);
            }
        } else {
            if (is_matrix(provider)) {  // component to matrix
                m.connect<OneToMany<OneToMany<Connector>>>(user, provider, args...);
            } else if (m.is_composite(provider)) {  // component to vector
                m.connect<OneToMany<Connector>>(user, provider, args...);
            } else {  // component to component
                m.connect<Connector>(user, provider, args...);
            }
        }
    }
};

template <typename ValueType>
struct ConnectMove : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress move, tc::Address model, tc::Address target) {
        // getting digraph representation of graphical model
        auto& gmref = m.get_composite(model);
        auto digraph =
            gmref.get_digraph();  // digraph of all components inside graphical model composite
        auto vertices = digraph.first;
        auto edges = digraph.second;

        // getting address of target in model (TODO: should probably be refactored into
        // tinycompo)
        auto tmp = model;
        auto target_name_in_model = target;
        auto f = [&]() {
            if (target_name_in_model.first() == tmp.first()) {
                target_name_in_model = target_name_in_model.rest();
                tmp = tmp.rest();
            } else {
                std::cerr << "Error in ConnectMove: target is not in model!\n";
                exit(1);
            }
        };
        while (model.is_composite()) { f(); }
        f();
        NodeName target_name_str = target_name_in_model.to_string();

        // FIXME considering a simple case where there are no deterministic nodes
        // in this case the markov blanket is simply nodes pointing to the target
        NameSet blanket;

        // Algorithm: blanket(targets, graph) =
        //   [prob nodes pointing to targets] U blanket([det nodes pointing to targets])
        std::function<NameSet(NameSet)> compute_blanket = [edges](NameSet targets) {
            NameSet partial_blanket, next_targets;
            for (auto e: edges) {
                if (targets.count(edge_dest(e)) > 1) { // points to a target
                    if ()
                }
            }
        }

        for (auto e : edges) {
            if (tc::Address(e.second).first() == target_name_str) {
                blanket.insert(tc::Address(e.first).first());
            }
        }

        // connect move to things
        m.connect<AdaptiveConnect<MoveToTarget<ValueType>>>(move, target);  // to target

        for (auto c : blanket) {
            m.connect<AdaptiveConnect<DirectedLogProb, LogProbSelector::Direction>>(
                tc::PortAddress("logprob", move.address), tc::Address(model, c),
                LogProbSelector::Full);
        }
    }
};
