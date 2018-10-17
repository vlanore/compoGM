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

template <class C>
struct AdaptiveOneToMany : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress user, tc::Address provider) {
        if (!m.is_composite(provider)) {
            m.connect<C>(user, provider);
        } else {
            auto addresses = m.get_composite(provider).all_addresses();
            for (auto address : addresses) { m.connect<C>(user, tc::Address(provider, address)); }
        }
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

template <typename ValueType>
struct ConnectIndividualMove : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress move, tc::Address model, tc::Address target,
        std::set<tc::Address> used_ss = {}) {  // use_ss is target->ss
        // getting digraph representation of graphical model
        auto& gmref = m.get_composite(model);
        auto digraph =
            gmref.get_digraph();  // digraph of all components inside graphical model composite
        auto vertices = digraph.first;
        auto edges = digraph.second;
        NodeName target_name_str = target.rebase(model).to_string();  // address of target in model

        // ss-related preparation
        auto is_supported = [used_ss](std::string n) {
            for (auto move : used_ss) {
                if (move.is_ancestor(tc::Address(n))) { return true; }
            }
            return false;
        };

        // Algorithm: blanket(targets, graph) =
        //   [prob nodes pointing to targets] U blanket([det nodes pointing to targets])
        std::function<NameSet(NameSet)> compute_blanket = [&compute_blanket, edges, gmref,
                                                              is_supported](NameSet targets) {
            if (targets.size() == 0) {
                return NameSet{};
            } else {
                NameSet partial_blanket, next_targets;
                for (auto e : edges) {
                    auto origin = edge_origin(e);
                    auto dest = edge_dest(e);
                    if (targets.count(dest) > 0) {  // points to a target
                        if (is_prob(origin, gmref)) {
                            if (is_supported(origin)) {
                                compoGM::p.message(
                                    "Skipping %s in blanket because it is supported by a suffstat",
                                    origin.c_str());
                            } else {
                                partial_blanket.insert(origin);
                            }
                        } else if (is_det(origin, gmref)) {
                            next_targets.insert(origin);
                        }
                    }
                }
                auto recursive_call = compute_blanket(next_targets);
                partial_blanket.insert(recursive_call.begin(), recursive_call.end());
                return partial_blanket;
            }
        };

        NameSet blanket = compute_blanket(NameSet{target_name_str});
        std::set<std::string> short_blanket;
        for (auto a : blanket) { short_blanket.insert(tc::Address(a).first()); }
        compoGM::p.message(
            "Blanket of %s is %s", target.c_str(), nameset_to_string(short_blanket).c_str());

        // connect move to things
        m.connect<MoveToTarget<ValueType>>(move, target);  // to target

        for (auto c : blanket) {
            m.connect<DirectedLogProb>(tc::PortAddress("logprob", move.address),
                tc::Address(model, tc::Address(c)), LogProbSelector::Full);
        }
    }
};

template <typename ValueType>
struct ConnectMove : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress move, tc::Address model, tc::Address target,
        std::set<tc::Address> used_ss = {}) {
        if (is_matrix(move.address, m)) {
            auto& tc = m.get_composite(target);
            auto element_addresses = tc.all_addresses();
            for (auto element_address : element_addresses) {
                m.connect<ConnectIndividualMove<ValueType>>(
                    tc::PortAddress(move.prop, tc::Address(move.address, element_address)), model,
                    tc::Address(target, element_address), used_ss);
            }

        } else if (is_array(move.address, m)) {
            auto target_adresses = m.get_composite(move.address).all_addresses();
            for (auto address : target_adresses) {
                m.connect<ConnectIndividualMove<ValueType>>(
                    tc::PortAddress(move.prop, tc::Address(move.address, address)), model,
                    tc::Address(target, address), used_ss);
            }
        } else {
            m.connect<ConnectIndividualMove<ValueType>>(move, model, target, used_ss);
        }
    }
};