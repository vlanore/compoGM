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
        m.connect<DirectedLogProb>(tc::PortAddress("logprob", move.address), target,
                                   LogProbSelector::X);
    }
};

// move that connects two components of possibly varying dimensions (component, array or matric for
// each) by using the correct connector
template <class Connector>
struct AdaptiveConnect : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress user, tc::Address provider) {
        auto is_matrix = [&m](tc::Address a) {
            if (m.is_composite(a)) {
                auto& contents = m.get_composite(a);
                if (contents.is_composite(contents.all_addresses().front())) {
                    return true;
                }
            }
            return false;
        };
        auto check = [](bool c) {
            if (!c) {
                std::cerr << "Error in AdaptiveConnect, connection malformed!\n";
                exit(1);
            }
        };
        auto dim1 = [&m](tc::Address a) { return m.get_composite(a).size(); };
        auto dim2 = [&m](tc::Address a) {
            auto& contents = m.get_composite(a);
            return contents.get_composite(contents.all_addresses().front()).size();
        };
        if (is_matrix(user)) {
            if (is_matrix(provider)) {  // matrix to matrix
                check(dim1(user) == dim1(provider) and dim2(user) == dim2(provider));
                m.connect<ManyToMany2D<Connector>>(user, provider);
            } else if (m.is_composite(provider)) {  // matrix to vector
                if (dim1(user) == dim1(provider)) {
                    m.connect<ManyToMany<ManyToOne<Connector>>>(user, provider);
                } else if (dim2(user) == dim1(provider)) {
                    m.connect<ManyToOne<ManyToMany<Connector>>>(user, provider);
                } else {
                    std::cerr << "Error in AdaptiveConnect: trying to connect matrix of size "
                              << dim1(user) << "x" << dim2(user) << " to vector of size "
                              << dim1(provider) << ". Cannot determine correct connection!\n";
                    exit(1);
                }
            } else {
                m.connect<ManyToOne<ManyToOne<Connector>>>(user, provider);
            }
        } else if (m.is_composite(user)) {
            if (is_matrix(provider)) {
                m.connect<ManyToMany2D<Connector>>(user, provider);
            } else if (m.is_composite(provider)) {
                m.connect<ManyToMany<Connector>>(user, provider);
            } else {
                m.connect<ManyToOne<Connector>>(user, provider);
            }
        } else {
            if (is_matrix(provider)) {
                m.connect<ManyToMany2D<Connector>>(user, provider);
            } else if (m.is_composite(provider)) {
                m.connect<ManyToMany<ManyToOne<Connector>>>(user, provider);
            } else {
                m.connect<ManyToOne<ManyToOne<Connector>>>(user, provider);
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
        while (model.is_composite()) {
            f();
        }
        f();
        std::string target_name_str = target_name_in_model.to_string();
        std::cout << "target is " << target_name_str << "\n";

        // FIXME considering a simple case where there are no deterministic nodes
        // in this case the markov blanket is simply nodes pointing to the target
        std::set<std::string> blanket;
        for (auto e : edges) {
            std::cout << e.first << ", " << e.second << "\n";
            if (tc::Address(e.second).first() == target_name_str) {
                blanket.insert(tc::Address(e.first).first());
            }
        }

        // connect move to things
        std::cout << "yolo\n";
        m.connect<MoveToTarget<ValueType>>(move, target);  // to target

        for (auto c : blanket) {
            std::cout << "Connect " << move.address.to_string() << " to "
                      << tc::Address(c).to_string() << "\n";
            m.connect<DirectedLogProb>(tc::PortAddress("logprob", move.address), c,
                                       LogProbSelector::Full);
        }
    }
};
