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

#include "tinycompo.hpp"

using NodeName = std::string;
using NameSet = std::set<NodeName>;
using Edge = std::pair<const NodeName, NodeName>;

NodeName edge_origin(Edge& edge) { return tc::Address(edge.first).to_string(); }

NodeName edge_dest(Edge& edge) { return tc::Address(edge.second).to_string(); }

NodeName first_part(NodeName name) { return tc::Address(name).first(); }

std::string nameset_to_string(NameSet s) {
    std::stringstream ss;
    ss << "{";
    for (auto e : s) { ss << e << " "; }
    ss << "}";
    return ss.str();
}

bool is_matrix(tc::Address a, const tc::Model& m) {
    if (m.is_composite(a)) {
        auto& contents = m.get_composite(a);
        if (contents.size() == 0) {
            compoGM::p.message("WARNING: %s is an empty composite!", a.to_string().c_str());
            return false;
        }
        auto element = contents.all_addresses().front().first();
        if (contents.is_composite(element)) { return true; }
    }
    return false;
}

bool is_array(tc::Address a, const tc::Model& m) { return !is_matrix(a, m) and m.is_composite(a); }

IndexSet get_array_indices(tc::Address a, const tc::Model& m) {
    auto names = m.get_composite(a).all_component_names();
    return make_index_set(names);
}

template <class T>
bool has_type(tc::Address a, const tc::Model& m) {
    if (m.is_composite(a)) {
        auto contents = m.get_composite(a);
        if (contents.size() == 0) {
            return false;
        } else {
            auto element = contents.all_addresses().front();
            if (is_matrix(a, m)) {
                auto subcontents = contents.get_composite(element.first());
                auto subelement = subcontents.all_addresses().front();
                return subcontents.has_type<T>(subelement);
            } else {
                return contents.has_type<T>(element);
            }
        }
    } else {
        return m.has_type<T>(a);
    }
}

bool is_prob(tc::Address a, const tc::Model& m) { return has_type<LogProb>(a, m); }

bool is_det(tc::Address a, const tc::Model& m) {
    return not has_type<LogProb>(a, m) and
           (has_type<Value<double>>(a, m) or has_type<Value<int>>(a, m));
}