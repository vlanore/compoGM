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

using IndexSet = std::set<std::string>;

/*
====================================================================================================
  ~*~ Named array ~*~
==================================================================================================*/
template <class Element>
struct NamedArray : public tc::Composite {
    template <class... Args>
    static void contents(tc::Model& m, IndexSet indices, Args... args) {
        for (auto&& index : indices) {
            m.component<Element>(index, args...);
        }
    }
};

/*
====================================================================================================
  ~*~ Named Matrix ~*~
==================================================================================================*/
template <class Element>
using NamedMatrix = NamedArray<NamedArray<Element>>;

/*
====================================================================================================
  ~*~ Matrix / array connectors ~*~
==================================================================================================*/
template <class P2PConnector>
struct ArrayMetaOneToOne {
    static void connect(tc::Model& m, tc::PortAddress user, tc::Address provider) {
        auto user_elements = m.get_composite(user.address).all_component_names(0, true);
        auto provider_elements = m.get_composite(provider).all_component_names(0, true);
        bool lengths_match = user_elements.size() == provider_elements.size();
        if (!lengths_match) {
            throw tc::TinycompoException("ArrayMetaOneToOne: lengths of composites " +
                                         user.address.to_string() + " and " + provider.to_string() +
                                         "  don't match.");
        }
        for (int i = 0; i < static_cast<int>(user_elements.size()); ++i) {
            m.connect<P2PConnector>(tc::PortAddress(user.prop, user.address, user_elements.at(i)),
                                    tc::Address(provider, provider_elements.at(i)));
        }
    }
};
