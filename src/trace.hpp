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

#include <fstream>
#include <iostream>
#include "interfaces.hpp"
#include "tinycompo.hpp"

template <class I>
class Trace {
    std::unique_ptr<std::ofstream> internal_stream{nullptr};
    std::ostream& os;
    tc::InstanceSet<I> components;

  public:
    Trace(tc::InstanceSet<I> components, std::ostream& os = std::cout)
        : os(os), components(components) {}
    Trace(tc::InstanceSet<I> components, std::string filename)
        : internal_stream(new std::ofstream(filename)),
          os(*internal_stream),
          components(components) {}

    void header() const {
        auto names = components.names();
        for (size_t i = 0; i < names.size() - 1; i++) {
            const auto& c = names[i];
            os << c.to_string() << '\t';
        }
        os << names.at(names.size() - 1).to_string() << std::endl;
    }

    void line() const {
        auto pointers = components.pointers();
        for (size_t i = 0; i < pointers.size() - 1; i++) {
            const auto& c = pointers[i];
            os << c->get_ref() << '\t';
        }
        os << pointers.at(pointers.size() - 1)->get_ref() << std::endl;
    }
};

template <class I, class Arg>
Trace<I> make_trace(tc::InstanceSet<I> components, Arg arg) {
    return Trace<I>(components, arg);
}
