/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2018).
Contributors:
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

This software is a component-based library to write bayesian inference programs based on the
graphical model.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

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

#include <tinycompo.hpp>

// C++11 integer_sequence implementation :/
template <int...>
struct seq {};

template <int N, int... S>
struct gens : gens<N - 1, N - 1, S...> {};

template <int... S>
struct gens<0, S...> {
    typedef seq<S...> type;
};

// class for the test
struct Hello : public tc::Component {
    void hello() { std::cout << "youpi tralala\n"; }
};

// candidate classes for tinycompo extension
struct Go {
    virtual void go() = 0;
};

template <class... Refs>
class Driver : public tc::Component {
    std::function<void(Refs...)> instructions;
    std::tuple<Refs...> refs;

    template <int... S>
    void call(seq<S...>) {
        instructions(std::get<S>(refs) ...);
    }

    void go() { call(typename gens<sizeof...(Refs)>::type()); }

  public:
    Driver(const std::function<void(Refs...)>& instructions) : instructions(instructions) {
        port("go", &Driver::go);
    }
};

int main() {
    tc::Model model;
    model.component<Driver<Hello*>>("test", [](Hello* r) { r->hello(); });

    tc::Assembly assembly(model);
    assembly.call("test", "go");
}
