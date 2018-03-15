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
    int i{17};
    void hello() { std::cout << i << " youpi tralala\n"; }
};

// candidate classes for tinycompo extension
struct Go {
    virtual void go() = 0;
};

struct _AbstractDriver {
    virtual void set_refs(std::vector<tc::Component*>) = 0;
};

template <class... Refs>
class Driver : public tc::Component, public _AbstractDriver {
    std::function<void(Refs...)> instructions;
    std::tuple<Refs...> refs;

    template <int... S>
    void call(seq<S...>) {
        instructions(std::get<S>(refs)...);
    }

    template <int i>
    void set_ref_helper(std::vector<tc::Component*>&) {}

    template <int i, class Head, class... Tail>
    void set_ref_helper(std::vector<tc::Component*>& ref_values) {
        std::get<i>(refs) = dynamic_cast<Head>(ref_values.at(i));
        set_ref_helper<i - 1, Tail...>(ref_values);
    }

    void set_refs(std::vector<tc::Component*> ref_values) {
        set_ref_helper<sizeof...(Refs) - 1, Refs...>(ref_values);
    }

    // void set_refs(Refs... ref_values) { set_ref_helper<sizeof...(Refs) - 1>(ref_values...); }

    void go() { call(typename gens<sizeof...(Refs)>::type()); }

  public:
    Driver(const std::function<void(Refs...)>& instructions) : instructions(instructions) {
        port("go", &Driver::go);
        port("refs", &Driver::set_refs);
    }
};

template <class... Addresses>
struct DriverConnect {
    static void helper(tc::Assembly&, std::vector<tc::Component*>&) {}

    template <class... Tail>
    static void helper(tc::Assembly& a, std::vector<tc::Component*>& result, tc::Address head,
                       Tail... tail) {
        result.push_back(&a.at(head));
        helper(a, result, tail...);
    }

    template <class... Tail>
    static void helper(tc::Assembly& a, std::vector<tc::Component*>& result, tc::PortAddress head,
                       Tail... tail) {
        auto provided_port = a.at(head.address).get<tc::Component>(head.prop);
        result.push_back(provided_port);
        helper(a, result, tail...);
    }

    static void _connect(tc::Assembly& a, tc::Address driver, Addresses... addresses) {
        std::vector<tc::Component*> result;
        helper(a, result, addresses...);
        a.at<_AbstractDriver>(driver).set_refs(result);
    }
};

int main() {
    tc::Model model;
    model.component<Driver<Hello*>>("test", [](Hello* r) { r->hello(); });
    model.component<Hello>("hello");
    model.connect<DriverConnect<tc::Address>>("test", tc::Address("hello"));

    tc::Assembly assembly(model);
    assembly.call("test", "go");
}
