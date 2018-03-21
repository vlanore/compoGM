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

#include <tinycompo.hpp>

using namespace std;
using tc::Model;
using tc::Assembly;
using tc::Component;
using tc::Composite;
using tc::Use;

//
// ==================================================================================================
template <class Type>
class Value : public Component {
    Type value;

  public:
    using type = Type;
    Value(Type value) : value(value) {}
    void set(Type new_value) { value = new_value; }
    Type get() { return value; }
};

//
// ==================================================================================================
template <class ValueType>
struct Node : public Composite {
    static void contents(Model& m, typename ValueType::type init) {
        m.component<ValueType>("value", init);
    }

    void after_construct() override { provide<ValueType>("port", "value"); }
};

//
// ==================================================================================================
template <class OperationType>
struct GraphOperation {};

//
// ==================================================================================================
template <class NodeType>
struct DiGraph : public Composite {
    template <class... Args>
    static void contents(Model& m, Args&&... args) {
        m.composite<NodeType>("node", std::forward<Args>(args)...);
    }

    void after_construct() override {
        port("child", &DiGraph::add_child);
    }

    vector<DiGraph*> children;
    void add_child(DiGraph* p) { children.emplace_back(p); }
};

//
// ==================================================================================================
int main() {
    using G = DiGraph<Node<Value<double>>>;

    Model model;
    model.composite<G>("a", 12).connect<Use<G>>("child", "b").connect<Use<G>>("child", "c");
    model.composite<G>("b", 13).connect<Use<G>>("child", "d");
    model.composite<G>("c", 14).connect<Use<G>>("child", "d");
    model.composite<G>("d", 15);

    model.dot_to_file();

    Assembly assembly(model);
}
