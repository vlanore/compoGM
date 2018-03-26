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
#include <cmath>
#include "interfaces.hpp"

using namespace std;
using namespace tc;

template <class Functor>
class UnaryNode : public Value<double>, public LogProb, public Component {
    double value{0};
    Value<double>* parent{nullptr};

  public:
    UnaryNode(double value) : value(value) { port("parent", &UnaryNode::parent); }
    void set(double v) final { value = v; }
    double get() final { return value; }
    double get_log_prob() final { return Functor()(value, parent->get()); }
};

template <class Functor>
class OrphanNode : public Value<double>, public LogProb, public Component {
    double value{0};
    function<double(double)> f;

  public:
    template <class... Args>
    OrphanNode(double value, Args... args)
        : value(value), f([args...](double v) { return Functor()(v, args...); }) {}
    void set(double v) final { value = v; }
    double get() final { return value; }
    double get_log_prob() final { return f(value); }
};

struct Exp {
    double operator()(double x, double lambda) {
        return lambda * exp(- lambda * x);
    }
};

int main() {
    Model m;

    m.component<OrphanNode<Exp>>("c1", 0.0, 1.0);
    m.component<UnaryNode<Exp>>("c2", 0.0).connect<Use<Value<double>>>("parent", "c1");

    cout << Exp()(1.0, 2.0) << "Hello world!\n";
}
