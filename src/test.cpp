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

#include <cmath>
#include <tinycompo.hpp>
#include "interfaces.hpp"

using namespace std;
using namespace tc;

template <class PDS>
class UnaryNode : public Value<double>, public LogProb, public Component {
    double value{0};
    Value<double>* parent{nullptr};

  public:
    UnaryNode(double value) : value(value) { port("parent", &UnaryNode::parent); }
    void set(double v) final { value = v; }
    double get() final { return value; }
    double get_log_prob() final { return PDS::full_log_prob(value, parent->get()); }
};

template <class PDS>
class OrphanNode : public Value<double>, public LogProb, public Component {
    double value{0};
    function<double(double)> f;

  public:
    template <class... Args>
    OrphanNode(double value, Args... args)
        : value(value), f([args...](double v) { return PDS::partial_x_log_prob(v, args...); }) {}
    void set(double v) final { value = v; }
    double get() final { return value; }
    double get_log_prob() final { return f(value); }
};

struct Exp {
    static double full_log_prob(double x, double lambda) { return log(lambda) - lambda * x; }

    static double partial_x_log_prob(double x, double lambda) { return -lambda * x; }

    static double partial_parent_log_prob(double x, double lambda) {
        return log(lambda) - lambda * x;
    }
};

template <class Move, class ValueType>
class SimpleMHMove : public Go, public Component {
    double tuning;
    Value<ValueType>* target;
    vector<LogProb*> log_probs;
    void add_log_prob(LogProb* ptr) { log_probs.push_back(ptr); }

  public:
    SimpleMHMove(double tuning = 1.0) : tuning(tuning) {
        port("target", &SimpleMHMove::target);
        port("logprob", &SimpleMHMove::add_log_prob);
    }

    void go() final {
        ValueType old_value = target->get();
        ValueType new_value = Move::move(tuning, old_value);
        double log_prob_before =
            accumulate(log_probs.begin(), log_probs.end(), 0.0,
                       [](double acc, LogProb* ptr) { return acc + ptr->get_log_prob(); });
        target->set(new_value);  // TODO TODO TODO change value so that it doesn't necessarily
                                 // creates a copy
                                 // TODO possibly new interface with backup inside it?
        double log_prob_after =
            accumulate(log_probs.begin(), log_probs.end(), 0.0,
                       [](double acc, LogProb* ptr) { return acc + ptr->get_log_prob(); });
        // TODO decide
    }
};

// TODO implement uniform draw in [0,1] (needed for move and scale)

struct Scale {
    static double move(double tuning, double value) {
        // TODO do the move :)
        return 1.0;
    }
};

int main() {
    Model m;
    m.component<OrphanNode<Exp>>("c1", 0.5, 1.0);
    m.component<UnaryNode<Exp>>("c2", 0.7).connect<Use<Value<double>>>("parent", "c1");

    m.component<SimpleMHMove<Scale, double>>("move1", 0.5)
        .connect<Use<Value<double>>>("target", "c1")
        .connect<Use<LogProb>>("logprob", "c1")
        .connect<Use<LogProb>>("logprob", "c2");

    Assembly a(m);
    cout << a.at<LogProb>("c2").get_log_prob() << endl;
}
