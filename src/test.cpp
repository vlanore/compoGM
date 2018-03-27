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
#include "distributions.hpp"
#include "interfaces.hpp"
#include "mcmc_moves.hpp"
#include "node_skeletons.hpp"

using namespace std;
using namespace tc;

struct Scale {
    static double move(double& value, double tuning = 1.0) {
        auto multiplier = tuning * (uniform(generator) - 0.5);
        value *= exp(multiplier);
        return multiplier;
    }
};

class Mean : public Component {
    double sum;
    int count;

  public:
    void add(double v) {
        sum += v;
        count++;
    }
    double mean() { return sum / count; }
};

int main() {
    Model m;
    m.component<OrphanNode<Exp>>("c1", 0.5, 1.0);
    m.component<UnaryNode<Exp>>("c2", 0.22).connect<Use<Value<double>>>("parent", "c1");

    m.component<SimpleMHMove<Scale, double>>("move1", 0.5)
        .connect<Use<Value<double>>>("target", "c1")
        .connect<Use<LogProb>>("logprob", "c1")
        .connect<Use<LogProb>>("logprob", "c2");

    m.component<Mean>("mean");

    m.driver("movescheduler", [](Go* move, Value<double>* c, Mean* mean) {

         for (int i = 0; i < 100000; i++) {
             move->go();
             mean->add(c->get_ref());
         }
         cout << 1 / mean->mean() << endl;

     }).connect("move1", "c1", "mean");

    Assembly a(m);
    a.call("movescheduler", "go");
}
