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
#include "compoGM.hpp"

using namespace std;
using tc::Assembly;
using tc::Component;
using tc::Model;
using tc::Use;

class MyMean : public Component {
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
    m.component<OrphanExp>("k", 0.5, 1.0);
    m.component<OrphanExp>("theta", 0.5, 1.0);
    m.component<tc::Array<Gamma>>("array", 10, -1)
        .connect<tc::MultiProvide<Value<double>>>("a", "k")
        .connect<tc::MultiProvide<Value<double>>>("b", "theta")
        .connect<tc::ArraySet<double>>(
            "x", vector<double>{1.5, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0, 0.5, 0.5, 1.2});

    m.component<GammaShapeScaleSuffstat>("gammasuffstat")
        .connect<UseValue>("a", "k")
        .connect<UseValue>("b", "theta")
        .connect<tc::MultiUse<Value<double>>>("values", "array");

    m.component<SimpleMHMove<Scale>>("move1")
        .connect<UseValue>("target", "k")
        .connect<Use<Backup>>("targetbackup", "k")
        .connect<DirectedLogProb>("logprob", "k", LogProbSelector::A)
        .connect<DirectedLogProb>("logprob", "gammasuffstat", LogProbSelector::Full);

    m.component<SimpleMHMove<Scale>>("move2")
        .connect<UseValue>("target", "theta")
        .connect<Use<Backup>>("targetbackup", "theta")
        .connect<DirectedLogProb>("logprob", "theta", LogProbSelector::B)
        .connect<DirectedLogProb>("logprob", "gammasuffstat", LogProbSelector::Full);

    m.component<MyMean>("meank");
    m.component<MyMean>("meantheta");

    m.driver("movescheduler", [](Move* move, Move* move2, Value<double>* k, Value<double>* theta,
                                  MyMean* mean, MyMean* mean2) {
         for (int i = 0; i < 100000; i++) {
             move->move(0.5);
             move2->move(0.5);
             mean->add(k->get_ref());
             mean2->add(theta->get_ref());
         }
         cout << "Mean: " << mean->mean() * mean2->mean() << endl;
     }).connect("move1", "move2", "k", "theta", "meank", "meantheta");

    Assembly a(m);
    a.at<Proxy>("gammasuffstat").acquire();
    a.call("movescheduler", "go");
}
