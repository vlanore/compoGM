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

#include <iomanip>
#include <thread>
#include "compoGM.hpp"
// #include "mpi_helpers.hpp"
#include "partition.hpp"

using namespace std;
using namespace tc;
using namespace compoGM_thread;
using DUse = Use<Value<double>>;

struct M0 : public Composite {
    static void contents(Model& m, IndexSet& experiments, IndexSet& samples,
                         map<string, map<string, int>>& data) {
        m.component<OrphanNode<Exp>>("alpha", 1, 1);

        m.component<NArray<BinaryNode<Gamma>>>("beta", experiments, 1)
            .connect<NArrayMultiprovide<DUse>>("a", "alpha")
            .connect<NArrayMultiprovide<DUse>>("b", "alpha");

        m.component<NMatrix<UnaryNode<Poisson>>>("lambda", experiments, samples, 0)
            .connect<NArrays1To1<NArrayMultiprovide<DUse>>>("a", "beta")
            .connect<SetNMatrix<int>>("x", data);
    }
};

void compute(int, char**) {
    IndexSet experiments{"e0", "e1"};
    IndexSet samples{"s0", "s1"};
    map<string, map<string, int>> data{{"e0", {{"s0", 12}, {"s1", 13}}},
                                       {"e1", {{"s0", 17}, {"s1", 19}}}};

    Model model;
    model.component<M0>("model", experiments, samples, data);

    model.component<GammaShapeScaleSuffstat>("beta_suffstats")
        .connect<DUse>("a", Address("model", "alpha"))
        .connect<DUse>("b", Address("model", "alpha"))
        .connect<NArrayMultiuse<DUse>>("values", Address("model", "beta"));

    model.component<SimpleMHMove<Scale>>("move_alpha")
        .connect<MoveToTarget<double>>("target", Address("model", "alpha"))
        .connect<DirectedLogProb>("logprob", "beta_suffstats", LogProbSelector::Full);

    model.component<NArray<SimpleMHMove<Scale>>>("move_beta", experiments)
        .connect<NArrays1To1<MoveToTarget<double>>>("target", Address("model", "beta"))
        .connect<NArrays1To1<NArrayMultiuse<DirectedLogProb>>>(
            "logprob", Address("model", "lambda"), LogProbSelector::A);

    model
        .driver("p0_driver",
                [](Move* move, Proxy* suffstats) {
                    suffstats->acquire();
                    move->move(1.0);
                    move->move(0.1);
                    move->move(0.01);
                    suffstats->release();
                })
        .connect("move_alpha", "beta_suffstats");

    Assembly assembly(model);
    auto& p0_driver = assembly.at<_AbstractDriver>("p0_driver");
    auto move_beta = assembly.get_all<Move>("move_beta");

    auto trace = make_trace(assembly.get_all<Value<double>>("model"), "tmp.dat");
    trace.header();

    for (int iteration = 0; iteration < 5000; iteration++) {
        p0_driver.go();
        for (auto& move : move_beta.pointers()) {
            move->move(1.0);
            move->move(0.1);
            move->move(0.01);
        }
        trace.line();
    }
}

int main(int argc, char** argv) { compute(argc, argv); }
