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

#include "compoGM.hpp"
#include "mpi_helpers.hpp"
#include "mpi_proxies.hpp"

using namespace std;
using namespace compoGM_thread;

struct M0 : public Composite {
    static void contents(Model& m, IndexSet& experiments, IndexSet& samples,
                         map<string, map<string, int>>& data) {
        m.component<OrphanExp>("alpha", 1, 1);
        m.component<OrphanExp>("mu", 1, 1);

        m.component<Array<Gamma>>("lambda", experiments, 1)
            .connect<ArrayToValue>("a", "alpha")
            .connect<ArrayToValue>("b", "mu");

        if (p.rank != 0) {  // slave only
            p.message("Got %d experiments", experiments.size());
            m.component<Matrix<Poisson>>("K", experiments, samples, 0)
                .connect<MatrixLinesToValueArray>("a", "lambda")
                .connect<SetMatrix<int>>("x", data);
        }
    }
};

void compute(int, char**) {
    IndexSet experiments_full{"e0", "e1"};
    IndexSet experiments = partition_slave(experiments_full, p);
    IndexSet samples{"s0", "s1"};
    map<string, map<string, int>> data{{"e0", {{"s0", 12}, {"s1", 13}}},
                                       {"e1", {{"s0", 17}, {"s1", 19}}}};

    Model m;
    m.component<M0>("model", experiments, samples, data);
    if (!p.rank) {
        // === master =============================================================================
        // mpi proxies
        for (int i = 1; i < p.size; i++) {
            m.component<ProbNodeProv>("alpha_proxy_" + std::to_string(i))
                .connect<UseValue>("target", Address("model", "alpha"));
            m.component<ProbNodeProv>("mu_proxy_" + std::to_string(i))
                .connect<UseValue>("target", Address("model", "mu"));
        }
        m.component<Array<ProbNodeUse>>("lambda_proxies", experiments)
            .connect<ArrayToValueArray>("target", Address("model", "lambda"));

        // moves
        m.component<SimpleMHMove<Scale>>("move_alpha")
            .connect<ConnectMove<double>>("target", "model", Address("model", "alpha"));

        m.component<SimpleMHMove<Scale>>("move_mu").connect<ConnectMove<double>>(
            "target", "model", Address("model", "mu"));

    } else {
        // === slaves =============================================================================
        // mpi proxies
        m.component<ProbNodeUse>("alpha_proxy")
            .connect<UseValue>("target", Address("model", "alpha"));
        m.component<ProbNodeUse>("mu_proxy").connect<UseValue>("target", Address("model", "mu"));

        m.component<Array<ProbNodeProv>>("lambda_proxies", experiments)
            .connect<ArrayToValueArray>("target", Address("model", "lambda"));

        // moves
        m.component<Array<SimpleMHMove<Scale>>>("move_lambda", experiments)
            .connect<ConnectMove<double>>("target", "model", Address("model", "lambda"));
    }

    // === mpi connections ========================================================================
    for (int i = 1; i < p.size; i++) {
        mpi_connect_master_slave(m, PortAddress("connection", "alpha_proxy_" + std::to_string(i)),
                                 PortAddress("connection", "alpha_proxy"), i);
        mpi_connect_master_slave(m, PortAddress("connection", "mu_proxy_" + std::to_string(i)),
                                 PortAddress("connection", "mu_proxy"), i);
    }
    for (auto experiment : experiments_full) {
        int dest_index = index_owner_slave(experiment, experiments_full, p);
        mpi_connect_master_slave(m, PortAddress("connection", "lambda_proxies", experiment),
                                 PortAddress("connection", "lambda_proxies", experiment),
                                 dest_index);
    }

    // std::stringstream ss;
    // m.print(ss);
    // p.message(ss.str());
    Assembly a(m);

    auto moves = a.get_all<Move>().pointers();
    p.message("Got %d moves", moves.size());
    auto proxies = a.get_all<Proxy>().pointers();
    p.message("Got %d proxies", proxies.size());

    auto trace = make_trace(a.get_all<Value<double>>("model"), "tmp.dat");
    trace.header();

    if (p.rank) {  // slaves broadcast their data
        for (auto proxy : proxies) {
            proxy->release();
        }
    }
    for (int iteration = 0; iteration < 50000; iteration++) {
        for (auto proxy : proxies) {
            proxy->acquire();
        }
        for (auto& move : moves) {
            move->move(1.0);
            move->move(0.1);
            move->move(0.01);
        }
        for (auto proxy : proxies) {
            proxy->release();
        }
        if (!p.rank) {
            trace.line();
        }
    }
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }
