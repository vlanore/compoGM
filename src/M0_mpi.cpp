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

#include "M0_data.hpp"
#include "compoGM_mpi.hpp"

using namespace std;
using namespace compoGM;

struct M0 : public Composite {
    static void contents(
        Model& m, IndexSet& experiments, IndexSet& samples, map<string, map<string, int>>& data) {
        m.component<OrphanExp>("alpha", 1, 10);
        m.component<OrphanExp>("mu", 1, 1);

        m.component<Array<Gamma>>("lambda", experiments, 10)
            .connect<ArrayToValue>("a", "alpha")
            .connect<ArrayToValue>("b", "mu");

        if (p.rank != 0) {  // slave only
            m.component<Matrix<Poisson>>("K", experiments, samples, 0)
                .connect<MatrixLinesToValueArray>("a", "lambda")
                .connect<SetMatrix<int>>("x", data);
        }
    }
};

void compute(int, char**) {
    IndexSet experiments = gen_indexset("e", (p.size - 1) * 10);
    IndexSet samples = gen_indexset("s", 5000);
    auto data = gen_data(experiments, samples);
    Partition experiment_partition(experiments, p.size - 1, 1);
    auto my_experiments = experiment_partition.my_partition();

    Model m;
    m.component<M0>("model", my_experiments, samples, data);

    m.component<Bcast>("alpha_mu_handler")
        .connect<UseValue>("target", Address("model", "alpha"))
        .connect<UseValue>("target", Address("model", "mu"));

    m.component<Gather>("lambda_handler", experiment_partition)
        .connect<OneToMany<UseValue>>("target", Address("model", "lambda"));

    MpiMCMC mcmc(m, "model");
    mcmc.master_add("alpha", scale);
    mcmc.master_add("mu", scale);
    mcmc.slave_add("lambda", scale);
    mcmc.declare_moves();

    mcmc.go(100, 10);
}

int main(int argc, char** argv) {
    // CE::silence_all = true;
    CE::master_and_ce1_only = true;
    mpi_run(argc, argv, compute);
}
