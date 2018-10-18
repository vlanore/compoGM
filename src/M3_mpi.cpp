/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2018).
Contributors:
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

This software is a component-based library to write bayesian inference programs based on the
graphical m.

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

#include "compoGM_mpi.hpp"

using namespace std;
using namespace compoGM;

struct M3 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
        map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
        map<string, double>& size_factors) {
        // global variables
        m.component<OrphanNormal>("a0", 1, -2, 2);
        m.component<OrphanNormal>("a1", 1, 0, 2);
        m.component<OrphanExp>("sigma_alpha", 1, 1);

        // has to be on master's ghost to compute logprob of log10(alpha) which is affected by globs
        m.component<Matrix<OrphanNormal>>("log10(q)", genes, conditions, 1, 2, 2);  // changed

        // deterministic nodes
        m.connect<MapPower10>("log10(q)", "q");
        m.component<Array<Sum>>("q_bar", genes).connect<ArrayToValueMatrixLines>("parent", "q");
        m.component<Array<DeterministicTernaryNode<double>>>("log10(alpha_bar)", genes,
             [](double a0, double a1, double q_bar) { return log10(a0 + a1 / q_bar); })
            .connect<ArrayToValue>("a", "a0")
            .connect<ArrayToValue>("b", "a1")
            .connect<ArrayToValueArray>("c", "q_bar");

        // has to be on master's ghost because it is in the balnket of a0/a1/sigma_alpha
        m.component<Array<Normal>>("log10(alpha)", genes, 1)
            .connect<ArrayToValueArray>("a", "log10(alpha_bar)")
            .connect<ArrayToValue>("b", "sigma_alpha");

        if (p.rank) {  // slave-only variables
            m.connect<MapInversePower10>("log10(alpha)", "1/alpha");

            m.component<Matrix<GammaSR>>("tau", genes, samples, 1)
                .connect<MatrixLinesToValueArray>("a", "1/alpha")
                .connect<MatrixLinesToValueArray>("b", "1/alpha");

            m.component<Array<Constant<double>>>("sf", samples, 0)
                .connect<SetArray<double>>("x", size_factors);

            m.component<Matrix<DeterministicTernaryNode<double>>>(
                 "lambda", genes, samples, [](double a, double b, double c) { return a * b * c; })
                .connect<MatrixColumnsToValueArray>("a", "sf")
                .connect<ManyToMany<ArraysMap<UseValue>>>("b", "q", condition_mapping)
                .connect<MatrixToValueMatrix>("c", "tau");

            m.component<Matrix<Poisson>>("K", genes, samples, 0)
                .connect<SetMatrix<int>>("x", counts)
                .connect<MatrixToValueMatrix>("a", "lambda");
        }
    }
};

void compute(int argc, char** argv) {
    if (argc < 2) {
        cerr << "usage:\n\tM3_bin <data_location>\n";
        exit(1);
    }

    Model m;

    // Parsing data files
    string data_location = argv[1];
    auto counts = parse_counts(data_location + "/counts.tsv");
    auto samples = parse_samples(data_location + "/samples.tsv");
    auto size_factors = parse_size_factors(data_location + "/size_factors.tsv");
    check_consistency(counts, samples, size_factors);

    // partitioning genes for slaves
    Partition gene_partition(counts.genes, p.size - 1, 1);
    auto my_genes = gene_partition.my_partition();

    // graphical model
    m.component<M3>("model", my_genes, samples.conditions, make_index_set(counts.samples),
        counts.counts, samples.condition_mapping, size_factors.size_factors);

    // MPI components
    m.component<Bcast>("globals_handler")
        .connect<UseValue>("target", Address("model", "a0"))
        .connect<UseValue>("target", Address("model", "a1"))
        .connect<UseValue>("target", Address("model", "sigma_alpha"));

    m.component<Gather>("log10(q)_log10(alpha)_handler", gene_partition)
        .connect<OneToMany<UseValue>>("target", Address("model", "log10(alpha)"));

    // suffstats and metropolis hastings moves
    MpiMCMC mcmc(m, "model");
    mcmc.master_add("a0", scale);
    mcmc.master_add("a1", scale);
    mcmc.master_add("sigma_alpha", scale);
    mcmc.slave_add("log10(q)", shift);
    mcmc.slave_add("tau", scale);
    mcmc.slave_add("log10(alpha)", shift);
    mcmc.declare_moves();

    // instantiating assembly
    Assembly a(m);

    // gathering pointers and preparing trace
    auto moves = a.get_all<Move>().pointers();
    auto proxies = a.get_all<Proxy>().pointers();
    auto trace = make_trace(a.get_all<Value<double>>("model"), "tmp" + to_string(p.rank) + ".dat");
    if (!p.rank) { trace.header(); }

    // main loop
    if (p.rank) {  // slaves broadcast their data
        for (auto proxy : proxies) { proxy->release(); }
    }
    p.message("Go!");
    Chrono total_time;
    Chrono computing_time;
    Chrono writing_time;
    int nb_iterations = 500;
    for (int iteration = 0; iteration < nb_iterations; iteration++) {
        for (auto proxy : proxies) { proxy->acquire(); }
        computing_time.start();
        for (int i = 0; i < 10; i++) {
            for (auto move : moves) {
                move->move(1.0);
                move->move(0.1);
                move->move(0.01);
            }
        }
        computing_time.end();
        for (auto proxy : proxies) { proxy->release(); }
        if (!p.rank) {
            writing_time.start();
            trace.line();
            writing_time.end();
        }
    }
    double elapsed_time = total_time.end();
    compoGM::p.message("MCMC chain has finished in %fms (%fms/iteration)", elapsed_time,
        elapsed_time / nb_iterations);
    compoGM::p.message("Average computing time is %fms", computing_time.mean());
    if (!p.rank) compoGM::p.message("Average writing time is %fms", writing_time.mean());
}

int main(int argc, char** argv) {
    CE::master_and_ce1_only = true;
    mpi_run(argc, argv, compute);
}
