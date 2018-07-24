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
using namespace compoGM_thread;
using DUse = Use<Value<double>>;

struct M2 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
                         map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
                         map<string, double>& size_factors) {
        m.component<Matrix<OrphanNode<Normal>>>("log10(q)", genes, conditions, 1, 3, 1.5);
        m.component<Matrix<DeterministicUnaryNode<double>>>("q", genes, conditions,
                                                            [](double a) { return pow(10, a); })
            .connect<ManyToMany2D<DUse>>("a", "log10(q)");

        m.component<Array<OrphanNode<Normal>>>("log10(alpha)", genes, 1, -2, 2);
        m.component<Array<DeterministicUnaryNode<double>>>(
             "1/alpha", genes, [](double a) { return 1. / double(pow(10, a)); })
            .connect<ManyToMany<DUse>>("a", "log10(alpha)");

        m.component<Matrix<BinaryNode<GammaShapeRate>>>("tau", genes, samples, 1)
            .connect<ManyToMany<ManyToOne<DUse>>>("a", "1/alpha")
            .connect<ManyToMany<ManyToOne<DUse>>>("b", "1/alpha");

        m.component<Array<Constant<double>>>("sf", samples, 0)
            .connect<SetArray<double>>("x", size_factors);

        m.component<Matrix<DeterministicTernaryNode<double>>>(
             "lambda", genes, samples, [](double a, double b, double c) { return a * b * c; })
            .connect<ManyToOne<ManyToMany<DUse>>>("a", "sf")
            .connect<ManyToMany<ArraysMap<DUse>>>("b", "q", condition_mapping)
            .connect<ManyToMany2D<DUse>>("c", "tau");

        m.component<Matrix<UnaryNode<Poisson>>>("K", genes, samples, 0)
            .connect<SetMatrix<int>>("x", counts)
            .connect<ManyToMany2D<DUse>>("a", "lambda");
    }
};

void compute(int argc, char** argv) {
    if (argc < 2) {
        cerr << "usage:\n\tM2_bin <data_location>\n";
        exit(1);
    }
    Assembly assembly;
    {
        Model& model = assembly.get_model();

        // Parsing data files
        string data_location = argv[1];
        auto counts = parse_counts(data_location + "/counts.tsv");
        auto samples = parse_samples(data_location + "/samples.tsv");
        auto size_factors = parse_size_factors(data_location + "/size_factors.tsv");
        check_consistency(counts, samples, size_factors);
        auto pgenes = partition(counts.genes, p);
        p.message("%d genes in partitioned gene list.", pgenes.size());

        // graphical model
        p.message("Creating component model...");
        model.component<M2>("model", pgenes, samples.conditions, make_index_set(counts.samples),
                            counts.counts, samples.condition_mapping, size_factors.size_factors);

        // suffstats and metropolis hastings moves
        model.component<Array<GammaShapeRateSuffstat>>("tau_suffstats", pgenes)
            .connect<ManyToMany<OneToMany<DUse>>>("values", Address("model", "tau"))
            .connect<ManyToMany<DUse>>("a", Address("model", "1/alpha"))
            .connect<ManyToMany<DUse>>("b", Address("model", "1/alpha"));

        model.component<Matrix<SimpleMHMove<Shift>>>("move_q", pgenes, samples.conditions)
            .connect<ManyToMany2D<MoveToTarget<double>>>("target", Address("model", "log10(q)"))
            .connect<ManyToMany<ArraysRevMap<DirectedLogProb>>>(
                "logprob", Address("model", "K"), samples.condition_mapping, LogProbSelector::A);

        model.component<Matrix<SimpleMHMove<Scale>>>("move_tau", pgenes, samples.samples)
            .connect<ManyToMany2D<MoveToTarget<double>>>("target", Address("model", "tau"))
            .connect<ManyToMany2D<DirectedLogProb>>("logprob", Address("model", "K"),
                                                    LogProbSelector::A);

        model.component<Array<SimpleMHMove<Shift>>>("move_alpha", pgenes)
            .connect<ManyToMany<MoveToTarget<double>>>("target", Address("model", "log10(alpha)"))
            .connect<ManyToMany<DirectedLogProb>>("logprob", "tau_suffstats",
                                                  LogProbSelector::Full);
        // .connect<ManyToMany<OneToMany<DirectedLogProb>>>(
        //     "logprob", Address("model", "tau"), LogProbSelector::Full);

        // assembly
        p.message("Instantiating assembly...");
        assembly.instantiate();
    }

    p.message("Preparations before running chain");
    auto moves_q_tau = assembly.get_all<Move>(std::set<Address>{"move_q", "move_tau"});
    auto moves_alpha = assembly.get_all<Move>("move_alpha");
    auto suffstats_tau = assembly.get_all<Proxy>("tau_suffstats");
    auto all_watched = assembly.get_all<Value<double>>(
        std::set<Address>{Address("model", "log10(q)"), Address("model", "log10(alpha)")}, "model");
    assembly.get_model() = Model();

    // trace header
    auto trace = make_trace(all_watched, "tmp" + to_string(p.rank) + ".dat");
    trace.header();

    p.message("Running the chain");
    for (int iteration = 0; iteration < 5000; iteration++) {
        for (int big_rep = 0; big_rep < 10; big_rep++) {
            for (auto&& move : moves_q_tau.pointers()) {
                move->move(1.0);
                move->move(0.1);
                move->move(0.01);
            }
            for (auto&& ss : suffstats_tau.pointers()) {
                ss->acquire();
            }
            for (auto&& move : moves_alpha.pointers()) {
                for (int rep = 0; rep < 10; rep++) {
                    move->move(1.0);
                    move->move(0.1);
                    move->move(0.01);
                }
            }
            for (auto&& ss : suffstats_tau.pointers()) {
                ss->release();
            }
        }
        trace.line();
    }
}

int main(int argc, char** argv) {
    auto threads = spawn(0, 1, compute, argc, argv);
    join(threads);
    // mpi_run(argc, argv, compute);
}
