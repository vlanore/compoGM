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
#include "compoGM.hpp"
#include "partition.hpp"
#include "thread_helpers.hpp"

using namespace std;
using namespace compoGM;

struct M1 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
                         map<string, map<string, int>>& counts,
                         map<string, string>& condition_mapping) {
        m.component<Matrix<OrphanNormal>>("lambda", genes, conditions, 1, 3, pow(1.5, 2));

        m.component<Matrix<DeterministicUnaryNode<double>>>("exp", genes, conditions,
                                                            [](double x) { return pow(10, x); })
            .connect<MatrixToValueMatrix>("a", "lambda");

        m.component<Matrix<Poisson>>("K", genes, samples, 0)
            .connect<SetMatrix<int>>("x", counts)
            .connect<ManyToMany<ArraysMap<UseValue>>>("a", "exp", condition_mapping);
    }
};

void compute(int argc, char** argv) {
    if (argc < 2) {
        cerr << "usage:\n\tM1_bin <data_location>\n";
        exit(1);
    }

    Assembly assembly;
    {
        Model model;

        // Parsing data files
        string data_folder = argv[1];
        auto counts = parse_counts(data_folder + "/counts.tsv");
        auto samples = parse_samples(data_folder + "/samples.tsv");
        check_consistency(counts, samples);

        Partition all_genes(counts.genes, p.size);
        IndexSet pgenes = all_genes.my_partition();
        p.message("%d genes in partitioned gene list", pgenes.size());

        // graphical model
        p.message("Creating component model...");
        model.component<M1>("model", pgenes, samples.conditions, make_index_set(counts.samples),
                            counts.counts, samples.condition_mapping);

        // suffstats and metropolis hastings moves
        model.component<Matrix<PoissonSuffstat>>("poissonsuffstats", pgenes, samples.conditions)
            .connect<ManyToMany2D<UseValue>>("lambda", Address("model", "exp"))
            .connect<ManyToMany<ArraysRevMap<Use<Value<int>>>>>("values", Address("model", "K"),
                                                                samples.condition_mapping);

        model.component<Matrix<SimpleMHMove<Scale>>>("moves", pgenes, samples.conditions)
            .connect<ManyToMany2D<MoveToTarget<double>>>("target", Address("model", "lambda"))
            .connect<ManyToMany2D<DirectedLogProb>>("logprob", "poissonsuffstats",
                                                    LogProbSelector::Full);

        // assembly
        p.message("Instantiating assembly...");
        assembly.instantiate_from(model);
    }

    p.message("Preparations before running chain");
    auto all_lambdas = assembly.get_all<OrphanNormal>("model");
    auto all_moves = assembly.get_all<SimpleMHMove<Scale>>();
    auto all_suffstats = assembly.get_all<PoissonSuffstat>();

    // gathering suff stats
    for (auto&& suffstat : all_suffstats.pointers()) {
        suffstat->acquire();
    }

    // trace header
    auto trace = make_trace(all_lambdas, "tmp" + to_string(p.rank) + ".dat");
    trace.header();

    p.message("Running the chain");
    for (int iteration = 0; iteration < 5000; iteration++) {
        for (int rep = 0; rep < 10; rep++) {
            for (auto&& move : all_moves.pointers()) {
                move->move(1.0);
                move->move(0.1);
                move->move(0.01);
            }
        }
        trace.line();
    }
    for (auto&& move : all_moves.pointers()) {
        cerr << setprecision(3) << "Accept rate" << setw(40) << move->get_name() << "  -->  "
             << move->accept_rate() * 100 << "%" << endl;
    }
}

int main(int argc, char** argv) {
    auto threads = spawn(0, 2, compute, argc, argv);
    join(threads);
}
