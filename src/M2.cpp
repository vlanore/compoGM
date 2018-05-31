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
#include "mpi_helpers.hpp"
#include "partition.hpp"

using namespace std;
using namespace tc;
using namespace compoGM_thread;
using DUse = Use<Value<double>>;

struct M2 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
                         map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
                         map<string, double>& size_factors) {
        m.component<NMatrix<OrphanNode<Normal>>>("log10(q)", genes, conditions, 1, 3, 1.5);
        m.component<NMatrix<DeterministicUnaryNode<double>>>("q", genes, conditions,
                                                             [](double a) { return pow(10, a); })
            .connect<NMatrices1To1<DUse>>("a", "log10(q)");

        m.component<NArray<OrphanNode<Normal>>>("log10(alpha)", genes, 1, -2, 2);
        m.component<NArray<DeterministicUnaryNode<double>>>(
             "1/alpha", genes, [](double a) { return 1. / double(pow(10, a)); })
            .connect<NArrays1To1<DUse>>("a", "log10(alpha)");

        m.component<NMatrix<BinaryNode<GammaShapeRate>>>("tau", genes, samples, 1)
            .connect<NArrays1To1<NArrayMultiprovide<DUse>>>("a", "1/alpha")
            .connect<NArrays1To1<NArrayMultiprovide<DUse>>>("b", "1/alpha");

        m.component<NArray<Constant<double>>>("sf", samples, 0)
            .connect<SetNArray<double>>("x", size_factors);

        m.component<NMatrix<DeterministicTernaryNode<double>>>(
             "lambda", genes, samples, [](double a, double b, double c) { return a * b * c; })
            .connect<NArrayMultiprovide<NArrays1To1<DUse>>>("a", "sf")
            .connect<NArrays1To1<NArraysMap<DUse>>>("b", "q", condition_mapping)
            .connect<NMatrices1To1<DUse>>("c", "tau");

        m.component<NMatrix<UnaryNode<Poisson>>>("K", genes, samples, 0)
            .connect<SetNMatrix<int>>("x", counts)
            .connect<NMatrices1To1<DUse>>("a", "lambda");
    }
};

void compute(int, char**) {
    Assembly assembly;
    {
        Model model;

        // Parsing data files
        auto counts = parse_counts("../data/rnaseq_mini/counts.tsv");
        auto samples = parse_samples("../data/rnaseq/samples.tsv");
        auto size_factors = parse_size_factors("../data/rnaseq/size_factors.tsv");
        check_consistency(counts, samples, size_factors);
        auto pgenes = partition(counts.genes, p);

        // graphical model
        std::cout << "-- Creating component model...\n";
        model.component<M2>("model", pgenes, samples.conditions, make_index_set(counts.samples),
                            counts.counts, samples.condition_mapping, size_factors.size_factors);

        // suffstats and metropolis hastings moves
        // model.component<NArray<GammaSuffstat>>("tau_suffstats", pgenes)
        //     .connect<NArrays1To1<NArrayMultiuse<DUse>>>("values", Address("model", "tau"))
        //     .connect<NArrays1To1<DUse>>("k", Address("model", "1/alpha"))
        //     .connect<NArrays1To1<DUse>>("theta", Address("model", "1/alpha"));

        model.component<NMatrix<SimpleMHMove<Shift>>>("move_q", pgenes, samples.conditions)
            .connect<NMatrices1To1<MoveToTarget<double>>>("target", Address("model", "log10(q)"))
            .connect<NArrays1To1<NArraysRevMap<Use<LogProb>>>>("logprob", Address("model", "K"),
                                                               samples.condition_mapping);

        model.component<NMatrix<SimpleMHMove<Scale>>>("move_tau", pgenes, samples.samples)
            .connect<NMatrices1To1<MoveToTarget<double>>>("target", Address("model", "tau"))
            .connect<NMatrices1To1<Use<LogProb>>>("logprob", Address("model", "K"));

        model.component<NArray<SimpleMHMove<Shift>>>("move_alpha", pgenes)
            .connect<NArrays1To1<MoveToTarget<double>>>("target", Address("model", "log10(alpha)"))
            .connect<NArrays1To1<NArrayMultiuse<Use<LogProb>>>>("logprob", Address("model", "tau"));

        // assembly
        std::cout << "-- Instantiating assembly...\n";
        assembly.instantiate_from(model);
    }

    std::cout << "-- Preparations before running chain\n";
    auto all_moves = assembly.get_all<Move>();
    auto all_watched = assembly.get_all<Value<double>>(
        std::set<Address>{Address("model", "log10(q)"), Address("model", "log10(alpha)")}, "model");

    // trace header
    auto trace = make_trace(all_watched, "tmp" + to_string(p.rank) + ".dat");
    trace.header();

    std::cout << "-- Running the chain\n";
    for (int iteration = 0; iteration < 500; iteration++) {
        for (int rep = 0; rep < 10; rep++) {
            for (auto&& move : all_moves.pointers()) {
                move->move(1.0);
                move->move(0.1);
                move->move(0.01);
            }
        }
        trace.line();
    }
    // for (auto&& move : all_moves.pointers()) {
    //     cerr << setprecision(3) << "Accept rate" << setw(40) << move->get_name() << "  -->  "
    //          << move->accept_rate() * 100 << "%" << endl;
    // }
}

int main(int argc, char** argv) {
    // int nb_threads = 1;
    // auto threads = spawn(0, nb_threads, compute);
    // join(threads);
    mpi_run(argc, argv, compute);
}
