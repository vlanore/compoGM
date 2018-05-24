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

using namespace std;
using namespace tc;
using DUse = Use<Value<double>>;

struct M2 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
                         map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
                         map<string, double>& size_factors) {
        m.component<NMatrix<OrphanNode<Normal>>>("q", genes, conditions, 1, -2, 2);
        m.component<NMatrix<DeterministicUnaryNode<double>>>("exp(q)", genes, conditions,
                                                             [](double a) { return pow(10, a); })
            .connect<NMatrices1To1<DUse>>("a", "q");

        m.component<NArray<OrphanNode<Normal>>>("alpha", genes, 1, 2, 2);
        m.component<NArray<DeterministicUnaryNode<double>>>(
             "1/exp(alpha)", genes, [](double a) { return 1. / double(pow(10, a)); })
            .connect<NArrays1To1<DUse>>("a", "alpha");

        m.component<NMatrix<BinaryNode<GammaShapeRate>>>("tau", genes, samples, 1)
            .connect<NArrays1To1<NArrayMultiprovide<DUse>>>("a", "1/exp(alpha)")
            .connect<NArrays1To1<NArrayMultiprovide<DUse>>>("b", "1/exp(alpha)");

        m.component<NArray<Constant<double>>>("sf", samples, 0)
            .connect<SetNArray<double>>("x", size_factors);

        m.component<NMatrix<DeterministicTernaryNode<double>>>(
             "lambda", genes, samples, [](double a, double b, double c) { return a * b * c; })
            .connect<NArrayMultiprovide<NArrays1To1<DUse>>>("a", "sf")
            .connect<NArrays1To1<NArraysMap<DUse>>>("b", "exp(q)", condition_mapping)
            .connect<NMatrices1To1<DUse>>("c", "tau");

        m.component<NMatrix<UnaryNode<Poisson>>>("K", genes, samples, 0)
            .connect<SetNMatrix<int>>("x", counts)
            .connect<NMatrices1To1<DUse>>("a", "lambda");
    }
};

int main() {
    Model model;

    // Parsing data files
    auto counts = parse_counts("../data/rnaseq_mini/counts.tsv");
    auto samples = parse_samples("../data/rnaseq/samples.tsv");
    auto size_factors = parse_size_factors("../data/rnaseq/size_factors.tsv");
    check_consistency(counts, samples, size_factors);

    // graphical model
    std::cout << "-- Creating component model...\n";
    model.component<M2>("model", counts.genes, samples.conditions, make_index_set(counts.samples),
                        counts.counts, samples.condition_mapping, size_factors.size_factors);

    // suffstats and metropolis hastings moves
    model.component<NMatrix<SimpleMHMove<Scale>>>("move_q", counts.genes, samples.conditions)
        .connect<NMatrices1To1<DUse>>("target", Address("model", "q"))
        .connect<NMatrices1To1<Use<Backup>>>("targetbackup", Address("model", "q"))
        .connect<NMatrices1To1<Use<LogProb>>>("logprob", Address("model", "q"))
        .connect<NArrays1To1<NArraysRevMap<Use<LogProb>>>>("logprob", Address("model", "K"),
                                                           samples.condition_mapping);

    model.component<NArray<SimpleMHMove<Scale>>>("move_alpha", counts.genes)
        .connect<NArrays1To1<DUse>>("target", Address("model", "alpha"))
        .connect<NArrays1To1<Use<Backup>>>("targetbackup", Address("model", "alpha"))
        .connect<NArrays1To1<Use<LogProb>>>("logprob", Address("model", "alpha"))
        .connect<NArrays1To1<NArrayMultiuse<Use<LogProb>>>>("logprob", Address("model", "K"));

    model.component<NMatrix<SimpleMHMove<Scale>>>("move_tau", counts.genes, samples.samples)
        .connect<NMatrices1To1<DUse>>("target", Address("model", "tau"))
        .connect<NMatrices1To1<Use<Backup>>>("targetbackup", Address("model", "tau"))
        .connect<NMatrices1To1<Use<LogProb>>>("logprob", Address("model", "tau"))
        .connect<NMatrices1To1<Use<LogProb>>>("logprob", Address("model", "K"));

    // model.dot_to_file();
    // model.print();
    // model.get_composite("model").get_composite("HRA1").dot_to_file();

    // assembly
    std::cout << "-- Instantiating assembly...\n";
    Assembly assembly(model);
    // assembly.print_all();

    std::cout << "-- Preparations before running chain\n";
    auto all_moves = assembly.get_all<SimpleMHMove<Scale>>();
    auto all_watched = assembly.get_all<Value<double>>(
        std::set<Address>{Address("model", "q"), Address("model", "alpha"),
                          Address("model", "tau")},
        "model");

    // trace header
    ofstream output("tmp.dat");
    for (auto&& name : all_watched.names()) {
        output << name.to_string() << '\t';
    }
    output << endl;

    std::cout << "-- Running the chain\n";
    for (int iteration = 0; iteration < 500; iteration++) {
        for (int rep = 0; rep < 10; rep++) {
            for (auto&& move : all_moves.pointers()) {
                move->move(1.0);
                move->move(0.1);
                move->move(0.01);
            }
        }
        for (auto&& pointer : all_watched.pointers()) {
            output << to_string(pointer->get_ref()) << '\t';
        }
        output << endl;
    }
    for (auto&& move : all_moves.pointers()) {
        cerr << setprecision(3) << "Accept rate" << setw(40) << move->get_name() << "  -->  "
             << move->accept_rate() * 100 << "%" << endl;
    }
}
