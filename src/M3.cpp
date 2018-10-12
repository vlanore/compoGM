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

#include <iomanip>
#include "compoGM.hpp"
#include "partition.hpp"
#include "thread_helpers.hpp"

using namespace std;
using namespace compoGM;

struct M3 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
        map<string, map<string, int>>& counts, IndexMapping& condition_mapping,
        map<string, double>& size_factors) {
        m.component<Matrix<OrphanNormal>>("log10(q)", genes, conditions, 1, 3, 1.5);
        m.component<Matrix<Power10>>("q", genes, conditions)
            .connect<MatrixToValueMatrix>("a", "log10(q)");

        m.component<OrphanNormal>("log10(a)", 1, -2, 2);
        m.component<OrphanNormal>("log10(b)", 1, 0, 2);
        m.component<OrphanExp>("sigma", 1, 1);

        m.component<Array<Sum>>("q_bar", genes).connect<ArrayToValueMatrixLines>("parent", "q");

        m.component<Array<DeterministicTernaryNode<double>>>("log10(alpha_bar)", genes,
             [](double a, double b, double q_bar) {
                 return log10(pow(10, a) + pow(10, b) / q_bar);
             })
            .connect<ArrayToValue>("a", "log10(a)")
            .connect<ArrayToValue>("b", "log10(b)")
            .connect<ArrayToValueArray>("c", "q_bar");

        m.component<Array<Normal>>("log10(alpha)", genes, 1)
            .connect<ArrayToValueArray>("a", "log10(alpha_bar)")
            .connect<ArrayToValue>("b", "sigma");
        m.component<Array<DeterministicUnaryNode<double>>>(
             "1/alpha", genes, [](double a) { return 1. / double(pow(10, a)); })
            .connect<ArrayToValueArray>("a", "log10(alpha)");

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

    // graphical model
    m.component<M3>("model", counts.genes, samples.conditions, make_index_set(counts.samples),
        counts.counts, samples.condition_mapping, size_factors.size_factors);

    // suffstats and metropolis hastings moves
    MoveSet ms(m, "model");
    ms.add("log10(a)", shift);
    ms.add("log10(b)", shift);
    ms.add("sigma", scale);
    ms.add("log10(q)", shift);
    ms.add("tau", scale);
    ms.add("log10(alpha)", shift);
    ms.declare_moves();

    // assembly
    Assembly a(m);

    auto moves = a.get_all<Move>().pointers();
    auto all_watched = a.get_all<Value<double>>(
        std::set<Address>{Address("model", "log10(q)"), Address("model", "log10(alpha)"),
            Address("model", "log10(a)"), Address("model", "log10(b)"), Address("model", "sigma")},
        "model");

    auto trace = make_trace(all_watched, "tmp" + to_string(p.rank) + ".dat");
    trace.header();

    for (int iteration = 0; iteration < 5000; iteration++) {
        for (auto&& move : moves) {
            move->move(1.0);
            move->move(0.1);
            move->move(0.01);
        }
        trace.line();
    }
}

int main(int argc, char** argv) { compute(argc, argv); }
