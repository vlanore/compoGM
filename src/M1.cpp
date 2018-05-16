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

#include <future>
#include <iomanip>
#include <map>
#include <thread>
#include <tinycompo.hpp>
#include "connectors.hpp"
#include "distributions.hpp"
#include "interfaces.hpp"
#include "mcmc_moves.hpp"
#include "moves.hpp"
#include "node_skeletons.hpp"
#include "parsing.hpp"
#include "suffstats.hpp"

using namespace std;
using namespace tc;

struct M1 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
                         map<string, map<string, int>>& counts,
                         map<string, string>& condition_mapping) {
        m.component<NMatrix<OrphanNode<Normal>>>("lambda", genes, conditions, 1, 3, pow(1.5, 2));

        m.component<NMatrix<DeterministicUnaryNode<double>>>("exp", genes, conditions,
                                                             [](double x) { return pow(10, x); })
            .connect<NMatrices1To1<Use<Value<double>>>>("a", "lambda");

        m.component<NMatrix<UnaryNode<Poisson>>>("K", genes, samples, 0)
            .connect<SetNMatrix<int>>("x", counts)
            .connect<NArrays1To1<NArraysMap<Use<Value<double>>>>>("a", "exp", condition_mapping);
    }
};

int main() {
    Model model;

    // Parsing data files
    auto counts = parse_counts("../data/rnaseq_mini/counts.tsv");
    auto samples = parse_samples("../data/rnaseq/samples.tsv");
    check_consistency(counts, samples);

    // graphical model
    model.component<M1>("model", counts.genes, samples.conditions, make_index_set(counts.samples),
                        counts.counts, samples.condition_mapping);

    // suffstats and metropolis hastings moves
    model.component<NMatrix<PoissonSuffstat>>("poissonsuffstats", counts.genes, samples.conditions)
        .connect<NMatrices1To1<Use<Value<double>>>>("lambda", Address("model", "exp"))
        .connect<NArrays1To1<NArraysMap<Use<Value<int>>>>>(
            "values", Address("model", "K"), reverse_index_mapping(samples.condition_mapping));

    model.component<NMatrix<SimpleMHMove<Scale>>>("moves", counts.genes, samples.conditions, 0.1)
        .connect<NMatrices1To1<Use<Value<double>>>>("target", Address("model", "lambda"))
        .connect<NMatrices1To1<Use<LogProb>>>("logprob", "poissonsuffstats")
        .connect<NMatrices1To1<Use<LogProb>>>("logprob", Address("model", "lambda"));

    // for (auto&& gene : counts.counts) {
    //     Model& gene_composite = model.get_composite("model").get_composite(gene.first);
    //     for (auto&& condition : samples.conditions) {  // creating moves connected to their
    //     targets
    //         gene_composite.component<PoissonSuffstat>("poissonsuffstat_" + condition)
    //             .connect<Use<Value<double>>>("lambda", "exp_" + condition);

    //         gene_composite.component<SimpleMHMove<Scale>>("move_lambda_" + condition, 0.1)
    //             .connect<Use<Value<double>>>("target", "lambda_" + condition)
    //             .connect<Use<Backup>>("targetbackup", "lambda_" + condition)
    //             .connect<Use<LogProb>>("logprob", "poissonsuffstat_" + condition)
    //             .connect<Use<LogProb>>("logprob", "lambda_" + condition);
    //     }
    //     for (auto&& sample : gene.second) {  // connecting moves to all children in model
    //         string condition = samples.condition_mapping.at(sample.first);

    //         // gene_composite.connect<Use<LogProb>>(PortAddress("logprob", "move_lambda_" +
    //         // condition), Address("K_" + sample.first));

    //         gene_composite.connect<Use<Value<int>>>(
    //             PortAddress("values", "poissonsuffstat_" + condition),
    //             Address("K_" + sample.first));
    //     }
    // }

    // model.dot_to_file();
    // // model.get_composite("model").get_composite("HRA1").dot_to_file();

    // // assembly
    model.print();
    Assembly assembly(model);
    assembly.print_all();

    // // preparations before running
    // auto all_lambdas = assembly.get_all<OrphanNode<Normal>>();
    // auto all_moves = assembly.get_all<SimpleMHMove<Scale>>();

    // // trace header
    // ofstream output("tmp.dat");
    // for (auto&& gene : counts.counts) {
    //     for (auto&& condition : samples.conditions) {
    //         output << "lambda_" + gene.first + "_" + condition + "\t";  // trace header
    //     }
    // }
    // output << endl;

    // // running the chain
    // // vector<future<void>> futures(all_moves.size());
    // for (int iteration = 0; iteration < 5000; iteration++) {
    //     for (int rep = 0; rep < 10; rep++) {
    //         for (auto&& move : all_moves) {
    //             move->go();
    //         }

    //         // for (size_t i = 0; i< all_moves.size(); ++i) {
    //         //     futures.at(i) = async(&Go::go, all_moves.at(i));
    //         // }
    //         // for (auto&& future : futures) {
    //         //     future.get();
    //         // }
    //     }
    //     for (auto&& lambda : all_lambdas) {
    //         output << lambda->get_ref() << "\t";
    //     }
    //     output << endl;
    // }
    // for (auto&& move : all_moves) {
    //     cerr << setprecision(3) << "Accept rate" << setw(40) << move->get_name() << "  -->  "
    //          << move->accept_rate() * 100 << "%" << endl;
    // }
}
