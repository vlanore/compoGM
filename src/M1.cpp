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

#include <map>
#include <tinycompo.hpp>
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
    static void contents(Model& model, map<string, map<string, int>>& counts,
                         std::set<string>& conditions, map<string, string>& condition_mapping) {
        for (auto&& gene : counts) {
            model.composite(gene.first);
            Model& gene_composite = model.get_composite(gene.first);
            for (auto&& condition : conditions) {  // lambda
                gene_composite.component<OrphanNode<Normal>>("lambda_" + condition, 1, 2, 2);
            }
            for (auto&& count : gene.second) {  // K (counts)
                gene_composite.component<UnaryNode<Poisson, int>>("K_" + count.first, count.second)
                    .connect<Use<Value<double>>>("parent",
                                                 "lambda_" + condition_mapping.at(count.first));
            }
        }
    }
};

int main() {
    Model model;

    // Parsing data files
    auto counts = parse_counts("/home/vlanore/git/data/rnaseq/mini-counts.tsv");
    auto samples = parse_samples("/home/vlanore/git/data/rnaseq/samples.tsv");
    check_consistency(counts, samples);

    // graphical model
    model.composite<M1>("model", counts.counts, samples.conditions, samples.condition_mapping);

    // suffstats and metropolis hastings moves
    for (auto&& gene : counts.counts) {
        Model& gene_composite = model.get_composite("model").get_composite(gene.first);
        for (auto&& condition : samples.conditions) {  // creating moves connected to their targets
            gene_composite.component<PoissonSuffstat>("poissonsuffstat_" + condition)
                .connect<Use<Value<double>>>("lambda", "lambda_" + condition);

            gene_composite.component<SimpleMHMove<Scale, double>>("move_lambda_" + condition)
                .connect<Use<Value<double>>>("target", "lambda_" + condition)
                .connect<Use<Backup>>("targetbackup", "lambda_" + condition)
                .connect<Use<LogProb>>("logprob", "lambda_" + condition)
                .connect<Use<LogProb>>("logprob", "poissonsuffstat_" + condition);
        }
        for (auto&& sample : gene.second) {  // connecting moves to all children in model
            string condition = samples.condition_mapping.at(sample.first);
            // gene_composite.connect<Use<LogProb>>(PortAddress("logprob", "move_lambda_" +
            // condition), Address("K_" + sample.first));
            gene_composite.connect<Use<Value<int>>>(
                PortAddress("values", "poissonsuffstat_" + condition),
                Address("K_" + sample.first));
        }
    }

    model.dot_to_file();
    // model.get_composite("model").get_composite("HRA1").dot_to_file();

    // assembly
    Assembly assembly(model);

    // preparations before running
    vector<Value<double>*> all_lambdas;  // list of lambda nodes for registration
    vector<Go*> all_moves;
    ofstream output("tmp.dat");
    for (auto&& gene : counts.counts) {
        for (auto&& condition : samples.conditions) {
            output << "lambda_" + gene.first + "_" + condition + "\t";  // trace header
            all_lambdas.push_back(
                &assembly.at<Value<double>>(Address("model", gene.first, "lambda_" + condition)));
            all_moves.push_back(
                &assembly.at<Go>(Address("model", gene.first, "move_lambda_" + condition)));
        }
    }
    output << endl;

    // running the chain
    for (int iteration = 0; iteration < 10000; iteration++) {
        for (auto&& move : all_moves) {
            move->go();
        }
        for (auto&& lambda : all_lambdas) {
            output << lambda->get_ref() << "\t";
        }
        output << endl;
    }
}
