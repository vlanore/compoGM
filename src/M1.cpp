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

#include <csv-parser.hpp>
#include <fstream>
#include <map>
#include <tinycompo.hpp>
#include "distributions.hpp"
#include "interfaces.hpp"
#include "mcmc_moves.hpp"
#include "node_skeletons.hpp"
#include "suffstats.hpp"

using namespace std;
using namespace tc;
using aria::csv::CsvParser;

struct M1 : public Composite {
    static void contents(Model& model, int nb_genes, int nb_samples, int nb_cond,
                         const vector<int>& conditions) {
        model.composite<Array<BinaryNode<Normal>>>("lambda", nb_genes * nb_cond, 1);

        // model.composite<Array<Array<UnaryNode<Poisson, int>>>>("K", i, j, -1);
        model.composite<Array<UnaryNode<Poisson, int>>>("K", nb_genes * nb_samples, -1)
            .connect<ArraySet<int>>("x", vector<int>{1, 2, 1, 3, 4, 5, 6, 1, 8, 9, 1, 2});

        for (int i = 0; i < nb_genes; i++) {
            for (int j = 0; j < nb_samples; j++) {
                int cond_j = conditions.at(j);
                int index_user = i * nb_samples + j;
                int index_provider = i * nb_cond + cond_j;
                model.connect<Use<Value<double>>>(PortAddress("parent", Address("K", index_user)),
                                                  Address("lambda", index_provider));
            }
        }
    }
};

int main() {
    Model model;

    // Files and parsers
    ifstream counts_file("/home/vlanore/git/data/rnaseq/counts.tsv");
    ifstream samples_file("/home/vlanore/git/data/rnaseq/samples.tsv");
    auto counts_parser = CsvParser(counts_file).delimiter(' ');
    auto samples_parser = CsvParser(samples_file).delimiter(' ');

    // Counts array
    vector<vector<int>> counts;  // bidim array of counts (first index is gene, second is sample)
    map<string, int> samples;    // sample -> sample index map
    map<string, int> genes;      // gene -> gene index map
    auto&& line = counts_parser.begin();
    for (int i = 1; i < static_cast<int>(line->size()); ++i) {  // first line of counts file
        samples.insert(make_pair((*line)[i], i));
    }
    int next_gene_index{0};
    for (++line; line != counts_parser.end(); ++line) {  // rest of the lines
        counts.emplace_back();
        auto&& field = line->begin();
        genes.insert(make_pair(*field, next_gene_index++));
        for (++field; field != line->end(); ++field) {
            counts.back().push_back(stoi(*field));
        }
    }

    // Conditions
    // for (auto& line : samples_parser) {
    // }

    int i = 4, j = 3, nb_cond = 2;
    vector<int> conditions{0, 1, 0};

    // graphical model
    model.composite<M1>("model", i, j, nb_cond, conditions);

    model.perform_meta();
    model.dot_to_file();

    // assembly
    Assembly assembly(model);
}
