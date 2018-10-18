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

struct M1 : public Composite {
    static void contents(Model& m, IndexSet& genes, IndexSet& conditions, IndexSet& samples,
        map<string, map<string, int>>& counts, map<string, string>& condition_mapping) {
        m.component<Matrix<OrphanNormal>>("log10(lambda)", genes, conditions, 1, 3, pow(1.5, 2));
        m.connect<MapPower10>("log10(lambda)", "lambda");

        m.component<Matrix<Poisson>>("K", genes, samples, 0)
            .connect<SetMatrix<int>>("x", counts)
            .connect<ManyToMany<ArraysMap<UseValue>>>("a", "lambda", condition_mapping);
    }
};

void compute(int argc, char** argv) {
    if (argc < 2) {
        cerr << "usage:\n\tM1_bin <data_location>\n";
        exit(1);
    }

    Model m;

    // Parsing data files
    string data_folder = argv[1];
    auto counts = parse_counts(data_folder + "/counts.tsv");
    auto samples = parse_samples(data_folder + "/samples.tsv");
    check_consistency(counts, samples);

    // graphical model
    m.component<M1>("model", counts.genes, samples.conditions, make_index_set(counts.samples),
        counts.counts, samples.condition_mapping);

    MCMC mcmc(m, "model");
    mcmc.move("log10(lambda)", scale);
    // TODO support array of suffstats and/or addresses for targets (instead of just strings)!
    // for (auto gene : counts.genes) {
    //     mcmc.suffstat("K__" + gene, {"log10(lambda)__" + gene}, poisson);
    // }
    mcmc.declare_moves();

    mcmc.go(5000, 10);
}

int main(int argc, char** argv) { compute(argc, argv); }
