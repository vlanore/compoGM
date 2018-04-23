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

#include <tinycompo.hpp>
#include "distributions.hpp"
#include "interfaces.hpp"
#include "mcmc_moves.hpp"
#include "node_skeletons.hpp"
#include "suffstats.hpp"

using namespace std;
using namespace tc;

template <class Interface>
struct IJtoIK_meta {
    static void connect(Model& m, PortAddress user, Address provider, int size_i, int size_j,
                        int nb_cond, const vector<int>& mapping) {
        for (int i = 0; i < size_i; i++) {
            for (int j = 0; j < size_j; j++) {
                int cond_j = mapping.at(i);
                int index_user = i * size_j + j;
                int index_provider = i * nb_cond + cond_j;
                m.connect<Use<Interface>>(PortAddress(user.prop, Address(user.address, index_user)),
                                          Address(provider, index_provider));
            }
        }
    }
};

struct M1 : public Composite {
    static void contents(Model& model, int i, int j, int nb_cond, const vector<int>& conditions) {
        model.composite<Array<BinaryNode<Normal>>>("lambda", i * nb_cond, 1);

        // model.composite<Array<Array<UnaryNode<Poisson, int>>>>("K", i, j, -1);
        model.composite<Array<UnaryNode<Poisson, int>>>("K", i * j, -1)
            .connect<ArraySet<int>>("x", vector<int>{1, 2, 1, 3, 4, 5, 6, 1, 8, 9, 1, 2});

        model.meta_connect<IJtoIK_meta<Value<double>>>(PortAddress("parent", "K"), "lambda", i, j,
                                                       nb_cond, conditions);
    }
};

int main() {
    Model model;

    // input data
    int i = 4, j = 3, nb_cond = 2;
    vector<int> conditions{0, 1, 0, 1};

    // graphical model
    model.composite<M1>("model", i, j, nb_cond, conditions);

    model.perform_meta();
    model.dot_to_file();

    // assembly
    Assembly assembly(model);
}
