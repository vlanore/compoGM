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

#pragma once

#include "distributions.hpp"
#include "partition.hpp"

IndexSet gen_indexset(std::string prefix, int amount) {
    IndexSet result;
    for (int i = 0; i < amount; i++) { result.insert(prefix + std::to_string(i)); }
    return result;
}

std::map<std::string, std::map<std::string, int>> gen_data(const IndexSet& s1, const IndexSet& s2) {
    std::map<std::string, std::map<std::string, int>> result;
    int i = 0;
    for (auto e_i : s1) {
        i++;
        int baseline_i = 10 + (i % 5);
        int j = 0;
        result[e_i] = {};
        for (auto e_j : s2) {
            j++;
            int value = baseline_i + (j % 5);
            result.at(e_i)[e_j] = value;
        }
    }
    return result;
}
