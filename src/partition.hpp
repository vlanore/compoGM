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

#include "computing_entity.hpp"
#include "index_set.hpp"

IndexSet partition(IndexSet s, CE p) {
    auto begin = s.begin();
    std::advance(begin, p.rank * s.size() / p.size);
    auto end = s.begin();
    std::advance(end, (p.rank + 1) * s.size() / p.size);
    return IndexSet(begin, end);
}

IndexSet partition_slave(IndexSet s, CE p) {
    if (p.rank) {
        auto begin = s.begin();
        std::advance(begin, (p.rank - 1) * s.size() / (p.size - 1));
        auto end = s.begin();
        std::advance(end, p.rank * s.size() / (p.size - 1));
        return IndexSet(begin, end);
    } else {       // called on master
        return s;  // returns full set
    }
}

int index_owner_slave(Index i, IndexSet s, CE p) {
    int int_index = std::distance(s.begin(), s.find(i));
    return (int_index * (p.size - 1) / s.size()) + 1;
}