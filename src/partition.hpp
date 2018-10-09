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

#include <map>
#include <set>
#include <string>
#include <vector>
#include "computing_entity.hpp"

// ================================================================================================
// classes and data types
using Index = std::string;
using IndexSet = std::set<Index>;
using IndexMapping = std::map<Index, Index>;

class Partition {
    size_t _offset, _size;
    std::vector<IndexSet> partition;

  public:
    Partition(IndexSet indexes, size_t size, size_t offset = 0) : _offset(offset), _size(size) {
        size_t nb_indexes = indexes.size();
        for (size_t i = 0; i < _size; i++) {
            auto begin = indexes.begin();
            std::advance(begin, i * nb_indexes / size);
            auto end = indexes.begin();
            std::advance(end, (i + 1) * nb_indexes / size);
            partition.emplace_back(begin, end);
            compoGM::p.message("Partition %d contains %d elements", i, partition.back().size());
        }
    }

    IndexSet get_partition(int i) const {
        int index = i - _offset;
        if (index >= 0 and index < int(_size)) {
            return partition.at(index);
        } else {  // if not in partition, returning everything
            IndexSet result;
            for (auto subpartition : partition) {
                result.insert(subpartition.begin(), subpartition.end());
            }
            return result;
        }
    }
    IndexSet my_partition() const { return get_partition(compoGM::p.rank); }

    size_t partition_size(int i) const {
        int index = i - _offset;
        if (index >= 0 and index < int(_size)) {
            return partition.at(i - _offset).size();
        } else {
            return partition_size_sum();
        }
    }
    size_t my_partition_size() const { return partition_size(compoGM::p.rank); }
    size_t partition_size_sum() const {
        return std::accumulate(partition.begin(), partition.end(), 0,
            [](int acc, IndexSet r) { return acc + r.size(); });
    }

    size_t size() const { return _size; }

    size_t max_partition_size() const {
        return std::accumulate(partition.begin(), partition.end(), 0,
            [](size_t max, IndexSet s) { return s.size() > max ? s.size() : max; });
    }

    int owner(Index index) const {
        for (size_t i = 0; i < _size; i++) {
            auto subpartition = partition.at(i);
            if (subpartition.find(index) != subpartition.end()) {
                // compoGM::p.message("Owner of %s is %d", index.c_str(), i+offset);
                return i + _offset;
            }
        }
        return -1;
    }
};

// ================================================================================================
// object creation and conversion
IndexSet make_index_set(std::vector<Index>& v) { return IndexSet(v.begin(), v.end()); }