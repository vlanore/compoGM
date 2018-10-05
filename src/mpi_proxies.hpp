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

#include <mpi.h>
#include "interfaces.hpp"
#include "tinycompo.hpp"
using tc::Component;

struct MPIConnection {
    int target_process{-1};
    int tag{-1};
    MPIConnection() = default;
    MPIConnection(int target, int tag) : target_process(target), tag(tag) {}
};

// assuming value type is double
class ProbNodeProv : public Component, public Proxy {
    Value<double>* target;
    MPIConnection connection;

  public:
    ProbNodeProv() {
        port("target", &ProbNodeProv::target);
        port("connection", &ProbNodeProv::connection);
    }

    void acquire() override {}

    void release() override {
        double buffer = target->get_ref();
        MPI_Send(&buffer, 1, MPI_DOUBLE, connection.target_process, connection.tag, MPI_COMM_WORLD);
        // compoGM::p.message("Sent value %f to %d", buffer, connection.target_process);
    }
};

// assuming value type is double
class ProbNodeUse : public Component, public Proxy {
    Value<double>* target;
    MPIConnection connection;

  public:
    ProbNodeUse() {
        port("target", &ProbNodeUse::target);
        port("connection", &ProbNodeUse::connection);
    }

    void acquire() override {
        double buffer = -1;
        MPI_Recv(&buffer, 1, MPI_DOUBLE, connection.target_process, connection.tag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // compoGM::p.message("Received value %f from %d", buffer, connection.target_process);
        target->get_ref() = buffer;
    }

    void release() override {}
};

// assuming value type is double
class MasterBcast : public Component, public Proxy {
    std::vector<Value<double>*> targets;
    void add_target(Value<double>* ptr) { targets.push_back(ptr); }
    std::vector<double> data;

  public:
    MasterBcast() { port("target", &MasterBcast::add_target); }

    void acquire() override {}

    void release() override {
        int n = targets.size();
        data.reserve(n);
        for (auto target : targets) {
            data.push_back(target->get_ref());
        }
        MPI_Bcast(data.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
};

// assuming value type is double
class SlaveBcast : public Component, public Proxy {
    std::vector<Value<double>*> targets;
    void add_target(Value<double>* ptr) { targets.push_back(ptr); }
    std::vector<double> data;

  public:
    SlaveBcast() { port("target", &SlaveBcast::add_target); }

    void acquire() override {
        size_t n = targets.size();
        data.assign(n, -1);
        MPI_Bcast(data.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (size_t i = 0; i < n; i++) {
            targets[i]->get_ref() = data[i];
        }
    }

    void release() override {}
};

// assuming value type is double
class MasterGather : public Component, public Proxy {
    std::vector<Value<double>*> targets;
    void add_target(Value<double>* ptr) { targets.push_back(ptr); }
    std::vector<double> data;
    Partition partition;

  public:
    MasterGather(Partition partition) : partition(partition) {
        if (partition.size() != size_t(compoGM::p.size - 1)) {
            std::cerr << "MasterGather error: number of partitions (" << partition.size()
                      << ") doesn't match number of workerss (" << compoGM::p.size - 1 << ")\n";
            exit(1);
        }
        port("target", &MasterGather::add_target);
    }

    void acquire() override {
        if (partition.partition_size_sum() != targets.size()) {
            std::cerr << "MasterGather error: number of targets (" << targets.size()
                      << ") doesn't match number of elements in partition ("
                      << partition.partition_size_sum() << ")\n";
            exit(1);
        }
        size_t nb_partitions = partition.size();
        size_t subbuffer_size = partition.max_partition_size();
        size_t buffer_size = nb_partitions * subbuffer_size;
        data.assign(buffer_size, -1);
        MPI_Gather(NULL, 0, MPI_DOUBLE, data.data(), buffer_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        size_t target_index = 0;
        for (size_t partition_index = 0; partition_index < nb_partitions; partition_index++) {
            for (size_t element_index = 0;
                 element_index < partition.partition_size(partition_index + 1); element_index++) {
                targets.at(target_index)->get_ref() =
                    data.at(subbuffer_size * partition_index + element_index);
            }
        }
    }

    void release() override {}
};

// assuming value type is double
class WorkerGather : public Component, public Proxy {
    std::vector<Value<double>*> targets;
    void add_target(Value<double>* ptr) { targets.push_back(ptr); }
    std::vector<double> data;
    Partition partition;

  public:
    WorkerGather(Partition partition) : partition(partition) {
        if (partition.size() != size_t(compoGM::p.size - 1)) {
            std::cerr << "WorkerGather error: number of partitions (" << partition.size()
                      << ") doesn't match number of workers (" << compoGM::p.size - 1 << ")\n";
            exit(1);
        }
        port("target", &WorkerGather::add_target);
    }

    void acquire() override {
        size_t buffer_size = partition.max_partition_size();
        size_t my_size = partition.my_partition_size();
        if (my_size != targets.size()) {
            std::cerr << "WorkerGather error: number of targets (" << targets.size()
                      << ") doesn't match number of elements in my partition (" << my_size << ")\n";
            exit(1);
        }
        data.assign(buffer_size, -1);  // filling buffer with -1s
        for (size_t i = 0; i < my_size; i++) {
            data[i] = targets[i]->get_ref();
        }
        MPI_Gather(data.data(), buffer_size, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void release() override {}
};