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

#include "mpi_proxies.hpp"
#include "partition.hpp"

template <class F, class... Args>
void mpi_run(int argc, char** argv, F f) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &compoGM::p.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &compoGM::p.size);
    compoGM::p.message("Started MPI process");
    f(argc, argv);
    compoGM::p.message("End of MPI process");
    MPI_Finalize();
}

int compoGM_mpi_tag = 0;

struct MasterSlaveConnect : tc::Meta {
    static void connect(tc::Model& m, tc::PortAddress port_master, tc::PortAddress port_slave,
        int slave_number, CE p = compoGM::p) {
        if (!p.rank) {  // master
            m.connect<tc::Set<MPIConnection>>(
                port_master, MPIConnection(slave_number, compoGM_mpi_tag));
        } else if (p.rank == slave_number) {
            m.connect<tc::Set<MPIConnection>>(port_slave, MPIConnection(0, compoGM_mpi_tag));
        }
        compoGM_mpi_tag++;
    }
};