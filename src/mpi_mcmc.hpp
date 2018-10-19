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

#include "mcmc.hpp"
#include "mpi_helpers.hpp"

class MpiMCMC : public MCMC {
  public:
    MpiMCMC(tc::Model& m, tc::Address gm) : MCMC(m, gm) {}

    template <class... Args>
    void master_add(Args&&... args) {
        if (compoGM::p.rank == 0) { move(std::forward<Args>(args)...); }
    }

    template <class... Args>
    void slave_add(Args&&... args) {
        if (compoGM::p.rank != 0) { move(std::forward<Args>(args)...); }
    }

    struct MoveExecData {
        std::vector<Move*> moves;
        int move_rep;
        double tuning_mult;
    };
    void go(int nb_iterations, int nb_rep_master, int nb_rep_slave, int burnin = 0) const {
        // instantiating assembly
        Assembly a(model);

        // gathering pointers and preparing trace
        std::vector<MoveExecData> move_exec;
        for (auto m : moves) {
            move_exec.emplace_back();
            move_exec.back().moves =
                a.get_all<Move>(
                     std::set<tc::Address>{tc::Address(m.target.to_string("-") + "_move")})
                    .pointers();
            move_exec.back().tuning_mult = m.tuning_mult;
            move_exec.back().move_rep = m.move_rep;
        }
        auto proxies = a.get_all<Proxy>().pointers();

        // main loop
        compoGM::p.message("Reaching go barrier");
        // MPI_Barrier(MPI_COMM_WORLD);
        compoGM::p.message("Go!");
        Chrono total_time, computing_time, acquire_time, release_time;
        // master ==================================================================================
        if (!compoGM::p.rank) {
            compoGM::p.message("Setting up trace");
            std::set<tc::Address> all_moved;
            for (auto m : MCMC::moves) { all_moved.insert(tc::Address(gm, m.target)); }
            std::string tracename =
                "trace_m3_" + std::to_string(compoGM::p.size) + "_processes.dat";
            auto trace = make_trace(a.get_all<Value<double>>(all_moved), tracename);
            trace.header();

            Chrono writing_time;
            for (int iteration = 0; iteration < nb_iterations; iteration++) {
                acquire_time.start();
                for (auto proxy : proxies) { proxy->acquire(); }
                acquire_time.end();
                computing_time.start();
                for (int i = 0; i < nb_rep_master; i++) {
                    for (auto move : move_exec) {
                        for (int r = 0; r < move.move_rep; r++) {
                            for (auto ptr : move.moves) {
                                ptr->move(move.tuning_mult * 1.0);
                                ptr->move(move.tuning_mult * 0.1);
                                ptr->move(move.tuning_mult * 0.01);
                            }
                        }
                    }
                }
                computing_time.end();
                release_time.start();
                for (auto proxy : proxies) { proxy->release(); }
                release_time.end();
                writing_time.start();
                if (iteration >= burnin) { trace.line(); }
                writing_time.end();
            }
            compoGM::p.message("Average writing time is %fms", writing_time.mean());
            // slaves ==============================================================================
        } else {
            for (auto proxy : proxies) { proxy->release(); }
            for (int iteration = 0; iteration < nb_iterations; iteration++) {
                acquire_time.start();
                for (auto proxy : proxies) { proxy->acquire(); }
                acquire_time.end();
                computing_time.start();
                for (int i = 0; i < nb_rep_slave; i++) {
                    for (auto move : move_exec) {
                        for (int r = 0; r < move.move_rep; r++) {
                            for (auto ptr : move.moves) {
                                ptr->move(move.tuning_mult * 1.0);
                                ptr->move(move.tuning_mult * 0.1);
                                ptr->move(move.tuning_mult * 0.01);
                            }
                        }
                    }
                }
                computing_time.end();
                release_time.start();
                for (auto proxy : proxies) { proxy->release(); }
                release_time.end();
            }
        }
        double elapsed_time = total_time.end();
        compoGM::p.message("MCMC chain has finished in %fms (%fms/iteration)", elapsed_time,
            elapsed_time / nb_iterations);
        compoGM::p.message("Average computing time is %fms", computing_time.mean());
        compoGM::p.message("Average acquire time is %fms", acquire_time.mean());
        compoGM::p.message("Average release time is %fms", release_time.mean());
    }
};