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

    void go(int nb_iterations, int nb_rep) const {
        // instantiating assembly
        Assembly a(model);

        // gathering pointers and preparing trace
        auto moves = a.get_all<Move>().pointers();
        auto proxies = a.get_all<Proxy>().pointers();

        // main loop
        compoGM::p.message("Go!");
        Chrono total_time;
        Chrono computing_time;
        // master ==================================================================================
        if (!compoGM::p.rank) {
            compoGM::p.message("Setting up trace");
            std::set<tc::Address> all_moved;
            for (auto m : MCMC::moves) { all_moved.insert(tc::Address(gm, m.target)); }
            std::string tracename = "trace_m3_" + std::to_string(compoGM::p.size) + "_processes.dat";
            auto trace = make_trace(a.get_all<Value<double>>(all_moved), tracename);
            trace.header();

            Chrono writing_time;
            for (int iteration = 0; iteration < nb_iterations; iteration++) {
                for (auto proxy : proxies) { proxy->acquire(); }
                computing_time.start();
                for (int i = 0; i < nb_rep; i++) {
                    for (auto move : moves) {
                        move->move(1.0);
                        move->move(0.1);
                        move->move(0.01);
                    }
                }
                computing_time.end();
                for (auto proxy : proxies) { proxy->release(); }
                writing_time.start();
                trace.line();
                writing_time.end();
            }
            double elapsed_time = total_time.end();
            compoGM::p.message("MCMC chain has finished in %fms (%fms/iteration)", elapsed_time,
                elapsed_time / nb_iterations);
            compoGM::p.message("Average computing time is %fms", computing_time.mean());
            compoGM::p.message("Average writing time is %fms", writing_time.mean());
            // slaves ==============================================================================
        } else {
            for (auto proxy : proxies) { proxy->release(); }
            for (int iteration = 0; iteration < nb_iterations; iteration++) {
                for (auto proxy : proxies) { proxy->acquire(); }
                computing_time.start();
                for (int i = 0; i < nb_rep; i++) {
                    for (auto move : moves) {
                        move->move(1.0);
                        move->move(0.1);
                        move->move(0.01);
                    }
                }
                computing_time.end();
                for (auto proxy : proxies) { proxy->release(); }
            }
            double elapsed_time = total_time.end();
            compoGM::p.message("MCMC chain has finished in %fms (%fms/iteration)", elapsed_time,
                elapsed_time / nb_iterations);
            compoGM::p.message("Average computing time is %fms", computing_time.mean());
        }
    }
};