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

#include "chrono.hpp"
#include "gm_connectors.hpp"
#include "introspection.hpp"
#include "mcmc_moves.hpp"
#include "moves.hpp"
#include "suffstats.hpp"
#include "tinycompo.hpp"
#include "trace.hpp"
using tc::Use;

class MCMC;

namespace compoGM {
    enum MoveType { scale, shift };
    enum DataType { integer, fp };
    enum SuffstatType { gamma_ss, gamma_sr, poisson };

    struct _MoveDecl {
        tc::Address target;
        compoGM::MoveType move_type;
        compoGM::DataType data_type;
    };

    struct _SuffstatDecl {
        tc::Address target;
        std::set<tc::Address> affected_moves;
        SuffstatType type;
    };
}  // namespace compoGM

class MCMC {
    tc::Model& model;
    tc::Address gm;
    std::vector<compoGM::_MoveDecl> moves;
    std::vector<compoGM::_SuffstatDecl> suffstats;
    std::map<tc::Address, std::pair<tc::Address, tc::Address>> ss_usage;  // move->(target, ss)

    template <class MoveType>
    void adaptive_create(tc::Address move_address, tc::Address target) const {
        if (is_matrix(target, model)) {
            auto indices = get_matrix_indices(target, model);
            model.component<Matrix<SimpleMHMove<MoveType>>>(
                move_address, indices.first, indices.second);
        } else if (is_array(target, model)) {
            auto indices = get_array_indices(target, model);
            model.component<Array<SimpleMHMove<MoveType>>>(move_address, indices);
        } else {
            model.component<SimpleMHMove<MoveType>>(move_address);
        }
    }

  public:
    MCMC(tc::Model& model, tc::Address gm) : model(model), gm(gm) {}

    void move(tc::Address target, compoGM::MoveType move_type,
        compoGM::DataType data_type = compoGM::fp) {
        moves.push_back({target, move_type, data_type});
    }

    void suffstat(
        tc::Address target, std::set<tc::Address> affected_moves, compoGM::SuffstatType type) {
        suffstats.push_back({target, affected_moves, type});
        for (auto m : affected_moves) { ss_usage[m] = {target, target.last() + "_suffstats"}; }
    }

    void declare_suffstat(tc::Address target, std::set<tc::Address> affected_moves,
        compoGM::SuffstatType type) const {
        compoGM::p.message("Adding sufftsat on %s in model %s", target.c_str(), gm.c_str());

        tc::Address ss_address(target.to_string("-") + "_suffstats");
        tc::Address target_glob(gm, target);
        switch (type) {  // declaring suffstat component
            case compoGM::gamma_sr: model.component<GammaShapeRateSuffstat>(ss_address); break;
            case compoGM::gamma_ss: model.component<GammaShapeScaleSuffstat>(ss_address); break;
            case compoGM::poisson: model.component<PoissonSuffstat>(ss_address); break;
        }

        // computing target's parents
        tc::Introspector i(model.get_composite(gm));
        auto edges = i.directed_binops();
        std::map<std::string, tc::Address> parents;  // key is the portname associated to the parent
        for (auto edge : edges) {
            if (target.is_ancestor(edge.first.address)) {
                parents[edge.first.prop] = tc::Address(gm, edge.second);
            }
        }
        for (auto p : parents) {
            compoGM::p.message(
                "Parent %s of %s is %s", p.first.c_str(), target.c_str(), p.second.c_str());
        }

        switch (type) {
            case compoGM::gamma_sr:
            case compoGM::gamma_ss:
                model.connect<AdaptiveOneToMany<Use<Value<double>>>>(
                    tc::PortAddress("values", ss_address), target_glob);
                break;
            case compoGM::poisson:
                model.connect<AdaptiveOneToMany<Use<Value<int>>>>(
                    tc::PortAddress("values", ss_address), target_glob);
                break;
        }
        for (auto p : parents) {
            model.connect<Use<Value<double>>>(tc::PortAddress(p.first, ss_address), p.second);
        }
        for (auto am : affected_moves) {
            tc::Address move_address(am.to_string("-") + "_move");
            model.connect<DirectedLogProb>(
                tc::PortAddress("logprob", move_address), ss_address, LogProbSelector::Full);
        }
    }

    void declare_move(tc::Address target, compoGM::MoveType move_type,
        compoGM::DataType data_type = compoGM::fp) const {
        compoGM::p.message("Adding move on %s in model %s", target.c_str(), gm.c_str());
        tc::Address target_glob(gm, target);
        tc::Address move_address(target.to_string("-") + "_move");
        tc::PortAddress mp("target", move_address);

        std::set<tc::Address> used_ss;
        for (auto ss : ss_usage) {
            if (ss.first == target) { used_ss.insert(ss.second.first); }
        }

        switch (move_type) {
            case compoGM::scale: adaptive_create<Scale>(move_address, target_glob); break;
            case compoGM::shift: adaptive_create<Shift>(move_address, target_glob); break;
        }
        switch (data_type) {
            case compoGM::integer:
                model.connect<ConnectMove<int>>(mp, gm, target_glob, used_ss);
                break;
            case compoGM::fp:
                model.connect<ConnectMove<double>>(mp, gm, target_glob, used_ss);
                break;
        }
    }

    void declare_moves() const {
        for (auto m : moves) { declare_move(m.target, m.move_type, m.data_type); }
        for (auto s : suffstats) { declare_suffstat(s.target, s.affected_moves, s.type); }
    }

    void go(int nb_iterations, int nb_rep) const {
        compoGM::p.message("Instantiating component assembly");
        tc::Assembly a(model);

        compoGM::p.message("Setting up trace");
        std::set<tc::Address> all_moved;
        for (auto m : moves) { all_moved.insert(tc::Address(gm, m.target)); }
        auto trace = make_trace(a.get_all<Value<double>>(all_moved), "tmp.dat");
        trace.header();

        // set of all moves, used to determine which ones are not covered by suffstats
        std::set<tc::Address> all_moves;
        for (auto m : moves) { all_moves.insert(tc::Address(m.target.to_string("-") + "_move")); }

        // debug
        std::stringstream schedule;

        compoGM::p.message("Gathering pointers to moves and suff stats");
        std::map<tc::Address, std::pair<Proxy*, std::vector<Move*>>> pointersets;
        for (auto ss : suffstats) {
            schedule << "\t* gather suff stats for " << ss.target
                     << "\n\t* perfom the following moves " << nb_rep << " times: ";
            pointersets[ss.target].first = &a.at<Proxy>(ss.target.to_string("-") + "_suffstats");
            for (auto m : ss.affected_moves) {
                tc::Address move_address(m.to_string("-") + "_move");
                schedule << move_address << " ";

                all_moves.erase(move_address);
                auto ptrs = a.get_all<Move>(std::set<tc::Address>{move_address}).pointers();
                auto& entry = pointersets.at(ss.target).second;
                entry.insert(entry.end(), ptrs.begin(), ptrs.end());
            }
            schedule << "\n\t* release suff stats for " << ss.target << "\n";
        }
        auto other_moves = a.get_all<Move>(all_moves).pointers();
        schedule << "\t* perfom the following moves: ";
        for (auto m : all_moves) { schedule << m << " "; }
        compoGM::p.message("Move schedule is:\n%s", schedule.str().c_str());

        compoGM::p.message("Starting MCMC chain for %d iterations", nb_iterations);
        Chrono total_time;
        for (int iteration = 0; iteration < nb_iterations; iteration++) {
            for (auto ps : pointersets) {
                ps.second.first->acquire();
                for (int rep = 0; rep < nb_rep; rep++) {
                    for (auto m : ps.second.second) {
                        m->move(1.0);
                        m->move(0.1);
                        m->move(0.01);
                    }
                }
                ps.second.first->release();
            }
            for (auto m : other_moves) {
                m->move(1.0);
                m->move(0.1);
                m->move(0.01);
            }
            trace.line();
        }
        double elapsed_time = total_time.end();
        compoGM::p.message("MCMC chain has finished in %fms (%fms/iteration)", elapsed_time,
            elapsed_time / nb_iterations);
    }
};