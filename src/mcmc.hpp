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
        std::string target_name;
        compoGM::MoveType move_type;
        compoGM::DataType data_type;
    };

    struct _SuffstatDecl {
        std::string target_name;
        std::set<std::string> affected_moves;
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
    void adaptive_create(std::string move_name, tc::Address target) const {
        if (is_matrix(target, model)) {
            auto& tc = model.get_composite(target);
            auto raw_indices_x = tc.all_component_names(0, true);
            auto indices_x = make_index_set(raw_indices_x);
            auto raw_indices_y = tc.get_composite(raw_indices_x.front()).all_component_names();
            auto indices_y = make_index_set(raw_indices_y);
            model.component<Matrix<SimpleMHMove<MoveType>>>(move_name, indices_x, indices_y);
        } else if (is_array(target, model)) {
            auto raw_indices = model.get_composite(target).all_component_names();
            auto indices = make_index_set(raw_indices);
            model.component<Array<SimpleMHMove<MoveType>>>(move_name, indices);
        } else {
            model.component<SimpleMHMove<MoveType>>(move_name);
        }
    }

  public:
    MCMC(tc::Model& model, tc::Address gm) : model(model), gm(gm) {}

    void move(std::string target_name, compoGM::MoveType move_type,
        compoGM::DataType data_type = compoGM::fp) {
        moves.push_back({target_name, move_type, data_type});
    }

    void suffstat(
        std::string target_name, std::set<std::string> affected_moves, compoGM::SuffstatType type) {
        suffstats.push_back({target_name, affected_moves, type});
        for (auto m : affected_moves) { ss_usage[m] = {target_name, target_name + "_suffstats"}; }
    }

    void declare_suffstat(std::string target_name, std::set<std::string> affected_moves,
        compoGM::SuffstatType type) const {
        compoGM::p.message(
            "Adding sufftsat on %s in model %s", target_name.c_str(), gm.to_string().c_str());

        auto ssname = target_name + "_suffstats";
        tc::Address target(gm, target_name);
        switch (type) {  // declaring suffstat component
            case compoGM::gamma_sr: model.component<GammaShapeRateSuffstat>(ssname); break;
            case compoGM::gamma_ss: model.component<GammaShapeScaleSuffstat>(ssname); break;
            case compoGM::poisson: model.component<PoissonSuffstat>(ssname); break;
        }

        // computing target's parents
        tc::Introspector i(model.get_composite(gm));
        auto target_in_gm = tc::Address(target_name);
        auto edges = i.directed_binops();
        std::map<std::string, tc::Address> parents;
        for (auto edge : edges) {
            if (target_in_gm.is_ancestor(edge.first.address)) {
                parents[edge.first.prop] = tc::Address(gm, edge.second);
            }
        }
        for (auto p : parents) {
            compoGM::p.message("Parent %s of %s is %s", p.first.c_str(), target.to_string().c_str(),
                p.second.to_string().c_str());
        }

        switch (type) {
            case compoGM::gamma_sr:
            case compoGM::gamma_ss:
                model.connect<AdaptiveOneToMany<Use<Value<double>>>>(
                    tc::PortAddress("values", ssname), target);
                break;
            case compoGM::poisson:
                model.connect<AdaptiveOneToMany<Use<Value<int>>>>(
                    tc::PortAddress("values", ssname), target);
                break;
        }
        for (auto p : parents) {
            model.connect<Use<Value<double>>>(tc::PortAddress(p.first, ssname), p.second);
        }
        for (auto am : affected_moves) {
            model.connect<DirectedLogProb>(
                tc::PortAddress("logprob", "move_" + am), ssname, LogProbSelector::Full);
        }
    }

    void declare_move(std::string target_name, compoGM::MoveType move_type,
        compoGM::DataType data_type = compoGM::fp) const {
        compoGM::p.message(
            "Adding move on %s in model %s", target_name.c_str(), gm.to_string().c_str());
        tc::Address target(gm, target_name);
        std::string move_name = "move_" + target_name;
        tc::PortAddress mp("target", move_name);

        std::map<tc::Address, tc::Address> use_ss;
        for (auto ss : ss_usage) {
            if (ss.first == target_name) {
                use_ss[tc::Address(ss.second.first)] = ss.second.second;
            }
        }

        switch (move_type) {
            case compoGM::scale: adaptive_create<Scale>(move_name, target); break;
            case compoGM::shift: adaptive_create<Shift>(move_name, target); break;
        }
        switch (data_type) {
            case compoGM::integer: model.connect<ConnectMove<int>>(mp, gm, target, use_ss); break;
            case compoGM::fp: model.connect<ConnectMove<double>>(mp, gm, target, use_ss); break;
        }
    }

    void declare_moves() const {
        for (auto m : moves) { declare_move(m.target_name, m.move_type, m.data_type); }
        for (auto s : suffstats) { declare_suffstat(s.target_name, s.affected_moves, s.type); }
    }

    void go(int nb_iterations, int nb_rep) const {
        tc::Assembly a(model);

        // trace
        std::set<tc::Address> all_moved;
        for (auto m : moves) { all_moved.insert(tc::Address(gm, m.target_name)); }
        auto trace = make_trace(a.get_all<Value<double>>(all_moved), "tmp.dat");
        trace.header();

        // gathering pointers to everything
        std::set<tc::Address> all_moves;
        for (auto m : moves) { all_moves.insert(tc::Address("move_" + m.target_name)); }
        std::map<std::string, std::pair<Proxy*, std::vector<Move*>>> pointersets;
        for (auto ss : suffstats) {
            pointersets[ss.target_name].first = &a.at<Proxy>(ss.target_name + "_suffstats");
            for (auto m : ss.affected_moves) {
                all_moves.erase(tc::Address("move_" + m));
                auto ptrs = a.get_all<Move>(std::set<tc::Address>{"move_" + m}).pointers();
                auto& entry = pointersets.at(ss.target_name).second;
                entry.insert(entry.end(), ptrs.begin(), ptrs.end());
            }
        }
        auto other_moves = a.get_all<Move>(all_moves).pointers();

        // go!
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
    }
};