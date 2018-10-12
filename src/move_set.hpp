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
#include "tinycompo.hpp"

class MoveSet;

namespace compoGM {
    enum MoveType { scale, shift };
    enum DataType { integer, fp };


    struct _MoveDecl {
        std::string target_name;
        compoGM::MoveType move_type;
        compoGM::DataType data_type;
    };
}  // namespace compoGM

class MoveSet {
    tc::Model& model;
    tc::Address gm;
    std::vector<compoGM::_MoveDecl> moves;

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
    MoveSet(tc::Model& model, tc::Address gm) : model(model), gm(gm) {}

    void add(std::string target_name, compoGM::MoveType move_type,
        compoGM::DataType data_type = compoGM::fp) {
        moves.push_back({target_name, move_type, data_type});
    }

    void declare_moves() const {
        for (auto m : moves) {
            compoGM::p.message(
                "Adding move on %s in model %s", m.target_name.c_str(), gm.to_string().c_str());
            tc::Address target(gm, m.target_name);
            std::string move_name = "move_" + m.target_name;
            tc::PortAddress move_port("target", move_name);
            switch (m.move_type) {
                case compoGM::scale: adaptive_create<Scale>(move_name, target); break;
                case compoGM::shift: adaptive_create<Shift>(move_name, target); break;
            }
            switch (m.data_type) {
                case compoGM::integer:
                    model.connect<ConnectMove<int>>(move_port, gm, target);
                    break;
                case compoGM::fp: model.connect<ConnectMove<double>>(move_port, gm, target); break;
            }
        }
    }
};