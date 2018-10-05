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

#include <cstdio>
#include <string>
#include <vector>

// CE = Computing Entity
// abstraction that groups mpi processes and threads
struct CE {
    static const std::vector<int> colors;
    int rank{0}, size{0};

    template <class... Args>
    void message(const std::string& format, Args&&... args) {
        std::string color = "\e[0m\e[" + std::to_string(colors[rank % colors.size()]) + "m";
        std::string bold = "\e[0m\e[1m";
        std::string normal = "\e[0m";
        std::string format2 = bold + "[" + color + "%d" + bold + "/" + color + "%d" + bold + "] " +
                              normal + format + "\n";
        printf(format2.c_str(), rank, size, std::forward<Args>(args)...);
    }

    template <class... Args>
    [[noreturn]] void fail(const std::string& format, Args&&... args) {
        std::string color = "\e[0m\e[" + std::to_string(colors[rank % colors.size()]) + "m";
        std::string bold = "\e[0m\e[1m";
        std::string redbold = "\e[0m\e[1m\e[31m";
        std::string normal = "\e[0m";
        std::string format2 = bold + "[" + color + "%d" + bold + "/" + color + "%d" + bold + "] " +
                              redbold + "Error: " + normal + format + "\n";
        printf(format2.c_str(), rank, size, std::forward<Args>(args)...);
        exit(1);
    }
};

const std::vector<int> CE::colors{31, 32, 33, 34, 35, 36, 91, 92, 93, 94, 95, 96};

namespace compoGM {
thread_local CE p;
}  // namespace compoGM

void set_ce(int rank, int size) {
    compoGM::p.rank = rank;
    compoGM::p.size = size;
}