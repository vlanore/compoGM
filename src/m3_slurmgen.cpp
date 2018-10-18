/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2018).
Contributors:
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

This software is a component-based library to write bayesian inference programs based on the
graphical m.

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

#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage:\n\tm3_slurmgen <nb_nodes>\n";
        exit(1);
    }
    int nb_nodes = atoi(argv[1]);

    std::ofstream f("m3_" + std::to_string(nb_nodes) + "_nodes.slurm");
    f << "#!/bin/bash\n#SBATCH -J m3\n#SBATCH --nodes=" << nb_nodes
      << "\n#SBATCH --ntasks=" << 24 * nb_nodes
      << "\n#SBTACH --ntasks-per-node=24\n#SBATCH --threads-per-core=1\n#SBATCH "
         "--time=00:05:00\n#SBATCH --output m3_"
      << nb_nodes
      << "_nodes.output\n#SBATCH --constraint=HSW24\n\nmodule purge\nmodule load "
         "intel/18.1 gcc/6.2.0 openmpi/gnu/2.0.2\nsrun -n $SLURM_NTASKS "
         "~/code_directory/compoGM/M3_mpi_bin ~/rnaseq";
}