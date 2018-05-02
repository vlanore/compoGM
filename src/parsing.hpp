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

#include <csv-parser.hpp>
#include <fstream>
#include <map>
#include <tinycompo.hpp>

using aria::csv::CsvParser;

/*
====================================================================================================
  ~*~ Counts parsing ~*~
==================================================================================================*/
struct CountParsingResult {
    std::map<std::string, std::map<std::string, int>> counts;
    std::vector<std::string> samples;  // list of samples in counts file
};

CountParsingResult parse_counts(std::string filename) {
    // Files and parsers
    std::ifstream file(filename);
    auto parser = CsvParser(file).delimiter('\t');

    // Result structure
    CountParsingResult result;

    // Counts array
    auto&& line = parser.begin();
    for (int i = 1; i < static_cast<int>(line->size()); ++i) {  // first line of counts file
        result.samples.push_back((*line)[i]);
        // cout << "Sample " << i << ": " << (*line)[i] << std::endl;
    }
    std::cerr << "-- Number of samples is " << result.samples.size() << std::endl;
    for (++line; line != parser.end(); ++line) {  // rest of the lines
        std::string gene = (*line)[0];
        for (int i = 1; i < static_cast<int>(line->size()); ++i) {
            result.counts[gene][result.samples.at(i - 1)] = stoi((*line)[i]);
        }
    }
    std::cerr << "-- Number of genes is " << result.counts.size() << std::endl;
    return result;
}

/*
====================================================================================================
  ~*~ Samples parsing ~*~
==================================================================================================*/
struct SamplesParsingResult {
    std::set<std::string> conditions;
    std::map<std::string, std::string> condition_mapping;  // sample -> condition
    std::set<std::string> samples;                         // set of samples in samples file
};

SamplesParsingResult parse_samples(std::string filename) {
    // Files and parsers
    std::ifstream file(filename);
    auto parser = CsvParser(file).delimiter('\t');

    // Result structure
    SamplesParsingResult result;

    for (auto line = ++parser.begin(); line != parser.end(); ++line) {
        result.conditions.insert((*line)[1]);
        result.samples.insert((*line)[0]);
        result.condition_mapping[(*line)[0]] = (*line)[1];
        // cout << (*line)[0] << ", " << (*line)[1] << endl;
    }
    std::cerr << "-- Number of conditions is " << result.conditions.size() << std::endl;
    return result;
}

/*
====================================================================================================
  ~*~ Checking consistency between counts and samples ~*~
==================================================================================================*/
void check_consistency(CountParsingResult counts, SamplesParsingResult samples) {
    // Checking that the two files samples identifiers match
    if (std::set<std::string>(counts.samples.begin(), counts.samples.end()) == samples.samples) {
        std::cerr << "-- List of samples in counts and samples match!\n";
    } else {
        std::cerr << "-- Mismatch between sample list in counts file (" << counts.samples.size()
                  << " samples) and samples files (" << samples.samples.size() << " samples)\n";
        exit(1);
    }
}
