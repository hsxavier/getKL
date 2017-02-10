#ifndef OUTPUT_H
#define OUTPUT_H 1

#include "ParameterList.hpp"
#include "Utilities.hpp"

// Print Table to file according to configuration file keyword:
template <typename type, typename intT>
void OutputTable(const ParameterList & config, std::string keyword, type **table, intT Nrows, intT Ncols) {
  std::ofstream outfile;
  std::string filename;

  // Output if requested:
  if (config.reads(keyword)!="0") {
    filename = config.reads(keyword);
    Announce(">> Writing "+keyword+" to "+filename+"...");
    outfile.open(filename.c_str());
    if (!outfile.is_open()) error("OutputTable: cannot open file "+filename);
    PrintVecs(table, Nrows, Ncols, &outfile);
    outfile.close();
    Announce();
  }
}

#endif
