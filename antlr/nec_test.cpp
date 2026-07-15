/*
  Test driver for the ANTLR 4 full NEC parser.

  Reads a complete NEC input file and runs it through the ANTLR 4 parser
  and NecBuildVisitor, which dispatches geometry cards to c_geometry and
  program cards via the existing card handler table.

  Build (inside Docker):
      antlr4 -Dlanguage=Cpp -o generated -visitor NECFull.g4
      g++ -std=c++17 ... -o nec_parse
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>

#include "antlr4-runtime.h"
#include "NECFullLexer.h"
#include "NECFullParser.h"
#include "nec_error_listener.h"
#include "nec_visitor.h"

#include "nec_context.h"
#include "c_geometry.h"
#include "nec_exception.h"

using namespace antlr4;

int main(int argc, char* argv[]) {
  std::string input_file;
  if (argc > 1) input_file = argv[1];

  // Read entire file into a string (needed for error listener source lookup)
  std::string source;
  std::unique_ptr<CharStream> input;
  if (input_file.empty()) {
    std::string line;
    while (std::getline(std::cin, line))
      source += line + "\n";
    input.reset(new ANTLRInputStream(source));
  } else {
    std::ifstream f(input_file);
    if (!f.is_open()) {
      std::cerr << "Cannot open: " << input_file << std::endl;
      return 1;
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    source = ss.str();
    input.reset(new ANTLRInputStream(source));
  }

  NECFullLexer lexer(input.get());
  CommonTokenStream tokens(&lexer);
  NECFullParser parser(&tokens);

  // Custom error listener with source-line lookup
  auto* err = new NecErrorListener(
      input_file.empty() ? "<stdin>" : input_file);
  err->setSource(source);
  lexer.removeErrorListeners();
  lexer.addErrorListener(err);
  parser.removeErrorListeners();
  parser.addErrorListener(err);

  nec_context context;
  nec_output_file s_output;
  nec_output_flags s_output_flags;
  // Connect the output sink.  Must use set_file (a FILE*) rather than
  // set_stream: nec_context::set_output caches m_output_fp via get_fp(),
  // and the internal fast-print routines (nec_printf/integer/real_out)
  // short-circuit when m_output_fp == NULL, dropping printf-style output.
  // Without this, all simulation output (freq header, currents, radiation
  // tables) is silently discarded.
  s_output.set_file(stdout);
  context.set_output(s_output, s_output_flags);

  c_geometry* geo = context.get_geometry();
  geo->set_context(&context);

  try {
    auto* tree = parser.necFile();

    if (parser.getNumberOfSyntaxErrors() > 0) {
      std::cerr << parser.getNumberOfSyntaxErrors()
                << " parse error(s)." << std::endl;
      return 1;
    }

    NecBuildVisitor visitor;
    visitor.nec = &context;
    visitor.geo = geo;
    tree->accept(&visitor);

    // calc_prepare() ran inside visitGeCard (see nec_visitor.h), so the
    // interaction buffers are sized before any program card was dispatched.
    // Emit the accumulated structured results (radiation patterns, antenna
    // input/impedance, near fields, currents) computed by RP/XQ/NE/NH/PT.
    // nec_results::write() clears each result's write_file flag as it goes,
    // so a single call here is correct and idempotent.
    context.write_results(std::cout);

    std::cout << "Parsed successfully. Wires: " << visitor.nwire
              << ", Segments: " << geo->n_segments << std::endl;
    return 0;

  } catch (const nec_exception& e) {
    std::cerr << "Error: " << e.get_message() << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
