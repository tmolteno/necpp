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
#include "common.h"      // nec_version macro

#include "cxxopts.hpp"

using namespace antlr4;

int main(int argc, char* argv[]) {
  // ---- Command-line options (cxxopts) ---------------------------------
  cxxopts::Options options(argv[0],
      "nec_parse — ANTLR 4 strict NEC-2 parser and simulation driver");
  options.add_options()
    ("i,input",   "input NEC file (reads from stdin if omitted)",
                   cxxopts::value<std::string>())
    ("o,output",  "output file (writes to stdout if omitted)",
                   cxxopts::value<std::string>())
    ("b,benchmark", "run the built-in NEC benchmark and exit")
    ("v,version", "print version and exit")
    ("h,help",    "print usage and exit")
    ;

  cxxopts::ParseResult opt;
  try {
    opt = options.parse(argc, argv);
  } catch (const cxxopts::exceptions::parsing& e) {
    std::cerr << "Error: " << e.what() << "\n\n"
              << options.help() << std::endl;
    return 1;
  }

  if (opt.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (opt.count("version")) {
    std::cout << "nec_parse " << nec_version << std::endl;
    return 0;
  }
  if (opt.count("benchmark")) {
    try {
      std::cout << "The nec2++ benchmark." << std::endl;
      std::cout << "nec_parse " << nec_version << std::endl << std::endl;
      nec_float score = nec_context::benchmark();
      std::cout << "Your computer's score is: " << score << " NEC's"
                << std::endl;
    } catch (const nec_exception& e) {
      std::cerr << "NEC++ Runtime Error:\n" << e.get_message() << std::endl;
      return 1;
    }
    return 0;
  }

  std::string input_file  = opt.count("input")  ? opt["input" ].as<std::string>() : "";
  std::string output_file = opt.count("output") ? opt["output"].as<std::string>() : "";

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
  //
  // Open a file if --output was given, else default to stdout.  The FILE*
  // is kept open for the lifetime of the program via a unique_ptr wrapper
  // so it is closed on any return path, including the success path below.
  FILE* out_fp = stdout;
  using file_closer_t = int (*)(FILE*);
  std::unique_ptr<FILE, file_closer_t> out_owner(stdout, static_cast<file_closer_t>(&fclose));
  if (!output_file.empty()) {
    FILE* f = fopen(output_file.c_str(), "w");
    if (!f) {
      std::cerr << "Cannot open output file: " << output_file << std::endl;
      return 1;
    }
    out_fp = f;
    out_owner.reset(f);
  }
  s_output.set_file(out_fp);
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
    //
    // As in nec2++, the structured summary goes to stdout while the detailed
    // NEC simulation log is directed to --output (or stdout) via nec_output.
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
