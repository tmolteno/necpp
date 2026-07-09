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
#include <memory>
#include <string>

#include "antlr4-runtime.h"
#include "NECFullLexer.h"
#include "NECFullParser.h"
#include "nec_visitor.h"

#include "nec_context.h"
#include "c_geometry.h"
#include "nec_exception.h"

using namespace antlr4;

int main(int argc, char* argv[]) {
  std::string input_file;
  if (argc > 1) input_file = argv[1];

  std::unique_ptr<CharStream> input;
  if (input_file.empty()) {
    input.reset(new ANTLRInputStream(std::cin));
  } else {
    std::ifstream f(input_file);
    if (!f.is_open()) {
      std::cerr << "Cannot open: " << input_file << std::endl;
      return 1;
    }
    input.reset(new ANTLRInputStream(f));
  }

  NECFullLexer lexer(input.get());
  CommonTokenStream tokens(&lexer);
  NECFullParser parser(&tokens);

  nec_context context;
  nec_output_file s_output;
  nec_output_flags s_output_flags;
  context.set_output(s_output, s_output_flags);

  // Use the context's internal geometry — program cards reference it.
  c_geometry* geo = context.get_geometry();
  geo->set_context(&context);

  try {
    auto* tree = parser.necFile();

    if (parser.getNumberOfSyntaxErrors() > 0) {
      std::cerr << "Parse failed with " << parser.getNumberOfSyntaxErrors()
                << " errors." << std::endl;
      return 1;
    }

    NecBuildVisitor visitor;
    visitor.nec = &context;
    visitor.geo = geo;
    tree->accept(&visitor);

    context.calc_prepare();

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
