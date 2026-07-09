/*
  Test driver for the ANTLR 4 NEC geometry parser.

  Reads geometry from stdin or a file, runs the ANTLR 4 parser,
  and uses GeometryBuildVisitor to call the c_geometry API.

  Build (inside Docker):
      antlr4 -Dlanguage=Cpp -o generated -visitor NECGeometry.g4
      g++ -std=c++17 -I generated -I ../src -isystem ../src/eigen \
          -I ../build/simple -I /usr/include/antlr4-runtime \
          generated/NECGeometryLexer.cpp generated/NECGeometryParser.cpp \
          generated/NECGeometryBaseVisitor.cpp \
          geometry_test.cpp \
          ../src/c_geometry.cpp ../src/c_ggrid.cpp ../src/libNEC.cpp \
          ../src/nec_exception.cpp ../src/nec_context.cpp ../src/nec_ground.cpp \
          ../src/nec_output.cpp ../src/nec_radiation_pattern.cpp \
          ../src/nec_results.cpp ../src/matrix_algebra.cpp ../src/misc.cpp \
          ../src/electromag.cpp ../src/c_evlcom.cpp ../src/c_plot_card.cpp \
          ../src/nec_structure_currents.cpp ../src/XGetopt.cpp \
          -lantlr4-runtime -lm -lstdc++ -o geometry_parse
*/

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include "antlr4-runtime.h"
#include "NECGeometryLexer.h"
#include "NECGeometryParser.h"
#include "geometry_visitor.h"

#include "nec_context.h"
#include "c_geometry.h"
#include "nec_exception.h"

using namespace antlr4;

int main(int argc, char* argv[]) {
  std::string input_file;
  if (argc > 1) input_file = argv[1];

  // Set up input stream
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

  // Create lexer and parser
  NECGeometryLexer lexer(input.get());
  CommonTokenStream tokens(&lexer);
  NECGeometryParser parser(&tokens);

  // Set up geometry
  nec_context context;
  c_geometry geo;
  geo.set_context(&context);

  nec_output_file s_output;
  nec_output_flags s_output_flags;
  context.set_output(s_output, s_output_flags);

  // Parse
  try {
    auto* tree = parser.geometryFile();

    if (parser.getNumberOfSyntaxErrors() > 0) {
      std::cerr << "Parse failed with " << parser.getNumberOfSyntaxErrors()
                << " errors." << std::endl;
      return 1;
    }

    // Walk with visitor
    GeometryBuildVisitor visitor;
    visitor.nec = &context;
    visitor.geo = &geo;
    tree->accept(&visitor);

    std::cout << "Parsed successfully. Wires: " << visitor.nwire
              << ", Segments: " << geo.n_segments << std::endl;
    return 0;

  } catch (const nec_exception& e) {
    std::cerr << "Geometry error: " << e.get_message() << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
