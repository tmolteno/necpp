/*
 * WASM wrapper for nec2++ — exposes a minimal C API for Emscripten builds.
 *
 * Build with:
 *   emcc -std=c++17 -O2 -I src -isystem src/eigen -I build/simple \
 *        -s WASM=1 -s EXPORTED_FUNCTIONS='[\"_nec_create_context\",\"_nec_delete_context\",\"_nec_process_input\",\"_nec_get_output\",\"_nec_get_output_length\",\"_nec_free\"]' \
 *        -s EXPORTED_RUNTIME_METHODS='[\"ccall\",\"cwrap\",\"UTF8ToString\",\"lengthBytesUTF8\"]' \
 *        -s ALLOW_MEMORY_GROWTH=1 \
 *        -o nec2pp.js src/nec_wasm.cpp ...
 */

#include "nec_context.h"
#include "nec_exception.h"

#include <string>
#include <sstream>
#include <cstring>

extern "C" {

/* Opaque handle for a NEC simulation context. */
struct nec_wasm_context {
	nec_context ctx;
	std::string output_buffer;
	bool initialized = false;
	bool has_results = false;
};

nec_wasm_context* nec_create_context(void)
{
	return new nec_wasm_context();
}

void nec_delete_context(nec_wasm_context* c)
{
	delete c;
}

/*
 * Process a complete NEC input file (as a C string).
 * Returns 0 on success, or a negative error code:
 *   -1 : null context or input
 *   -2 : parse/execution error (message stored in output)
 */
int nec_process_input(nec_wasm_context* c, const char* input_text)
{
	if (!c || !input_text)
		return -1;

	try {
		std::istringstream input(input_text);
		std::ostringstream output;

		nec_output_file s_output;
		s_output.set_stream(output);
		nec_output_flags s_output_flags;

		c->ctx.set_output(s_output, s_output_flags);
		c->ctx.initialize();

		/* Parse geometry from the stream. */
		/* TODO: Wire up stream-based geometry parsing when available. */
		/* For now, this is a stub showing the API shape. */
		c->ctx.get_geometry()->parse_geometry(&c->ctx, stdin);

		c->ctx.calc_prepare();

		c->output_buffer = output.str();
		c->has_results = true;
		return 0;

	} catch (const nec_exception& e) {
		c->output_buffer = std::string("Error: ") + e.get_message();
		return -2;
	} catch (const std::exception& e) {
		c->output_buffer = std::string("Error: ") + e.what();
		return -2;
	} catch (...) {
		c->output_buffer = "Unknown error";
		return -2;
	}
}

/* Returns the output buffer (caller must NOT free). */
const char* nec_get_output(nec_wasm_context* c)
{
	if (!c)
		return "";
	return c->output_buffer.c_str();
}

/* Length of the output string, for JS interop. */
int nec_get_output_length(nec_wasm_context* c)
{
	if (!c)
		return 0;
	return static_cast<int>(c->output_buffer.size());
}

/* Free a string returned by the API (for future use). */
void nec_free(void* ptr)
{
	/* Currently unused — strings are owned by nec_wasm_context. */
	(void)ptr;
}

} /* extern "C" */
