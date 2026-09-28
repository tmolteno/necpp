// Runtime smoke test for the Emscripten build (nec2pp.js + nec2pp.wasm).
//
// Usage: node tests/wasm_smoke.js [path/to/nec2pp.js]
//
// Feeds nec_process_input a deck that makes the geometry parser throw a
// nec_exception and checks it comes back as rc -2 with the message in the
// output buffer. If the build lacks C++ exception catching, the throw
// escapes into JS instead and this test fails.
//
// nec_process_input is still a stub that parses geometry from stdin rather
// than from its input string, so stdin is set to EOF here.
'use strict';

const path = require('path');

const jsPath = path.resolve(process.argv[2] || 'nec2pp.js');
const M = require(jsPath);

function fail(msg) {
  console.error('FAIL: ' + msg);
  process.exit(1);
}

M.stdin = () => null;
M.onAbort = (what) => fail('module aborted: ' + what);
M.onRuntimeInitialized = () => {
  const ctx = M.ccall('nec_create_context', 'number', [], []);
  if (!ctx)
    fail('nec_create_context returned null');

  let rc;
  try {
    rc = M.ccall('nec_process_input', 'number', ['number', 'string'],
                 [ctx, 'CE\nEN\n']);
  } catch (e) {
    fail('nec_process_input threw into JS (exception catching disabled?): ' + e);
  }

  const out = M.UTF8ToString(M.ccall('nec_get_output', 'number', ['number'], [ctx]));
  M.ccall('nec_delete_context', null, ['number'], [ctx]);

  if (rc !== -2)
    fail('expected rc -2, got ' + rc + ' (output: ' + JSON.stringify(out) + ')');
  if (!out.startsWith('Error: GEOMETRY DATA CARD ERROR'))
    fail('unexpected output: ' + JSON.stringify(out));

  console.log('WASM smoke test passed: ' + out);
};
