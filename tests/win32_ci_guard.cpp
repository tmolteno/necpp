/*
  Windows-only CI guard for the nec2++ unit-test runner.

  On MSVC, an unhandled SEH fault (e.g. access violation), abort(), an invalid
  CRT parameter, or std::terminate all default to popping a modal dialog
  (WerFault.exe / the "retry/abort/ignore" abort box). On a headless GitHub
  Actions windows-latest runner there is no desktop to dismiss it, so the test
  process blocks until the 6-hour CI ceiling — observed as an hour+ hang with
  no output after "Start 1: necpp_unit".

  This translation unit is linked only into the nec2++_tests target. A static
  initializer (runs before main()) reconfigures the CRT so every fault class
  instead calls _exit() immediately with a distinct code, which CTest reports
  as a normal (fast) failure. The production nec2++/nec2diff binaries are
  untouched, so end-user crash diagnostics are preserved.

  On non-Windows this file compiles to nothing.
*/
#if defined(_WIN32)

#include <cstdlib>
#include <exception>
#include <crtdbg.h>
#include <windows.h>

struct WinCIGuard {
    WinCIGuard() {
        // Suppress the critical-error / WER / GPF hard-error message boxes.
        SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
        // abort() must not show its "retry/abort/ignore" dialog or report to WER.
        _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
        // Invalid CRT arguments (e.g. bad printf format) must not invoke the
        // debugger. _exit() avoids running atexit/DTOR chains that could re-fault.
        _set_invalid_parameter_handler([](const wchar_t*, const wchar_t*,
                                          const wchar_t*, unsigned int, uintptr_t) {
            _exit(134);
        });
        // std::terminate (e.g. an exception escaping during stack unwinding) —
        // exit immediately rather than letting the default handler call abort(),
        // which would otherwise re-trigger the dialog path above.
        std::set_terminate([]() { _exit(133); });
    }
} static g_win_ci_guard;

#endif /* defined(_WIN32) */
