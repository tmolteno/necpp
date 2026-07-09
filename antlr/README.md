# ANTLR 4 Parser for NEC-2

## Comparison with NEC-2 Standard (nec2prt3.pdf)

The ANTLR grammar uses free-form whitespace-separated fields (`NUM?`)
rather than the fixed-width FORTRAN columns of the original NEC-2 spec.
This is consistent with the production parser's `parse_nec_card()` in
`src/nec_card_parser.h` — both treat blank fields as zero and accept
space or comma separators.

| NEC-2 Feature | Grammar Handling |
|---|---|
| Fixed-width fields (I2, I3, F10.5, etc.) | Free-form `NUM?` tokens |
| Blank fields → 0 | Optional tokens → `dval(nullptr) = 0.0` in visitor |
| CM/CE comment block before geometry | `CMT` lexer rule skips entire lines |
| GW rad=0 signals mandatory GC taper | `_gw` state in `NecBuildVisitor` |
| GE ends geometry, enters program section | `geCard` in `programCard` dispatch, calls `geometry_complete` |
| SP/SM may chain multiple SC cards | Single `scCard?` — multi-patch chaining not yet implemented |
| GF (NGF) card | Intentionally excluded (unsupported by necpp) |
| XT debug card | Not included (production parser handles it separately) |

## Building

```bash
make        # generate, compile
make test   # run test suite
make clean  # remove generated files
```

Requires Docker (no system ANTLR or C++ dependencies needed).
