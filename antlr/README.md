## Modern nec2++ grammar for antenna geometry.

**DEPRECATED — this is experimental, abandoned code.** The ANTLR-based parser
was a proof-of-concept and is far from complete. The production parser lives in
`src/nec2cpp.cpp` and uses a simple `load_line()` + card dispatch approach.

If revived, it would need:
- ANTLR 4 (ANTLR 2.x C++ target is discontinued)
- Full card coverage (SP, GM, GR, NE, NH, KH, etc.)
- Integration with the Eigen-backed math layer

---

Original README follows:

This is a parser for NEC antenna geometry files written for the ANTLR parser
generator. At the moment it is rudimentary.

The goal is to update the NEC parser to handle a more modern grammar, and
to allow useful things like multi-line comments e.t.c.

In addition, this will replace the terrible parser code that currently exists,
and will hopefully catch errors like wires accidentally intersecting as
at parse time.

Requires ANTLR (legacy v2):

    aptitude install libantlr-dev libantlr-java

### Testing

Just build in the current directory

    make 
    make test_clean
    make test_all

## New parser for NEC cards

The grammar is in the file nec.g. This is an ANTLR grammar for the existing NEC2 card deck for describing
antennas.


## New Language

The new language for nec++ will be more explicit and easy to read. It will also allow comments anywhere 
inside the code. In addition it will not depend on whitespace to skip over missing parameters, and therefore
will be far more robust.

    /* This is a comment */
    geometry {
      w0 = wire(start=[1,2,3], end=[2,3,4], r=0.01, n=5);
      w1 = wire(start=w0.end,  end=[2,3,5], r=0.01, n=5);
      arc(origin=[0,1,2], arc_radius=1.0, r=0.01 );
      helix();
      patch();
      scale(3.0);
    }

    ground {
      type=PERFECT;
    }

    excitation {
      type=VOLTAGE;
      segment=w0.0; /* Segment indexing starts at 0 */
      freq=range(start=1.575GHz, end=1.675GHz, n=5);
      extended_thin_wire_kernel = false; /* default */
    }
        
    radiation_pattern {
      mode = "normal";
      theta = range(start=0.0, end=180, n=30);
      phi = range(start=0.0, end=90, n=30);
    }

    execute();

Expression Syntax:

    { id = } TYPE([param=expression]*) ;