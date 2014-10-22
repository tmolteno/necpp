## Modern nec2++ grammar for antenna geometry.

This is a parser for NEC antenna geometry files written for the ANTLR parser
generator. At the moment it is rudimentary.

The goal is to update the NEC parser to handle a more modern grammar, and
to allow useful things like multi-line comments e.t.c.

In addition, this will replace the terrible parser code that currently exists,
and will hopefully catch errors like wires accidentally intersecting as
at parse time.

Requires ANTLR

	aptitude install libantlr-dev libantlr-java

Also requires lapack libraries.

	aptitude install math-atlas-dev

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