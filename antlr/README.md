## Modern nec2++ grammar for antenna geometry.

This is a parser for NEC antenna geometry files written for the ANTLR parser
generator. At the moment it is rudimentary.

The goal is to update the NEC parser to handle a more modern grammar, and
to allow useful things like multi-line comments e.t.c.

In addition, this will replace the terrible parser code that currently exists,
and will hopefully catch errors like wires accidentally intersecting as
at parse time.

Requires ANTLR

	aptitude install libantlr-dev

Also requires lapack libraries.

	aptitude install math-atlas-dev

## Overview of new grammar

The grammar is in the file nec.g
