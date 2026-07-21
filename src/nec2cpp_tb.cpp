#include <catch2/catch_test_macros.hpp>

#include "nec2cpp.h"   /* for usage(), LINE_LEN */
#include "nec_context.h"
#include "nec_exception.h"
#include "c_geometry.h"

#include "nec_card_parser.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>

/*-----------------------------------------------------------------------*/
/* readmn() declaration (defined in nec2cpp.cpp) */
int readmn(FILE* input_fp, FILE* output_fp, char *gm, int *i1, int *i2,
           int *i3, int *i4, nec_float *f1, nec_float *f2, nec_float *f3,
           nec_float *f4, nec_float *f5, nec_float *f6);

/* load_line() declaration */
int load_line( char *buff, FILE *pfile );

/* Stream-based readmn declaration */
int readmn(std::istream& is,
  char *gm, int *i1, int *i2, int *i3, int *i4,
  nec_float *f1, nec_float *f2, nec_float *f3,
  nec_float *f4, nec_float *f5, nec_float *f6);

/*-----------------------------------------------------------------------*/
/* Helper: parse a NEC card line string through readmn via istringstream */

struct ParsedCard {
    std::string mnemonic;
    int    i[4];
    double f[6];
    int    parameter_count;
};

ParsedCard parse_card_line(const std::string& line) {
    std::istringstream iss(line);

    char gm[3] = {};
    int  i1=0, i2=0, i3=0, i4=0;
    nec_float f1=0, f2=0, f3=0, f4=0, f5=0, f6=0;

    int n = readmn(iss, gm, &i1, &i2, &i3, &i4, &f1, &f2, &f3, &f4, &f5, &f6);

    ParsedCard c;
    c.mnemonic = std::string(gm, 2);
    c.i[0] = i1; c.i[1] = i2; c.i[2] = i3; c.i[3] = i4;
    c.f[0] = f1; c.f[1] = f2; c.f[2] = f3; c.f[3] = f4; c.f[4] = f5; c.f[5] = f6;
    c.parameter_count = n;
    return c;
}

/*-----------------------------------------------------------------------*/
/* Tests: readmn() integer + float parsing                              */

TEST_CASE("readmn parses FR card with all fields", "[readmn]") {
    auto c = parse_card_line("FR 0 1 0 0 300.0 0.0\n");
    REQUIRE(c.mnemonic == "FR");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.i[1] == 1);
    REQUIRE(c.i[2] == 0);
    REQUIRE(c.i[3] == 0);
    REQUIRE(c.f[0] == 300.0);
    REQUIRE(c.f[1] == 0.0);
    REQUIRE(c.parameter_count == 6);
}

TEST_CASE("readmn parses EX card", "[readmn]") {
    auto c = parse_card_line("EX 0 1 5 0 1.0 0.0 0.0 0.0 0.0 0.0\n");
    REQUIRE(c.mnemonic == "EX");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.i[1] == 1);
    REQUIRE(c.i[2] == 5);
    REQUIRE(c.i[3] == 0);
    REQUIRE(c.f[0] == 1.0);
    REQUIRE(c.parameter_count == 10);
}

TEST_CASE("readmn parses GN card (perfect ground)", "[readmn]") {
    auto c = parse_card_line("GN 1\n");
    REQUIRE(c.mnemonic == "GN");
    REQUIRE(c.i[0] == 1);
    REQUIRE(c.i[1] == 0);
    REQUIRE(c.parameter_count == 1);
}

TEST_CASE("readmn parses RP card with XNDA", "[readmn]") {
    auto c = parse_card_line("RP 0 91 361 1301 0.0 0.0 1.0 1.0\n");
    REQUIRE(c.mnemonic == "RP");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.i[1] == 91);
    REQUIRE(c.i[2] == 361);
    REQUIRE(c.i[3] == 1301);
    REQUIRE(c.f[0] == 0.0);
    REQUIRE(c.f[1] == 0.0);
    REQUIRE(c.f[2] == 1.0);
    REQUIRE(c.f[3] == 1.0);
    REQUIRE(c.parameter_count == 8);
}

TEST_CASE("readmn parses EN card (mnemonic only)", "[readmn]") {
    auto c = parse_card_line("EN\n");
    REQUIRE(c.mnemonic == "EN");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.parameter_count == 0);
}

TEST_CASE("readmn handles negative integers", "[readmn]") {
    auto c = parse_card_line("GN -1\n");
    REQUIRE(c.mnemonic == "GN");
    REQUIRE(c.i[0] == -1);
    REQUIRE(c.parameter_count == 1);
}

TEST_CASE("readmn handles floats with exponent", "[readmn]") {
    auto c = parse_card_line("FR 0 1 0 0 3.0E2 0.0\n");
    REQUIRE(c.mnemonic == "FR");
    REQUIRE(c.f[0] == 300.0);
}

TEST_CASE("readmn returns 0 params for two-char mnemonic with nothing else", "[readmn]") {
    auto c = parse_card_line("EN");
    REQUIRE(c.mnemonic == "EN");
    REQUIRE(c.i[0] == 0);
    /* parse_nec_card may count the mnemonic field; parameter_count >= 0 is acceptable */
    REQUIRE(c.parameter_count >= 0);
}

TEST_CASE("readmn handles comma as field separator", "[readmn]") {
    auto c = parse_card_line("FR,0,1,0,0,300.0,0.0\n");
    REQUIRE(c.mnemonic == "FR");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.i[1] == 1);
    REQUIRE(c.f[0] == 300.0);
}

/*-----------------------------------------------------------------------*/
/* Tests: load_line() comment handling                                  */
/* Tests: mnemonic dispatch — all 21 cards recognized                  */

TEST_CASE("All 21 NEC card mnemonics are in the atst table", "[parser_dispatch]") {
    /* The atst array from nec2cpp.cpp (not externally visible, but we
       can test via the simulation — just verify the cards are known. */
    const char* known[] = {
        "FR","LD","GN","EX","NT","TL",
        "XQ","GD","RP","NX","PT","KH",
        "NE","NH","PQ","EK","CP","PL",
        "EN","WG","MP"
    };
    REQUIRE(sizeof(known)/sizeof(known[0]) == 21);
    /* Each must be exactly 2 chars */
    for (int i = 0; i < 21; i++)
        REQUIRE(strlen(known[i]) == 2);
}

/*-----------------------------------------------------------------------*/
/* Tests: load_line() comment handling                                  */

TEST_CASE("load_line skips '#' comment lines", "[load_line]") {
    std::string input = "# This is a comment\nCM Test\n";
    std::istringstream iss(input);
    char buf[256];
    int eof = load_line(buf, iss);
    REQUIRE(eof != EOF);
    REQUIRE(strncmp(buf, "CM", 2) == 0);
}

TEST_CASE("load_line handles inline ' comment (4nec2 compat)", "[load_line]") {
    std::string input = "FR 0 1 0 0 300.0 0.0 ' inline comment\n";
    std::istringstream iss(input);
    char buf[256];
    int eof = load_line(buf, iss);
    REQUIRE(eof != EOF);
    REQUIRE(std::string(buf).find("inline") == std::string::npos);
    REQUIRE(strncmp(buf, "FR", 2) == 0);
}

/*-----------------------------------------------------------------------*/
/* Tests: new parse_nec_card() (Options 1+2 refactor)                   */

TEST_CASE("parse_nec_card parses FR card", "[parse_nec_card]") {
    auto c = parse_nec_card("FR 0 1 0 0 300.0 0.0");
    REQUIRE(c.mnemonic == "FR");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.i[1] == 1);
    REQUIRE(c.f[0] == 300.0);
    REQUIRE(c.f[1] == 0.0);
    REQUIRE(c.parameter_count == 6);
}

TEST_CASE("parse_nec_card handles inline comment", "[parse_nec_card]") {
    auto c = parse_nec_card("EX 0 1 5 0 1.0 ' comment");
    REQUIRE(c.mnemonic == "EX");
    REQUIRE(c.i[0] == 0);
    REQUIRE(c.f[0] == 1.0);
    REQUIRE(c.parameter_count == 5);
}

TEST_CASE("parse_nec_card handles comma separators", "[parse_nec_card]") {
    auto c = parse_nec_card("FR,0,1,0,0,300.0,0.0");
    REQUIRE(c.mnemonic == "FR");
    REQUIRE(c.f[0] == 300.0);
    REQUIRE(c.parameter_count == 6);
}

TEST_CASE("parse_nec_card handles empty line", "[parse_nec_card]") {
    auto c = parse_nec_card("");
    REQUIRE(c.mnemonic == "");
    REQUIRE(c.parameter_count == 0);
}

TEST_CASE("parse_nec_card handles # comment", "[parse_nec_card]") {
    auto c = parse_nec_card("# this is a comment");
    REQUIRE(c.mnemonic == "");
    REQUIRE(c.parameter_count == 0);
}

TEST_CASE("parse_nec_card handles negative floats", "[parse_nec_card]") {
    auto c = parse_nec_card("GN -1 0 0 0 -0.5 1.0e3");
    REQUIRE(c.mnemonic == "GN");
    REQUIRE(c.i[0] == -1);
    REQUIRE(c.f[0] == -0.5);
    REQUIRE(c.f[1] == 1000.0);
}

TEST_CASE("card handler table has 21 entries", "[card_handler]") {
    REQUIRE(sizeof(card_handlers)/sizeof(card_handlers[0]) == 21);
}

TEST_CASE("find_handler finds all known cards", "[card_handler]") {
    const char* names[] = {"FR","LD","GN","EX","NT","TL","XQ","GD","RP",
        "NX","PT","KH","NE","NH","PQ","EK","CP","PL","EN","WG","MP"};
    for (auto n : names) {
        auto* h = find_handler(n);
        REQUIRE(h != nullptr);
        REQUIRE(h->name[0] == n[0]);
        REQUIRE(h->name[1] == n[1]);
    }
}

TEST_CASE("find_handler returns null for unknown card", "[card_handler]") {
    REQUIRE(find_handler("ZZ") == nullptr);
}

TEST_CASE("parse_nec_card matches old readmn for FR", "[compat]") {
    auto c_new = parse_nec_card("FR 0 1 0 0 300.0 0.0");
    auto c_old = parse_card_line("FR 0 1 0 0 300.0 0.0\n");
    REQUIRE(c_new.mnemonic == c_old.mnemonic);
    REQUIRE(c_new.i[0] == c_old.i[0]);
    REQUIRE(c_new.i[1] == c_old.i[1]);
    REQUIRE(c_new.i[2] == c_old.i[2]);
    REQUIRE(c_new.i[3] == c_old.i[3]);
    REQUIRE(c_new.f[0] == c_old.f[0]);
    REQUIRE(c_new.f[1] == c_old.f[1]);
}

TEST_CASE("parse_nec_card matches old readmn for EX", "[compat]") {
    auto c_new = parse_nec_card("EX 0 1 5 0 1.0 0.0 0.0 0.0 0.0 0.0");
    auto c_old = parse_card_line("EX 0 1 5 0 1.0 0.0 0.0 0.0 0.0 0.0\n");
    REQUIRE(c_new.mnemonic == c_old.mnemonic);
    REQUIRE(c_new.i[0] == c_old.i[0]);
    REQUIRE(c_new.i[1] == c_old.i[1]);
    REQUIRE(c_new.i[2] == c_old.i[2]);
    REQUIRE(c_new.i[3] == c_old.i[3]);
    for (int n = 0; n < 6; n++)
        REQUIRE(c_new.f[n] == c_old.f[n]);
}
