/*
  Custom ANTLR 4 error listener for NEC files.

  Reports parse errors with file name, line number, the card mnemonic
  on the offending line, and a human-readable explanation.

  Copyright (C) 2025  Timothy C.A. Molteno
  SPDX-License-Identifier: GPL-3.0-or-later
*/

#pragma once

#include <iostream>
#include <string>
#include "antlr4-runtime.h"

class NecErrorListener : public antlr4::BaseErrorListener {
public:
  explicit NecErrorListener(const std::string& filename)
    : _filename(filename) {}

  void syntaxError(antlr4::Recognizer* recognizer,
                   antlr4::Token* offendingSymbol,
                   size_t line, size_t charPositionInLine,
                   const std::string& msg,
                   std::exception_ptr e) override {
    (void)recognizer; (void)e;

    std::string card = _cardAt(line);
    std::string mnemonic = card.substr(0, std::min(card.size(), size_t(2)));

    std::cerr << _filename << ":" << line << ":" << charPositionInLine
              << ": ";

    if (!mnemonic.empty()) {
      std::cerr << mnemonic << " card: ";
    }

    if (offendingSymbol && offendingSymbol->getType() == antlr4::Token::EOF) {
      std::cerr << "unexpected end of file";
    } else if (offendingSymbol) {
      std::string extra = offendingSymbol->getText();
      // Detect "extra field" pattern — NUM token where the card
      // mnemonic matches but there are too many fields.
      std::cerr << "unexpected extra field '" << extra << "'"
                << " — card has more fields than expected by the "
                << mnemonic << " specification";
    } else {
      std::cerr << msg;
    }

    std::cerr << std::endl;
  }

  // Store file content so we can look up card mnemonics by line.
  void setSource(const std::string& source) {
    _lines.clear();
    size_t pos = 0;
    while (pos < source.size()) {
      size_t nl = source.find('\n', pos);
      if (nl == std::string::npos) {
        _lines.push_back(source.substr(pos));
        break;
      }
      _lines.push_back(source.substr(pos, nl - pos));
      pos = nl + 1;
    }
  }

private:
  std::string _filename;
  std::vector<std::string> _lines;

  std::string _cardAt(size_t line) {
    if (line == 0 || line > _lines.size()) return "";
    return _lines[line - 1];
  }
};
