/*
  Copyright (C) 2004,2015  Timothy C.A. Molteno
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*
  necDiff

  A program to compare two output files from NEC style 
  antenna simulations. The purpose is to support a
  testharness for the development of a C++ version of
  NEC-2.

  Author: Tim Molteno tim@physics.otago.ac.nz
*/

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#include "common.h"
#include "AntennaInput.h"
#include "RadiationInput.h"
#include "PowerBudget.h"


int main(int argc, char** argv) {

  if (argc < 2) {
    cout << "Usage: nec2diff fila, file b" << endl;
  }
  cout << "NecDiff 0.12" << endl;

  string file1(argv[1]);
  string file2(argv[2]);

  cout << "File1: " << file1 << endl;
  cout << "File2: " << file2 << endl;

  {
    AntennaInput ai1(file1);
    AntennaInput ai2(file2);

    cout << "AntennaInput Difference = " << ai1.difference(ai2) << endl;
    if (!ai1.equalto(ai2)) {
      cerr << file1 << "::" << file2;
      cerr << ". Input Parameters Different " << ai1.difference(ai2) <<  endl;
    }
  }
  
  // Compare Power Budget

  {
    PowerBudget r1(file1);
    PowerBudget r2(file2);

    cout << "PowerBudget Difference = " << r1.difference(r2) << endl;
    if (!r1.equalto(r2)) {
      cerr << file1 << "::" << file2;
      cerr << ". Power Budgets Different " << r1.difference(r2) << endl;
      // return 1;
    }
  }

  // Compare Radiation Patterns
  {
    RadiationInput r1(file1);
    RadiationInput r2(file2);

    cout << "RadiationInput Difference = " << r1.difference(r2) << endl;
    if (!r1.equalto(r2)) {
      cerr << file1 << "::" << file2;
      cerr << ". Radiation Patterns Different " << r1.difference(r2) << endl;
      // return 1;
    }
  }
  
  return 0;
}


