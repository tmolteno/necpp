/*
  Copyright (C) 2015  Timothy C.A. Molteno
  
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

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "safe_array.h"

TEST_CASE( "Resizing Array", "[safe_array]") {

    SECTION( "Capacity is bigger or equal to size" ) {
        safe_array<float> v( 5 );
        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 5 );
    }

    SECTION( "resizing bigger changes size and capacity" ) {
        safe_array<float> v( 5 );
        v.resize( 10 );

        REQUIRE( v.size() == 10 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "resizing smaller changes size but not capacity" ) {
        safe_array<float> v( 5 );
        v.resize( 0 );

        REQUIRE( v.size() == 0 );
        REQUIRE( v.capacity() >= 5 );
    }
}

TEST_CASE( "Indexing Array", "[safe_array]") {
  safe_array<float> v( 5 );

  
  SECTION( "Filling works" ) {
    v.fill(0,5,5);
    for (int i=0;i<5;i++)
      REQUIRE( v[i] == 5 );
  }

  SECTION( "Negative Indices" ) {
    REQUIRE_THROWS( v[-1] );      
  }
  SECTION( "Out of Bound Indices" ) {
    REQUIRE_THROWS( v[5] );      
  }
}


TEST_CASE( "Filling Array", "[safe_array]") {
    safe_array<float> v( 5 );

    v.fill(0,5,5);
    
    for (int i=0;i<5;i++)
       REQUIRE( v[i] == 5 );

    SECTION( "Filling from the Middle" ) {
        v.resize( 10 );
        v.fill(5,5,0);
        
        for (int i=0;i<5;i++)
          REQUIRE( v[i] == 5 );      
        for (int i=5;i<10;i++)
          REQUIRE( v[i] == 0 );      
    }
}

TEST_CASE( "Segments", "[safe_array]") {
  safe_array<int> v( 5 );
  for (int i=0;i<5;i++)
      v[i] = i;

  const safe_array<int>& f = v.segment(0,3);
  
  REQUIRE(f.size() == 4);
  
  for (int i=0;i<5;i++)
      REQUIRE( v[i] == i );

  for (int i=0;i<3;i++)
      REQUIRE( f[i] == i );

  SECTION( "modifying segment" ) {
      
      v[2] = 1;
      
      REQUIRE( f[2] == 1 );      
      REQUIRE( v[2] == 1 );      
  }
  
  SECTION( "Zero length segment" ) {
      const safe_array<int>& g = v.segment(1,1);
              
      REQUIRE( g[0] == 1 );      
  }
}


TEST_CASE( "2D", "[safe_array]") {
    safe_matrix<int> v( 5,5 );
    for (int i=0;i<5;i++)
       for (int j=0;j<5;j++)
         v(i,j) = i*j;

    REQUIRE(v.size() == 25);
    REQUIRE(v(1,1) == 1);
    REQUIRE(v(4,1) == 4);
    
    REQUIRE_THROWS(v(5,1));
    REQUIRE_THROWS(v(5,0));
}

