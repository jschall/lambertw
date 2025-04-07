/*
  Command-line test of accuracy of Lambert W function implementations

  Copyright (C) 2015 Darko Veberic, darko.veberic@ijs.si

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <LambertW.h>
#include <FukushimaLambertW.h>
#include "lambw_grok.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;


inline
bool
CloseTo(const double a, const double b, const double eps)
{
  return fabs(a - b) < eps;
}


int
main()
{
  const double eps = 3e-5;
  const double step1 = 0.00001;
  const double max2 = 30;
  const double step2 = 0.00001;

  int count = 0;
  for (int branch = -1; branch <= 0; ++branch) {
    for (double x = -1/M_E; x < 0; x += step1) {
      count++;
      const float w = branch==0?lambert_w0(x):lambert_wm1(x);
      const float fw = Fukushima::LambertW(branch, x);
      if (!CloseTo(w, fw, eps*5)) {
        cout << count << " fail " << " f" << branch << '(' << x << ") = " << (w - fw) << endl;
//         return EXIT_FAILURE;
      }
    }
      cout << "success branch " << branch << " count " << count << endl;
      count = 0;
  }

  for (double x = 0; x < max2; x += step2) {
    count++;
    const double w = lambert_w0(x);
    const double fw = Fukushima::LambertW(0, x);
    if (!CloseTo(w, fw, eps)) {
      cout << count << " fail " << "f0(" << x << ") = " << (w - fw) << endl;
    }
  }
  cout << "success branch 0 count " << count << endl;
  count = 0;

  return EXIT_SUCCESS;
}
