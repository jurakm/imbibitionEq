/*
 * test_functions.cc
 *
 *  Created on: 24 Apr 2015
 *      Author: jurak
 */
#include "../config.h"

#include <fstream>
#include <iostream>

#include "imbibitionFunctions.hh"
#include "vangenuchten.hh"

int main()
{
  double muw = 1e-3;
  double mun = 0.8e-3;

  double PentryVG = 2.0;  // scaled
  double nVG = 2.0;
  Dumux::VanGenuchten::Params params(1.0/PentryVG, nVG);
  RealImbibitionFunctions<Dumux::VanGenuchten> vg_functions(params, muw, mun);

  double a = 1.0;
  ArtifImbibitionFunctions artif_functions(a);

//  plot("ab_vg.txt", vg_functions, 5);
  plot("ab_vg.txt", vg_functions, vg_functions.NofPts());
//  std::cout << artif_functions.NofPts() << std::endl;
//  plot("ab_artif.txt", artif_functions, 5);
  plot("ab_artif.txt", artif_functions, artif_functions.NofPts());

  return 0;
}
