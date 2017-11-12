/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "HistogramBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Pbc.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDCALC SPHERICAL_KDE
/*
Accumulate the average probability density on a spherical grid.

\par Examples

*/
//+ENDPLUMEDOC

class SphericalKDE : public HistogramBase {
private:
  double hh;
  std::vector<double> center;
  double von_misses_norm;
  double von_misses_concentration;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit SphericalKDE(const ActionOptions&ao);
  void setupNeighborsVector();
  void getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& out_nbin, 
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void buildSingleKernel( std::vector<unsigned>& tflags, const double& height, std::vector<double>& args );
  double calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const ;
  void addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const ;
};

PLUMED_REGISTER_ACTION(SphericalKDE,"SPHERICAL_KDE")
PLUMED_REGISTER_SHORTCUT(SphericalKDE,"SPHERICAL_KDE")

void SphericalKDE::shortcutKeywords( Keywords& keys ){
  HistogramBase::shortcutKeywords( keys );
}

void SphericalKDE::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions ) {
  HistogramBase::resolveNormalizationShortcut( lab, words, keys, actions );
}

void SphericalKDE::registerKeywords( Keywords& keys ) {
  HistogramBase::registerKeywords( keys );
  keys.add("compulsory","GRID_BIN","the number of points on the fibonacci sphere");
  keys.add("compulsory","CONCENTRATION","the concentration parameter for Von Mises-Fisher distributions");
}

SphericalKDE::SphericalKDE(const ActionOptions&ao):
  Action(ao),
  HistogramBase(ao),
  hh(0),
  center(getNumberOfDerivatives(),0)
{
  if( getNumberOfDerivatives()!=3 ) error("should have three coordinates in input to this action"); 

  unsigned nbins; parse("GRID_BIN",nbins);
  log.printf("  setting number of bins to %d \n", nbins );
  parse("CONCENTRATION",von_misses_concentration);
  von_misses_norm = von_misses_concentration / ( 4*pi*sinh( von_misses_concentration ) );
  log.printf("  setting concentration parameter to %f \n", von_misses_concentration );

  // Create a value
  std::vector<bool> ipbc( getNumberOfDerivatives(), false ); 
  double fib_cutoff = std::log( epsilon / von_misses_norm ) / von_misses_concentration;
  gridobject.setup( "fibonacci", ipbc, nbins, fib_cutoff ); checkRead();
  
  // Setup the grid 
  std::vector<unsigned> shape(1); shape[0]=nbins; 
  addValueWithDerivatives( shape ); setupNeighborsVector(); 
}

void SphericalKDE::setupNeighborsVector() { }

void SphericalKDE::getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min,
                                         std::vector<std::string>& max, std::vector<unsigned>& out_nbin, 
                                         std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  error("cannot print spherical grids in this way");
} 

void SphericalKDE::buildSingleKernel( std::vector<unsigned>& tflags, const double& height, std::vector<double>& args ) {
  unsigned num_neigh; std::vector<unsigned> neighbors, nneigh;
  hh=height; for(unsigned i=0;i<args.size();++i) center[i]=args[i];
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors ); 
  for(unsigned i=0;i<num_neigh;++i) tflags[ neighbors[i] ] = 1;
}

double SphericalKDE::calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const {
  double dot=0; for(unsigned i=0; i<der.size(); ++i){ dot += args[i]*center[i]; } 
  double newval = hh*von_misses_norm*exp( von_misses_concentration*dot );
  for(unsigned i=0; i<der.size(); ++i) der[i] = von_misses_concentration*newval*args[i];
  return newval;
}

void SphericalKDE::addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const {
  std::vector<double> gpoint( args.size() ), der( args.size() );
  unsigned num_neigh; std::vector<unsigned> neighbors, nneigh; 
  gridobject.getNeighbors( args, nneigh, num_neigh, neighbors );
  for(unsigned i=0;i<num_neigh;++i) {
      gridobject.getGridPointCoordinates( neighbors[i], gpoint );
      double dot=0; for(unsigned j=0; j<der.size(); ++j) dot += args[j]*gpoint[j];
      double newval = height*von_misses_norm*exp( von_misses_concentration*dot );
      buffer[ bufstart + neighbors[i]*(1+der.size()) ] += newval;
      for(unsigned j=0; j<der.size(); ++j) buffer[ bufstart + neighbors[i]*(1+der.size()) + 1 + j ] += von_misses_concentration*newval*gpoint[j];
  } 
}

}
}
