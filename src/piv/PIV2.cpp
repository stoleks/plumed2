/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2018 of Pipolo Silvio, Alexandre Jedrecy and Fabio Pietrucci.

The piv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The piv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

#include <string>
#include <cmath>
#include <iostream>

namespace PLMD
{
namespace piv
{

//+PLUMEDOC COLVAR PIV2
/*
Calculates the PIV2-distance: the squared Cartesian distance between the PIV2 \cite gallet2013structural,pipolo2017navigating
associated to the configuration of the system during the dynamics and a reference configuration provided
as input (PDB file format).
PIV2 can be used together with \ref FUNCPATHMSD to define a path in the PIV2 space.
\par Examples

The following example calculates PIV2-distances from three reference configurations in Ref1.pdb, Ref2.pdb and Ref3.pdb
and prints the results in a file named colvar.
Three atoms (PIV2ATOMS=3) with names (pdb file) A B and C are used to construct the PIV2 and all PIV2 blocks (AA, BB, CC, AB, AC, BC) are considered.
SFACTOR is a scaling factor that multiplies the contribution to the PIV2-distance given by the single PIV2 block.
NLIST sets the use of neighbor lists for calculating atom-atom distances.
The SWITCH keyword specifies the perameters of the switching function that transforms atom-atom distances.
SORT=1 meand that the PIV2 block elements are sorted (SORT=0 no sorting.)
Values for SORT, SFACTOR and Neighborlist parameters have to be specified for each block.
The order is the following: AA,BB,CC,AB,AC,BC. If ONLYDIRECT (ONLYCROSS) is used the order is AA,BB,CC (AB,AC,BC).
The sorting operation within each PIV2 block is performed using the counting sort algorithm, PRECISION specifies the size of the counting array.
\plumedfile
PIV2 ...
LABEL=Pivd1
PRECISION=1000
NLIST
REF_FILE=Ref1.pdb
PIV2ATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV2
PIV2 ...
LABEL=Pivd2
PRECISION=1000
NLIST
REF_FILE=Ref2.pdb
PIV2ATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV2
PIV2 ...
LABEL=Pivd3
PRECISION=1000
NLIST
REF_FILE=Ref3.pdb
PIV2ATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV2

PRINT ARG=Pivd1,Pivd2,Pivd3 FILE=colvar
\endplumedfile

WARNING:
Both the "CRYST" and "ATOM" lines of the PDB files must conform precisely to the official pdb format, including the width of each alphanumerical field:

\verbatim
CRYST1   31.028   36.957   23.143  89.93  92.31  89.99 P 1           1
ATOM      1  OW1 wate    1      15.630  19.750   1.520  1.00  0.00
\endverbatim

In each pdb frame, atoms must be numbered in the same order and with the same element symbol as in the input of the MD program.

The following example calculates the PIV2-distances from two reference configurations Ref1.pdb and Ref2.pdb
and uses PIV2-distances to define a Path Collective Variable (\ref FUNCPATHMSD) with only two references (Ref1.pdb and Ref2.pdb).
With the VOLUME keyword one scales the atom-atom distances by the cubic root of the ratio between the specified value and the box volume of the initial step of the trajectory file.

\plumedfile
PIV2 ...
LABEL=c1
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref1.pdb
PIV2ATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.5 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV2
PIV2 ...
LABEL=c2
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref2.pdb
PIV2ATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV2

p1: FUNCPATHMSD ARG=c1,c2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=c1,c2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIV2 please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class PIV2 : public Colvar
{
public:
  static void registerKeywords( Keywords& keys );
  PIV2(const ActionOptions&);
  ~PIV2();
  // active methods:
  virtual void calculate();
  virtual void computeReferencePIV();
  void checkFieldsAllowed() {}
private:
  Vector distanceAB(Vector A, Vector B);
private:
  ForwardDecl<Stopwatch> stopwatch_fwd;
  /// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch=*stopwatch_fwd;
  bool m_pbc, m_serial, m_timer, m_referenceIsInitialized, m_neighborIsInitialized;
  bool m_doScaleVolume, m_cross, m_direct, m_doNeighbor, m_test, m_computeDerivatives, m_centerOfMass;
  unsigned m_nPrecision, m_nAtomTypes, m_numberBlocks, m_neighborlistSize, m_numberReferences;
  int m_updatePIV2;
  double m_volumeFactor, m_volume0, m_PIV2distance;
  Tensor m_virial;
  NeighborList *m_neighborListAll;
  std::vector<bool> m_doSort;
  std::vector<double> m_blockScaling, m_r00;
  std::vector<double> m_neighborListSkin;
  std::vector<double> m_massFactor;
  std::vector<Vector> m_positionCOM;
  std::vector<Vector> m_deriv;
  std::vector<std::string> sw;
  std::vector<NeighborList*> m_neighborList;
  std::vector<NeighborList*> m_neighborlistCOM;
  std::vector<SwitchingFunction> m_switchingFunctions;
  std::vector<std::vector<double> > m_referencePIV2;
};

PLUMED_REGISTER_ACTION(PIV2,"PIV2")

void PIV2::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered", "SWITCH", "The switching functions parameter."
           "You should specify a Switching function for all PIV2 blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("compulsory", "PRECISION", "the precision for approximating reals with integers in sorting.");
  keys.add("compulsory", "REF_FILE", "PDB file name that contains the i-th reference structure.");
  keys.add("compulsory", "PIV2ATOMS", "Number of atoms to use for PIV2.");
  keys.add("compulsory", "SORT","Whether to sort or not the PIV2 block.");
  keys.add("compulsory", "ATOMTYPES", "The atomtypes to use for PIV2.");
  keys.add("optional", "SFACTOR", "Scale the PIV2-distance by such block-specific factor");
  keys.add("optional", "VOLUME", "Scale atom-atom distances by the cubic root of the cell volume. The input volume is used to scale the R_0 value of the switching function. ");
  keys.add("optional", "UPDATEPIV2", "Frequency (timesteps) at which the PIV2 is updated.");
  keys.add("optional", "NL_CUTOFF", "Neighbour lists cutoff.");
  keys.add("optional", "NL_STRIDE", "Update neighbour lists every NL_STRIDE steps.");
  keys.add("optional", "NL_SKIN", "The maximum atom displacement tolerated for the neighbor lists update.");
  keys.addFlag("TEST", false, "Print the actual and reference PIV2 and exit");
  keys.addFlag("COM", false, "Use centers of mass of groups of atoms instead of atoms as secified in the Pdb file");
  keys.addFlag("ONLYCROSS", false, "Use only cross-terms (A-B, A-C, B-C, ...) in PIV2");
  keys.addFlag("ONLYDIRECT", false, "Use only direct-terms (A-A, B-B, C-C, ...) in PIV2");
  keys.addFlag("DERIVATIVES", false, "Activate the calculation of the PIV2 for every class (needed for numerical derivatives).");
  keys.addFlag("NLIST", false, "Use a neighbour list for distance calculations.");
  keys.addFlag("SERIAL", false, "Perform the calculation in serial - for debug purpose");
  keys.addFlag("TIMER", false, "Permorm timing analysis on heavy loops.");
  keys.reset_style("SWITCH", "compulsory");
}

//---------------------------------------------------------
// CONSTRUCTOR
//---------------------------------------------------------

PIV2::PIV2(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  m_pbc(true),
  m_serial(false),
  m_timer(false),
  m_referenceIsInitialized(false),
  m_neighborIsInitialized(false),
  m_doScaleVolume(false),
  m_cross(true),
  m_direct(true),
  m_doNeighbor(false),
  m_test(false),
  m_computeDerivatives(false),
  m_centerOfMass(false),
  m_nPrecision(1000),
  m_nAtomTypes(1),
  m_neighborlistSize(1),
  m_numberReferences(1),
  m_updatePIV2(1),
  m_volumeFactor(1.),
  m_volume0(0.),
  m_PIV2distance(0.) 
{
  log << "Starting PIV2 Constructor\n";
  bool onlyCross = false, onlyDirect = false;
  std::string referenceFile;

  parse("VOLUME", m_volume0);
  parse("PIV2ATOMS", m_nAtomTypes);
  parse("PRECISION", m_nPrecision);
  //parse("REF_FILE", referenceFile);
  parse("NUMBER_REF_FILE", m_numberReferences);

  parseFlag("TEST", m_test);
  parseFlag("NOPBC", m_pbc);
  parseFlag("TIMER", m_timer);
  parseFlag("SERIAL", m_serial);
  parseFlag("NLIST", m_doNeighbor);
  parseFlag("ONLYCROSS", onlyCross);
  parseFlag("ONLYDIRECT", onlyDirect);
  parseFlag("COM", m_centerOfMass);
  parseFlag("DERIVATIVES", m_computeDerivatives);

  std::vector<std::string> atomTypes(m_nAtomTypes);
  parseVector("ATOMTYPES", atomTypes);

  if (keywords.exists("UPDATEPIV2")){
    parse("UPDATEPIV2", m_updatePIV2);
  }
 
  // Precision on the real-to-integer transformation for the sorting
  if (m_nPrecision < 2) { 
    error("Precision must be => 2");
  }
  // PBC
  if (m_pbc) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }
  // Serial or parallel
  if (m_serial) {
    log << "Serial PIV2 construction\n";
  } else     {
    log << "Parallel PIV2 construction\n";
  }
  // Derivatives
  if (m_computeDerivatives) {
    log << "Computing Derivatives\n";
  }
  // Timing
  if(m_timer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }
  // Center of Mass
  if(m_centerOfMass){
    log << "Building PIV2 using COMs\n";
  }
  // Volume Scaling
  if (m_volume0 > 0) {
    m_doScaleVolume = true;
  }
  // PIV2 direct and cross blocks
  if (onlyCross && onlyDirect) {
    error("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }

  m_numberBlocks = 0;
  if(onlyCross) {
    m_direct = false;
    log << "Using only CROSS-PIV2 blocks\n";
  } 
  if(onlyDirect) {
    m_cross = false;
    log << "Using only DIRECT-PIV2 blocks\n";
  }
  if (m_cross) {
    m_numberBlocks += unsigned(double(m_nAtomTypes * (m_nAtomTypes - 1)) / 2.);
  }
  if (m_direct) {
    m_numberBlocks += unsigned(m_nAtomTypes);
  }
 
  // resizing all class vector according to m_numberBlocks
  m_referencePIV2.resize(m_numberBlocks);
  m_doSort.resize(m_numberBlocks);
  m_blockScaling.resize(m_numberBlocks);
  m_r00.resize(m_numberBlocks);
  sw.resize(m_numberBlocks); 
  m_switchingFunctions.resize(m_numberBlocks);
  m_neighborList.resize(m_numberBlocks);
  m_neighborListSkin.resize(m_numberBlocks);

  // Sorting
  std::vector<int> yesNoSort(m_numberBlocks);
  parseVector("SORT", yesNoSort);
  for (unsigned i = 0; i < m_numberBlocks; i++) {
    if(yesNoSort[i] == 0 || m_computeDerivatives) {
      m_doSort[i] = false;
    } else {
      m_doSort[i] = true;
    }
  }

  // PIV2 scaled option
  for(unsigned j = 0; j < m_numberBlocks; j++) {
    m_blockScaling[j] = 1.;
  }
  if(keywords.exists("SFACTOR")) {
    parseVector("SFACTOR", m_blockScaling);
  }
  
  // read parameters and set-up switching functions here only if computing derivatives
  for (unsigned i = 0; i < m_numberBlocks; i++) {
    if( !parseNumbered("SWITCH", i+1, sw[i]) ) break;
  }
  if (m_computeDerivatives) {
    log << "Switching Function Parameters \n";
    m_switchingFunctions.resize(m_numberBlocks);
    std::string errors;
    for (unsigned i = 0; i < m_numberBlocks; i++) {
      std::string num;
      Tools::convert(i + 1, num);
      m_switchingFunctions[i].set(sw[i], errors);
      if (errors.length() != 0){
        error("problem reading switch" + num + " keyword : " + errors );
      }
      m_r00[i] = m_switchingFunctions[i].get_r0();
      log << "  Swf: " << j << "  r0=" << (m_switchingFunctions[i].description()).c_str() << "\n";
    }
  }

  // Reference PDB file
  
  std::vector<std::string> referenceFiles;
  for (unsigned i = 0; i < m_numberBlocks; i++) {
    if ( !parseNumbered("REF_FILE", i+1, referenceFiles[i]) ) break;
  }
  PDB myPDB;
  FILE* fp = fopen(referenceFile.c_str(), "r");
  if (fp != NULL) {
    log << "Opening PDB file with reference frame: " << referenceFile.c_str() << "\n";
    myPDB.readFromFilepointer(fp, plumed.getAtoms().usingNaturalUnits(), 
                              0.1 / atoms.getUnits().getLength());
    fclose (fp);
  } else {
    error("Error in reference PDB file");
  }

  // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
  std::vector<std::vector<AtomNumber> > pointList(m_nAtomTypes);
  std::vector<std::vector<AtomNumber> > comAtoms(1);
  m_neighborlistSize = myPDB.getAtomNumbers().size();
  unsigned residueNumber = 0;
  // residueIndex (array size: number of residues) contains the contains the residue number
  // this is because the residue numbers may not alwyas be ordered from 1 to residueNumber
  std::vector<unsigned> residueIndex;
  // Build residueIndex
  for (unsigned i = 0; i < myPDB.getAtomNumbers().size(); i++) {
    unsigned residueNumber = myPDB.getResidueNumber(myPDB.getAtomNumbers()[i]);
    if (std::find_if(residueIndex.begin(), residueIndex.end(),
                     [residueNumber] (const unsigned p) 
                                     { return p == residueNumber; }) 
        == residueIndex.end()) residueIndex.push_back(residueNumber);
  }
  residueNumber = residueIndex.size();

  // Pind0 is the atom/COM used in Nlists (for COM Pind0 is the first atom in the pdb belonging to that COM)
  unsigned Pind0size;
  if(m_centerOfMass) {
    Pind0size = residueNumber;
  } else {
    Pind0size = m_neighborlistSize;
  }
  std::vector<unsigned> Pind0(Pind0size);
  // If COM resize important arrays
  comAtoms.resize(m_neighborlistSize);
  if(m_centerOfMass) {
    m_neighborlistCOM.resize(m_neighborlistSize);
    m_positionCOM.resize(m_neighborlistSize);
    m_massFactor.resize(m_neighborlistSize, 0.);
  }
  log << "Total COM/Atoms: " << m_nAtomTypes * residueNumber << " \n";
  // Build lists of Atoms/COMs for NLists
  //   comAtoms filled also for non_COM calculation for analysis purposes
  for (unsigned j = 0; j < m_nAtomTypes; j++) {
    unsigned oind;
    std::fill(Pind0.begin(), Pind0.end(), 0);
    for (unsigned i = 0; i < myPDB.getAtomNumbers().size(); i++) {
      // Residue/Atom AtomNumber: used to build NL for COMS/Atoms pairs.
      AtomNumber anum = myPDB.getAtomNumbers()[i];
      // Index associated to residue/atom: used to separate COM-lists
      unsigned residueNumber = myPDB.getResidueNumber(anum);
      unsigned aind = anum.index();
      // This builds lists for NL
      std::string Pname;
      unsigned Pind;
      if(m_centerOfMass) {
        Pname = myPDB.getAtomName(anum);
        for(unsigned l = 0; l < residueNumber; l++) {
          if(residueNumber == residueIndex[l]) {
            Pind = l;
          }
        }
      } else {
        Pname = myPDB.getAtomName(anum);
        Pind = aind;
      }
      if(Pname == atomTypes[j]) {
        if(Pind0[Pind] == 0) {
          // adding the atomnumber to the atom/COM list for pairs
          pointList[j].push_back(anum);
          Pind0[Pind] = aind+1;
          oind = Pind;
        }
        // adding the atomnumber to list of atoms for every COM/Atoms
        comAtoms[Pind0[Pind]-1].push_back(anum);
      }
    }
    // Output Lists
    log << "  Groups of type  " << j << ": " << pointList[j].size() << " \n";
    std::string gname;
    unsigned gsize;
    if(m_centerOfMass) {
      gname = myPDB.getResidueName(comAtoms[Pind0[oind] - 1][0]);
      gsize = comAtoms[Pind0[oind] - 1].size();
    } else {
      gname = myPDB.getAtomName(comAtoms[Pind0[oind] - 1][0]);
      gsize = 1;
    }
    log.printf("    %6s %3s %13s %10i %6s\n", "type  ", gname.c_str(),"   containing ",gsize," atoms");
  }

  // This is to build the list with all the atoms
  std::vector<AtomNumber> listall;
  for (unsigned i = 0; i < myPDB.getAtomNumbers().size(); i++) {
    listall.push_back(myPDB.getAtomNumbers()[i]);
  }

  if(m_doNeighbor) {
    std::vector<double> m_neighborListCut(m_numberBlocks,0.);
    std::vector<int> m_neighborListStride(m_numberBlocks,0);
    parseVector("NL_CUTOFF", m_neighborListCut);
    parseVector("NL_STRIDE", m_neighborListStride);
    parseVector("NL_SKIN", m_neighborListSkin);
    for (unsigned j = 0; j < m_numberBlocks; j++) {
      if (m_neighborListCut[j] <= 0.0) error("NL_CUTOFF should be explicitly specified and positive");
      if (m_neighborListStride[j] <= 0)    error("NL_STRIDE should be explicitly specified and positive");
      if (m_neighborListSkin[j] <= 0.) error("NL_SKIN should be explicitly specified and positive");
      m_neighborListCut[j] = m_neighborListCut[j] + m_neighborListSkin[j];
    }
    log << "Creating Neighbor Lists \n";
    // WARNING: is m_neighborListCut meaningful here?
    m_neighborListAll= new NeighborList(listall, m_pbc, getPbc(), m_neighborListCut[0], m_neighborListStride[0]);
    log << "neighborlist size = " << m_neighborList.size() << "\n";
    if(m_centerOfMass) {
      log << "center of mass \n";
      //Build lists of Atoms for every COM
      for (unsigned i = 0; i < m_positionCOM.size(); i++) {
        // WARNING: is m_neighborListCut meaningful here?
        m_neighborlistCOM[i] = new NeighborList(comAtoms[i], m_pbc, getPbc(), m_neighborListCut[0], m_neighborListStride[0]);
      }
    }
    unsigned ncnt = 0;
    // Direct blocks AA, BB, CC, ...
    if(m_direct) {
      log << m_neighborList.size() << " " << m_nAtomTypes << "\n";
      for (unsigned j = 0; j < m_nAtomTypes; j++) {
        m_neighborList[ncnt] = new NeighborList(pointList[j], m_pbc, getPbc(), m_neighborListCut[j], m_neighborListStride[j]);
        ncnt += 1;
      }
    }
    // Cross blocks AB, AC, BC, ...
    if(m_cross) {
      log << m_neighborList.size() << " " << m_nAtomTypes << "\n";            
      for (unsigned j = 0; j < m_nAtomTypes; j++) {
        for (unsigned i = j + 1; i < m_nAtomTypes; i++) {
          m_neighborList[ncnt] = new NeighborList(pointList[i], pointList[j], false, m_pbc, getPbc(), m_neighborListCut[ncnt], m_neighborListStride[ncnt]);
          ncnt += 1;
        }
      }
    }
  } else {
    log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
    m_neighborListAll = new NeighborList(listall, m_pbc, getPbc());
    for (unsigned j = 0; j < m_numberBlocks; j++) {
      m_neighborList[j] = new NeighborList(pointList[j], pointList[j], true, m_pbc, getPbc());
    }
  }
  // Output neighborlist
  log << "Total Nlists: " << m_numberBlocks << " \n";
  for (unsigned j = 0; j < m_numberBlocks; j++) {
    log << "  list " << j + 1 << "   size " << m_neighborList[j]->size() << " \n";
  }
  // Calculate COM masses once and for all from lists
  if(m_centerOfMass) {
    //log << "Computing COM masses  \n";
    for (unsigned j = 0; j<m_positionCOM.size(); j++) {
      double massCOM = 0.;
      for (auto& atom : m_neighborlistCOM[j]->getFullAtomList()) {
        massCOM += myPDB.getOccupancy()[atom.index()];
      }
      for(auto& atom : m_neighborlistCOM[j]->getFullAtomList()) {
        if(massCOM > 0.) {
          m_massFactor[atom.index()] = myPDB.getOccupancy()[atom.index()] / massCOM;
        } else {
          m_massFactor[atom.index()] = 1.;
        }
      }
    }
  }

  //build box vectors and correct for pbc
  log << "Building the box from PDB data ... \n";
  Tensor box = myPDB.getBoxVec();
  log << "  Done! A,B,C vectors in Cartesian space:  \n";
  log.printf("  A:  %12.6f%12.6f%12.6f\n", box[0][0], box[0][1], box[0][2]);
  log.printf("  B:  %12.6f%12.6f%12.6f\n", box[1][0], box[1][1], box[1][2]);
  log.printf("  C:  %12.6f%12.6f%12.6f\n", box[2][0], box[2][1], box[2][2]);
  log << "Changing the PBC according to the new box \n";
  Pbc myPbc;
  myPbc.setBox(box);
  log << "The box volume is " << myPbc.getBox().determinant() << " \n";

  //Compute scaling factor
  if(m_doScaleVolume) {
    m_volumeFactor = cbrt(m_volume0 / myPbc.getBox().determinant());
    log << "Scaling atom distances by  " << m_volumeFactor << " \n";
  } else {
    log << "Using unscaled atom distances \n";
  }


  // build COMs from positions if requested
  if(m_centerOfMass) {
    for(unsigned j = 0; j < m_positionCOM.size(); j++) {
      m_positionCOM[j][0] = 0.;
      m_positionCOM[j][1] = 0.;
      m_positionCOM[j][2] = 0.; /*
      for(unsigned i = 0; i < m_neighborlistCOM[j]->getFullAtomList().size(); i++) {
        unsigned atomIndex = m_neighborlistCOM[j]->getFullAtomList()[i].index();
        m_positionCOM[j] += m_massFactor[atomIndex] * myPDB.getPositions()[atomIndex];
      } */
      for (auto& atom : m_neighborlistCOM[j]->getFullAtomList()) {
        m_positionCOM[j] += m_massFactor[atom.index()] * myPDB.getPosition(atom.index());
      }
    }
  }
  // build the rPIV2 distances (transformation and sorting is done afterwards)
  if (m_computeDerivatives) {
    log << "  PIV2  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
  }
  for(unsigned j = 0; j < m_numberBlocks; j++) {
    for(unsigned i = 0; i < m_neighborList[j]->size(); i++) {
      unsigned i0 = (m_neighborList[j]->getClosePairAtomNumber(i).first).index();
      unsigned i1 = (m_neighborList[j]->getClosePairAtomNumber(i).second).index();
      //calculate/get COM position of centers i0 and i1
      Vector position0, position1;
      if(m_centerOfMass) {
        //if(m_pbc) makeWhole();
        position0 = m_positionCOM[i0];
        position1 = m_positionCOM[i1];
      } else {
        position0 = myPDB.getPositions()[i0];
        position1 = myPDB.getPositions()[i1];
      }
      Vector ddist;
      if(m_pbc) {
        ddist = myPbc.distance(position0, position1);
      } else {
        ddist = delta(position0, position1);
      }
      double df = 0.;
      // Transformation and sorting done at the first timestep to solve the r0 definition issue
      if (m_computeDerivatives) {
        m_referencePIV2[j].push_back(m_switchingFunctions[j].calculate(ddist.modulo() * m_volumeFactor, df));
      } else {
        m_referencePIV2[j].push_back(ddist.modulo() * m_volumeFactor);
      }
    }
    if (m_computeDerivatives) {
      if(m_doSort[j]) {
        std::sort(m_referencePIV2[j].begin(), m_referencePIV2[j].end());
      }
      int lmt0 = 0;
      int lmt1 = 0;
      for (unsigned i = 0; i < m_referencePIV2[j].size(); i++) {
        if (int(m_referencePIV2[j][i] * double(m_nPrecision - 1)) == 0) {
          lmt0 += 1;
        }
        if (int(m_referencePIV2[j][i] * double(m_nPrecision - 1)) == 1) {
          lmt1 += 1;
        }
      }
      log.printf("       |%10i|%15i|%15i|%15i|\n", j, m_referencePIV2[j].size(), lmt0, lmt1);
    }
  }

  checkRead();
  
  // add N distance for the N reference states
  for (unsigned i = 0; i < m_numberReferences; i++) {
    std::string name = 
    addComponentWithDerivatives("d1");
    componentIsNotPeriodic("d1");
  }
  // add the lambda component for the path collective variables
  addComponent("lambda");
  componentIsNotPeriodic("lambda");
  /* // From the plumed manual on how to build-up a new Colvar
  addValueWithDerivatives();
  requestAtoms(m_neighborListAll->getFullAtomList());
  setNotPeriodic();
  // getValue()->setPeridodicity(false); */

  // set size of derivative vector
  m_deriv.resize(getNumberOfAtoms()); 
}

//---------------------------------------------------------
// DESTRUCTOR
//---------------------------------------------------------

PIV2::~PIV2()
{
  for (unsigned j = 0; j < m_numberBlocks; j++) {
    delete m_neighborList[j];
  }
  if(m_centerOfMass) {
    for (unsigned j = 0; j < m_numberBlocks; j++) {
      delete m_neighborlistCOM[j];
    }
  }
  delete m_neighborListAll;
}

//---------------------------------------------------------
// DISTANCE 
//---------------------------------------------------------

Vector PIV2::distanceAB(Vector A, Vector B)
{
  if(m_pbc) {
    return pbcDistance(A, B);
  } else {
    return delta(A, B);
  }
}

void PIV2::computeReferencePIV()
{}

//---------------------------------------------------------
// CALCULATE
//---------------------------------------------------------

void PIV2::calculate()
{
  // Local varaibles
  static int prev_stp = -1;
  // static int init_stp = 1;
  static std::vector<std::vector<Vector> > prev_pos(m_numberBlocks);
  static std::vector<std::vector<double> > cPIV2(m_numberBlocks);
  static std::vector<std::vector<int> > Atom0(m_numberBlocks);
  static std::vector<std::vector<int> > Atom1(m_numberBlocks);
  std::vector<std::vector<int> > A0(m_nPrecision);
  std::vector<std::vector<int> > A1(m_nPrecision);
  unsigned stride = 1;
  unsigned rank = 0;

  if(!m_serial) {
    stride = comm.Get_size();
    rank = comm.Get_rank();
  }

  // Transform (and sort) the rPIV2 before starting the dynamics
  /*if (((prev_stp == -1) || (init_stp == 1)) && !m_computeDerivatives) {
    if(prev_stp != -1) {
      init_stp = 0;
    } */
  if (!m_referenceIsInitialized && !m_computeDerivatives) {
    // Calculate the volume scaling factor
    if (m_doScaleVolume) {
      m_volumeFactor = cbrt(m_volume0 / getBox().determinant());
    }
    //Set switching function parameters
    log << "\n";
    log << "REFERENCE PDB # " << prev_stp + 2 << " \n";
    // now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    m_switchingFunctions.resize(m_numberBlocks);
    std::string errors;
    for (unsigned j = 0; j < m_numberBlocks; j++) {
      if(m_doScaleVolume) {
        double r0;
        std::vector<std::string> data = Tools::getWords(sw[j]);
        data.erase(data.begin());
        bool tmp = Tools::parse(data, "R_0", r0);
        std::string old_r0; 
        Tools::convert(r0, old_r0);
        r0 *= m_volumeFactor;
        std::string new_r0; 
        Tools::convert(r0, new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos + 4, old_r0.size(), new_r0);
      }
      m_switchingFunctions[j].set(sw[j], errors);
      std::string num;
      Tools::convert(j + 1, num);
      if (errors.length() != 0){
        error("problem reading SWITCH" + num + " keyword : " + errors );
      }
      m_r00[j] = m_switchingFunctions[j].get_r0();
      log << "  Swf: " << j << "  r0 = " << (m_switchingFunctions[j].description()).c_str() << " \n";
    }
    //Transform and sort
    log << "Building Reference PIV2 Vector \n";
    log << "  PIV2  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    double df = 0.;
    for (unsigned j = 0; j < m_numberBlocks; j++) {
      for (unsigned i = 0; i < m_referencePIV2[j].size(); i++) {
        m_referencePIV2[j][i] = m_switchingFunctions[j].calculate(m_referencePIV2[j][i], df);
      }
      if(m_doSort[j]) {
        std::sort(m_referencePIV2[j].begin(), m_referencePIV2[j].end());
      }
      int lmt0 = 0;
      int lmt1 = 0;
      for (unsigned i = 0; i < m_referencePIV2[j].size(); i++) {
        if (int(m_referencePIV2[j][i] * double(m_nPrecision - 1)) == 0) {
          lmt0 += 1;
        }
        if (int(m_referencePIV2[j][i] * double(m_nPrecision - 1)) == 1) {
          lmt1 += 1;
        }
      }
      log.printf("       |%10i|%15i|%15i|%15i|\n", j, m_referencePIV2[j].size(), lmt0, lmt1);
    }
    m_referenceIsInitialized = true;
    log << "\n";
  } // building of the reference PIV2 

  // Do the sorting only once per timestep to avoid building the PIV2 N times for N referencePIV2 PDB structures!
  if ((getStep() > prev_stp && getStep() % m_updatePIV2 == 0) || m_computeDerivatives) {
    if (m_computeDerivatives) {
      log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV2 \n";
    }
    // build COMs from positions if requested
    if(m_centerOfMass) {
      if(m_pbc) makeWhole();
      for(unsigned j = 0; j < m_positionCOM.size(); j++) {
        m_positionCOM[j][0] = 0.;
        m_positionCOM[j][1] = 0.;
        m_positionCOM[j][2] = 0.;
        for(auto& atom : m_neighborlistCOM[j]->getFullAtomList()) {
          m_positionCOM[j] += m_massFactor[atom.index()] * getPosition(atom.index());
        }
      }
    }
    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (m_doNeighbor) {
      bool doUpdate = false;
      // For the first step build previous positions = actual positions
      //if (prev_stp == -1) {
      if (!m_neighborIsInitialized) {
        for (unsigned j = 0; j < m_numberBlocks; j++) {
          for (unsigned i = 0; i < m_neighborList[j]->getFullAtomList().size(); i++) {
            Vector position;
            if(m_centerOfMass) {
              position = m_positionCOM[i];
            } else {
              position = getPosition(m_neighborList[j]->getFullAtomList()[i].index());
            }
            prev_pos[j].push_back(position);
          }
        }
        m_neighborIsInitialized = true;
        doUpdate = true;
      }
      // Decide whether to update lists based on atom displacement, every stride
      std::vector<std::vector<Vector> > tmp_pos(m_numberBlocks);
      if (getStep() % m_neighborListAll->getStride() ==0) {
        for (unsigned j = 0; j < m_numberBlocks; j++) {
          for (unsigned i=0; i < m_neighborList[j]->getFullAtomList().size(); i++) {
            Vector position;
            if(m_centerOfMass) {
              position = m_positionCOM[i];
            } else {
              position = getPosition(m_neighborList[j]->getFullAtomList()[i].index());
            }
            tmp_pos[j].push_back(position);
            if (pbcDistance(tmp_pos[j][i], prev_pos[j][i]).modulo() >= m_neighborListSkin[j]) {
              doUpdate = true;
            }
          }
        }
      }
      // Update Nlists if needed
      if (doUpdate == true) {
        for (unsigned j = 0; j < m_numberBlocks; j++) {
          for (unsigned i = 0; i < m_neighborList[j]->getFullAtomList().size(); i++) {
            prev_pos[j][i] = tmp_pos[j][i];
          }
          m_neighborList[j]->update(prev_pos[j]);
          log << " Step " << getStep() << "  Neighbour lists updated " << m_neighborList[j]->size() << " \n";
        }
      }
    } // if do neighbor
    // Calculate the volume scaling factor
    if (m_doScaleVolume) {
      m_volumeFactor = cbrt(m_volume0 / getBox().determinant());
    }
    Vector ddist;

    // Global to local variables
    bool doserial = m_serial;
    // Build "Nlist" PIV2 blocks
    for(unsigned j = 0; j < m_numberBlocks; j++) {
      if(m_doSort[j]) {
        // from global to local variables to speedup the for loop with if statements
        bool docom = m_centerOfMass;
        bool dopbc = m_pbc;
        // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
        std::vector<int> OrdVec(m_nPrecision, 0);
        cPIV2[j].resize(0);
        Atom0[j].resize(0);
        Atom1[j].resize(0);

        // Building distances for the PIV2 vector at time t
        if(m_timer) stopwatch.start("1 Build cPIV2");

        // previous loop jump in memory according to the number of cores 
        // for(unsigned i = rank; i < m_neighborList[j]->size(); i += stride) {
        // new loop through memory without jump
        // If we have N cores and Na atoms, we need (Na/N + 1) process by cores
        unsigned processByCore = unsigned(m_neighborList[j]->size() / stride + 1);
        for (unsigned i = rank * processByCore; (i < (rank + 1) * processByCore) && (i < m_neighborList[j]->size()); i++) {
          unsigned i0 = (m_neighborList[j]->getClosePairAtomNumber(i).first).index();
          unsigned i1 = (m_neighborList[j]->getClosePairAtomNumber(i).second).index();
          Vector position0, position1;
          if(docom) {
            position0 = m_positionCOM[i0];
            position1 = m_positionCOM[i1];
          } else {
            position0 = getPosition(i0);
            position1 = getPosition(i1);
          }
          ddist = distanceAB(position0, position1); /*
          if(dopbc) {
            ddist = pbcDistance(position0, position1);
          } else {
            ddist = delta(position0, position1);
          } */

          double df = 0.;
          //Integer sorting ... faster!
          //Transforming distances with the Switching function + real to integer transformation
          int Vint = int(m_switchingFunctions[j].calculate(ddist.modulo() * m_volumeFactor, df) * double(m_nPrecision - 1) + 0.5);
          //Integer transformed distance values as index of the Ordering Vector OrdVec
          OrdVec[Vint] += 1;
          //Keeps track of atom indices for force and virial calculations
          A0[Vint].push_back(i0);
          A1[Vint].push_back(i1);
        }

        if(m_timer) stopwatch.stop("1 Build cPIV2");
        if(m_timer) stopwatch.start("2 Sort cPIV2");

        if(!doserial) {
          // Vectors keeping track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          std::vector<int> Vdim(stride,0);
          std::vector<int> Vpos(stride,0);
          // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
          std::vector<int> OrdVecAll(stride * m_nPrecision);
          // Big vectors contining all Atom indexes for every occupancy (Atom0O(Nprec,n) and Atom1O(Nprec,n) matrices in one vector)
          std::vector<int> Atom0F;
          std::vector<int> Atom1F;
          // Vector used to reconstruct arrays
          std::vector<unsigned> k(stride,0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for(unsigned i = 1; i < m_nPrecision; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            // Can this be avoided ???
            Atom0F.insert(Atom0F.end(), A0[i].begin(), A0[i].end());
            Atom1F.insert(Atom1F.end(), A1[i].begin(), A1[i].end());
            A0[i].resize(0);
            A1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV2 block
          A0[0].resize(0);
          A1[0].resize(0);
          A0[m_nPrecision - 1].resize(0);
          A1[m_nPrecision - 1].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          OrdVec[0] = 0;
          OrdVec[m_nPrecision - 1] = 0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim=Atom0F.size();
          // Vdim and Vpos keep track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim, 1, &Vdim[0], 1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int Fdim = 0;
          for(unsigned i = 1; i < stride; i++) {
            Vpos[i] = Vpos[i-1] + Vdim[i-1];
            Fdim += Vdim[i];
          }
          Fdim += Vdim[0];
          // build big vectors for atom pairs on all ranks for all ranks
          std::vector<int> Atom0FAll(Fdim);
          std::vector<int> Atom1FAll(Fdim);
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          //   Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV2
          comm.Allgather(&OrdVec[0], m_nPrecision, &OrdVecAll[0], m_nPrecision);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&Atom0F[0],Atom0F.size(),&Atom0FAll[0],&Vdim[0],&Vpos[0]);
          comm.Allgatherv(&Atom1F[0],Atom1F.size(),&Atom1FAll[0],&Vdim[0],&Vpos[0]);

          if(m_timer) stopwatch.stop("2 Sort cPIV2");
          if(m_timer) stopwatch.start("3 Reconstruct cPIV2");

          // Reconstruct the full vectors from collections of Allgathered parts (this is a serial step)
          // This is the tricky serial step, to assemble toghether PIV2 and atom-pair info from head-tail big vectors
          // Loop before on l and then on i would be better but the allgather should be modified
          // Loop on blocks
          // for(unsigned m = 0; m < Nlist; m++) {
          // Loop on Ordering Vector size excluding zeros (i=1)
          for(unsigned i = 1; i < m_nPrecision; i++) {
            // Loop on the ranks
            for(unsigned l = 0; l < stride; l++) {
              // Loop on the number of head-to-tail pieces
              for(unsigned m = 0; m < OrdVecAll[i + l * m_nPrecision]; m++) {
                // cPIV2 is the current PIV2 at time t
                cPIV2[j].push_back(double(i)/double(m_nPrecision - 1));
                Atom0[j].push_back(Atom0FAll[k[l]+Vpos[l]]);
                Atom1[j].push_back(Atom1FAll[k[l]+Vpos[l]]);
                k[l]+=1;
              }
            }
          }
          if(m_timer) stopwatch.stop("3 Reconstruct cPIV2");
        } else {
          for(unsigned i = 1; i < m_nPrecision; i++) {
            for(unsigned m = 0; m < OrdVec[i]; m++) {
              cPIV2[j].push_back(double(i) / double(m_nPrecision - 1));
              Atom0[j].push_back(A0[i][m]);
              Atom1[j].push_back(A1[i][m]);
            }
          }
        }
      }
    }
  }

  Vector distance;
  double dfunc=0.;
  // Calculate volume scaling factor
  if (m_doScaleVolume) {
    m_volumeFactor = cbrt(m_volume0 / getBox().determinant());
  }

  // This test may be run by specifying the TEST keyword as input, it pritnts referencePIV2 and cPIV2 and quits
  if(m_test) {
    unsigned limit = 0;
    for(unsigned j = 0; j < m_numberBlocks; j++) {
      if(m_doSort[j]) {
        limit = cPIV2[j].size();
      } else {
        limit = m_referencePIV2[j].size();
      }
      log.printf("PIV2 Block:  %6i %12s %6i \n", j, "      Size:", limit);
      log.printf("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV2   ","    r-PIV2   ","   i-j distance vector       ");
      for(unsigned i = 0; i < limit; i++) {
        unsigned i0 = 0;
        unsigned i1 = 0;
        if(m_doSort[j]) {
          i0 = Atom0[j][i];
          i1 = Atom1[j][i];
        } else {
          i0 = (m_neighborList[j]->getClosePairAtomNumber(i).first).index();
          i1 = (m_neighborList[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector position0,position1;
        if(m_centerOfMass) {
          position0 = m_positionCOM[i0];
          position1 = m_positionCOM[i1];
        } else {
          position0 = getPosition(i0);
          position1 = getPosition(i1);
        }
        distance = distanceAB(position0, position1);/*
        if(m_pbc) {
          distance = pbcDistance(position0,position1);
        } else {
          distance = delta(position0,position1);
        }*/

        dfunc=0.;
        double cP, rP;
        if(m_doSort[j]) {
          cP = cPIV2[j][i];
          rP = m_referencePIV2[j][m_referencePIV2[j].size() - cPIV2[j].size()+i];
        } else {
          double dm = distance.modulo();
          cP = m_switchingFunctions[j].calculate(dm * m_volumeFactor, dfunc);
          rP = m_referencePIV2[j][i];
        }
        log.printf("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
      }
    }
    log.printf("This was a test, now exit \n");
    exit();
  }

  if(m_timer) stopwatch.start("4 Build For Derivatives");
  
  // non-global variables Nder and Scalevol defined to speedup if structures in cycles
  bool Nder = m_computeDerivatives;
  bool Scalevol = m_doScaleVolume;
  if(getStep() % m_updatePIV2 == 0) {
    // set to zero PIV2distance, derivatives and virial when they are calculated
    for(unsigned j = 0; j < m_deriv.size(); j++) {
      for(unsigned k = 0; k < 3; k++) {
        m_deriv[j][k] = 0.;
      }
    }
    for(unsigned j = 0; j < 3; j++) {
      for(unsigned k = 0; k < 3; k++) {
        m_virial[j][k] = 0.;
      }
    }
    m_PIV2distance = 0.;
    // Re-compute atomic distances for derivatives and compute PIV2-PIV2 distance
    for(unsigned j = 0; j < m_numberBlocks; j++) {
      unsigned limit=0;
      // dosorting definition is to speedup if structure in cycles with non-global variables
      bool dosorting = m_doSort[j];
      bool docom = m_centerOfMass;
      bool dopbc = m_pbc;
      if(dosorting) {
        limit = cPIV2[j].size();
      } else {
        limit = m_referencePIV2[j].size();
      }
      unsigned processByCore = unsigned(limit / stride + 1);
      for (unsigned i = rank * processByCore; (i < (rank + 1) * processByCore) && (i < limit); i++) {
        unsigned i0 = 0;
        unsigned i1 = 0;
        if(dosorting) {
          i0 = Atom0[j][i];
          i1 = Atom1[j][i];
        } else {
          i0 = (m_neighborList[j]->getClosePairAtomNumber(i).first).index();
          i1 = (m_neighborList[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector position0, position1;
        if(docom) {
          position0 = m_positionCOM[i0];
          position1 = m_positionCOM[i1];
        } else {
          position0 = getPosition(i0);
          position1 = getPosition(i1);
        }
        distance = distanceAB(position0, position1);/*
        if(dopbc) {
          distance = pbcDistance(position0,position1);
        } else {
          distance = delta(position0,position1);
        }*/
        dfunc = 0.;
        // this is needed for dfunc and dervatives
        double dm = distance.modulo();
        double tPIV2 = m_switchingFunctions[j].calculate(dm * m_volumeFactor, dfunc);
        // PIV2 distance
        double coord = 0.;
        if(!dosorting || Nder) {
          coord = tPIV2 - m_referencePIV2[j][i];
        } else {
          coord = cPIV2[j][i] - m_referencePIV2[j][m_referencePIV2[j].size() - cPIV2[j].size()+i];
        }
        // Calculate derivatives, virial, and variable = sum_j (scaling[j] *(cPIV2-rPIV2)_j^2)
        // WARNING: dfunc=dswf/(m_volumeFactor * dm)  (this may change in future Plumed versions)
        double tmp = 2. * m_blockScaling[j] * coord * m_volumeFactor * m_volumeFactor * dfunc;
        Vector tmpder = tmp*distance;
        // 0.5*(x_i-x_k)*f_ik         (force on atom k due to atom i)
        if(docom) {
          Vector dist;
          for(auto& atom0 : m_neighborlistCOM[i0]->getFullAtomList()) {
            unsigned x0 = atom0.index();
            m_deriv[x0] -= tmpder * m_massFactor[x0];
            for(unsigned l = 0; l < 3; l++) {
              dist[l] = 0.;
            }
            Vector position0 = getPosition(x0);
            for(auto& atom1 : m_neighborlistCOM[i0]->getFullAtomList()) {
              dist += distanceAB(position0, getPosition(atom1.index()) );/*
              unsigned x1 = atom1.index();
              Vector position1 = getPosition(x1);
              if(dopbc) {
                dist += pbcDistance(position0, position1);
              } else {
                dist += delta(position0, position1);
              }*/
            }
            for(auto& atom1 : m_neighborlistCOM[i1]->getFullAtomList()) {
              dist += distanceAB(position0, getPosition(atom1.index()) );/*
              unsigned x1 = atom1.index();
              Vector position1 = getPosition(x1);
              if(dopbc) {
                dist += pbcDistance(position0, position1);
              } else {
                dist += delta(position0, position1);
              }*/
            }
            m_virial -= 0.25 * m_massFactor[x0] * Tensor(dist, tmpder);
          } // loop on first atom of each pair
          for(auto& atom1 : m_neighborlistCOM[i1]->getFullAtomList()) {
            unsigned x1 = atom1.index();
            m_deriv[x1] += tmpder * m_massFactor[x1];
            for(unsigned l = 0; l < 3; l++) {
              dist[l] = 0.;
            }
            Vector position1 = getPosition(x1);
            for(auto& atom0 : m_neighborlistCOM[i1]->getFullAtomList()) {
              dist += distanceAB(position1, getPosition(atom0.index()) );/*
              unsigned x0 = atom0.index();
              Vector position0 = getPosition(x0);
              if(dopbc) {
                dist += pbcDistance(position1, position0);
              } else {
                dist += delta(position1, position0);
              }*/
            }
            for(auto& atom0 : m_neighborlistCOM[i0]->getFullAtomList()) {
              dist += distanceAB(position1, getPosition(atom0.index()) );/*
              unsigned x0 = atom0.index();
              Vector position0 = getPosition(x0);
              if(dopbc) {
                dist += pbcDistance(position1, position0);
              } else {
                dist += delta(position1, position0);
              }*/
            }
            m_virial += 0.25 * m_massFactor[x1] * Tensor(dist, tmpder);
          } // loop on second atom of each pair
        } else { // if do center of mass
          m_deriv[i0] -= tmpder;
          m_deriv[i1] += tmpder;
          m_virial    -= tmp * Tensor(distance, distance);
        }
        if(Scalevol) {
          m_virial += 1./3. * tmp * dm * dm * Tensor::identity();
        }
        m_PIV2distance += m_blockScaling[j] * coord * coord;
      }
    } // loop on all atoms

    if (!m_serial) {
      comm.Barrier();
      comm.Sum(&m_PIV2distance, 1);
      if(!m_deriv.empty()) {
        comm.Sum(&m_deriv[0][0], 3*m_deriv.size());
      }
      comm.Sum(&m_virial[0][0], 9);
    }
  } // if update_piv

  prev_stp = getStep();

  //Timing
  if(m_timer) stopwatch.stop("4 Build For Derivatives");
  if(m_timer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log << stopwatch;
  }

  // Update derivatives, virial, and variable (PIV2-distance^2)
  for(unsigned i = 0; i < m_deriv.size(); ++i) setAtomsDerivatives(i, m_deriv[i]);
  setValue           (m_PIV2distance);
  setBoxDerivatives  (m_virial);
} // end of calculate

} // close namespace piv
} // close namespace PLMD
