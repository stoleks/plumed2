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
  void checkFieldsAllowed() {}
private:
  Vector distanceAB(Vector A, Vector B);
  void setCutOff();
private:
  ForwardDecl<Stopwatch> stopwatch_fwd;
  /// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch = *stopwatch_fwd;
  bool m_pbc, m_serial, m_timer, m_refIsIni, m_neighborIsIni; 
  bool m_doCutOff, m_doScaleVolume, m_cross, m_direct, m_doNeighbor;
  bool m_doTest, m_computeDerivatives, m_centerOfMass;
  unsigned m_nPrecision, m_nAtomTypes, m_numberBlocks; //, m_neighborListSize, 
  unsigned m_numberReferences, m_updatePIV2;
  double m_volumeFactor, m_volume0, m_lambda;
  Tensor m_virial;
  NeighborList* m_neighborListAll;
  std::vector<bool> m_doSort;
  std::vector<double> m_blockScaling, m_r00;
  std::vector<double> m_neighborListSkin;
  std::vector<double> m_massFactor;
  std::vector<double> m_cutOff;
  std::vector<double> m_distancePIV2;
  std::vector<Vector> m_positionCOM;
  std::vector<Vector> m_deriv;
  std::vector<std::string> m_switchData;
  std::vector<NeighborList*> m_neighborList;
  std::vector<NeighborList*> m_neighborlistCOM;
  std::vector<SwitchingFunction> m_switchFunc;
  std::vector<std::vector<Vector>> m_prevPosition;
  // first vector for reference, second for block, third for atoms
  std::vector<std::vector<std::vector<double>>> m_refPIV2;
};

PLUMED_REGISTER_ACTION(PIV2,"PIV2")

void PIV2::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered", "SWITCH", "The switching functions parameter."
           "You should specify a Switching function for all PIV2 blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("numbered", "REF_FILE", "PDB file name that contains the i-th reference structure.");
  keys.add("compulsory", "N_REFERENCE", "The number of reference structure."
           "It should match with the number of REF_FILE");
  keys.add("compulsory", "PRECISION", "the precision for approximating reals with integers in sorting.");
  keys.add("compulsory", "PIV2ATOMS", "Number of atoms to use for PIV2.");
  keys.add("compulsory", "SORT", "Whether to sort or not the PIV2 block.");
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
  keys.addFlag("CUTOFF", false, "Use cut-off to reduce the computational cost of the PIV.");
  keys.reset_style("SWITCH", "compulsory");
  componentsAreNotOptional(keys);
  // output
  keys.addOutputComponent("lambda", "default", "optimal lambda needed for the pathCV");
  keys.addOutputComponent("di", "default", "PIV distance between the i-th"
                          " reference state and the current state");
}

//---------------------------------------------------------
// CONSTRUCTOR
//---------------------------------------------------------

PIV2::PIV2(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  m_pbc(true),
  m_serial(false),
  m_timer(false),
  m_refIsIni(false),
  m_neighborIsIni(false),
  m_doScaleVolume(false),
  m_cross(true),
  m_direct(true),
  m_doNeighbor(false),
  m_doTest(false),
  m_computeDerivatives(false),
  m_centerOfMass(false),
  m_nPrecision(1000),
  m_nAtomTypes(1),
  // m_neighborListSize(1),
  m_numberReferences(1),
  m_updatePIV2(1),
  m_volumeFactor(1.),
  m_volume0(1.),
  m_lambda(1.)
{
  log << "Starting PIV2 Constructor\n";
  bool onlyCross = false, onlyDirect = false, noPbc = !m_pbc;

  // parse all the mandatory inputs that are not vector
  parse("VOLUME", m_volume0);
  parse("PIV2ATOMS", m_nAtomTypes);
  parse("PRECISION", m_nPrecision);
  parse("N_REFERENCE", m_numberReferences);
  // parse all the options
  parseFlag("TEST", m_doTest);
  parseFlag("NOPBC", noPbc);
  parseFlag("TIMER", m_timer);
  parseFlag("SERIAL", m_serial);
  parseFlag("NLIST", m_doNeighbor);
  parseFlag("CUTOFF", m_doCutOff);
  parseFlag("ONLYCROSS", onlyCross);
  parseFlag("ONLYDIRECT", onlyDirect);
  parseFlag("COM", m_centerOfMass);
  parseFlag("DERIVATIVES", m_computeDerivatives);
  // parse the atom names 
  std::vector<std::string> atomTypes(m_nAtomTypes);
  parseVector("ATOMTYPES", atomTypes);

  // Stride for which the PIV are computed
  if (keywords.exists("UPDATEPIV2")){
    parse("UPDATEPIV2", m_updatePIV2);
  }
  // Precision on the real-to-integer transformation for the sorting
  if (m_nPrecision < 2) { 
    error("Precision must be => 2");
  }
  // PBC
  m_pbc = !noPbc;
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
  if (m_timer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }
  // Center of Mass
  if (m_centerOfMass){
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
  if (onlyCross) {
    m_direct = false;
    log << "Using only CROSS-PIV2 blocks\n";
  } 
  if (onlyDirect) {
    m_cross = false;
    log << "Using only DIRECT-PIV2 blocks\n";
  }
  m_numberBlocks = 0;
  if (m_cross) {
    m_numberBlocks += unsigned(double(m_nAtomTypes * (m_nAtomTypes - 1)) / 2.);
  }
  if (m_direct) {
    m_numberBlocks += unsigned(m_nAtomTypes);
  }

  // resizing all class vector according to m_numberBlocks
  m_distancePIV2.resize(m_numberReferences);
  m_refPIV2.resize(m_numberReferences);
  for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
    m_refPIV2[iRef].resize(m_numberBlocks);
  }
  m_doSort.resize(m_numberBlocks);
  m_blockScaling.resize(m_numberBlocks);
  m_r00.resize(m_numberBlocks);
  m_switchData.resize(m_numberBlocks); 
  m_neighborList.resize(m_numberBlocks);
  m_neighborListSkin.resize(m_numberBlocks);
  m_prevPosition.resize(m_numberBlocks);

  // setting neighborlist parameters 
  std::vector<double> neighborListCut(m_numberBlocks, 0.);
  std::vector<int> neighborListStride(m_numberBlocks, 0);
  if (m_doNeighbor) {
    parseVector("NL_CUTOFF", neighborListCut);
    parseVector("NL_STRIDE", neighborListStride);
    parseVector("NL_SKIN", m_neighborListSkin);
    for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      if (neighborListCut[iBlock] <= 0.0) {
        error("NL_CUTOFF should be explicitly specified and positive");
      }
      if (neighborListStride[iBlock] <= 0) {
        error("NL_STRIDE should be explicitly specified and positive");
      }
      if (m_neighborListSkin[iBlock] <= 0.) {
        error("NL_SKIN should be explicitly specified and positive");
      }
      neighborListCut[iBlock] = neighborListCut[iBlock] + m_neighborListSkin[iBlock];
      log << "For Block " << iBlock + 1 
          << ", neighbor list cut-off=" << neighborListCut[iBlock]
          << ", stride=" << neighborListStride[iBlock]
          << ", shell=" << m_neighborListSkin[iBlock] << "\n";
    }
  }

  // Sorting
  std::vector<unsigned> yesNoSort(m_numberBlocks);
  parseVector("SORT", yesNoSort);
  for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
    if (yesNoSort[iBlock] == 0 || m_computeDerivatives) {
      m_doSort[iBlock] = false;
      log << "Not sorting block " << iBlock + 1 << ". ";
    } else {
      m_doSort[iBlock] = true;
      log << "Sorting block " << iBlock + 1 << ". ";
    }
  }
  log << "\n";

  // PIV2 scaled option
  for(unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
    m_blockScaling[iBlock] = 1.;
  }
  if (keywords.exists("SFACTOR")) {
    parseVector("SFACTOR", m_blockScaling);
  }

  // read parameters and set-up switching functions here only if computing derivatives
  for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
    if ( !parseNumbered("SWITCH", iBlock + 1, m_switchData[iBlock]) ){
      log << "Problem while reading the switching function parameters.\n";
      break;
    }
  }
  if (m_computeDerivatives) {
    log << "Switching Function Parameters \n";
    m_switchFunc.resize(m_numberBlocks);
    std::string errors;
    for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
     // std::string num;
     // Tools::convert(iBlock + 1, num);
      m_switchFunc[iBlock].set(m_switchData[iBlock], errors);
      if (errors.length() != 0){
        error("problem reading switch" + std::to_string(iBlock + 1) + " keyword : " + errors );
      }
      m_r00[iBlock] = m_switchFunc[iBlock].get_r0();
      log << "  Swf: " << iBlock + 1 << "  r0=" << (m_switchFunc[iBlock].description()).c_str() << "\n";
    }
  }

  // set the cut-off
  setCutOff();

  // Reference PDB file 
  for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
    std::string referenceFile; 
    if ( !parseNumbered("REF_FILE", iRef + 1, referenceFile) ){
      log << "Error while trying to read reference " << iRef << " " 
          << referenceFile << "\n";
      break;
    } 
    log << "\n----- Reference " << iRef + 1 << " -----\n"; 

    // opening of the reference PBD file
    PDB myPDB;
    FILE* fp = fopen(referenceFile.c_str(), "r");
    if (fp != NULL) {
      log << "Opening PDB file with reference frame: " << referenceFile << "\n";
      myPDB.readFromFilepointer(fp, plumed.getAtoms().usingNaturalUnits(), 
                                0.1 / atoms.getUnits().getLength());
      fclose (fp);
    } else {
      error("Error in reference PDB file " + referenceFile);
    }
 
    // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
    std::vector<std::vector<AtomNumber>> pairList(m_nAtomTypes);
    std::vector<std::vector<AtomNumber>> comAtoms(1);
    // number of atoms in the PDB cell
    unsigned neighborListSize = myPDB.getAtomNumbers().size();
    log << "Atoms numbers " << myPDB.getAtomNumbers().size() << "\n";

    // Build residueIndex
    std::vector<unsigned> residueIndex;
    for (unsigned iAtm = 0; iAtm < myPDB.getAtomNumbers().size(); iAtm++) {
      unsigned residueNumber = myPDB.getResidueNumber(myPDB.getAtomNumbers()[iAtm]);
      if (std::find(residueIndex.begin(), residueIndex.end(), residueNumber)
          == residueIndex.end()) {
        residueIndex.push_back(residueNumber);
      }
    }
    unsigned nResidues = residueIndex.size();
 
    // indexList is the index of atom/COM used in neighborlists (for COM indexList
    // is the index of the first atom in the pdb belonging to that COM)
    unsigned indexListSize;
    if (m_centerOfMass) {
      indexListSize = nResidues;
    } else {
      indexListSize = neighborListSize; 
    }
    std::vector<unsigned> indexList(indexListSize);

    // If COM resize important arrays
    comAtoms.resize(neighborListSize);
    if (m_centerOfMass) {
      m_neighborlistCOM.resize(neighborListSize);
      m_positionCOM.resize(neighborListSize);
      m_massFactor.resize(neighborListSize, 0.);
    }
    log << "Total COM/Atoms: " << m_nAtomTypes * nResidues << " \n";

    // Build lists of Atoms/COMs for neighbor lists
    // comAtoms filled also for non_COM calculation for analysis purposes
    for (unsigned iType = 0; iType < m_nAtomTypes; iType++) {
      unsigned outputID = 0;
      std::fill (indexList.begin(), indexList.end(), 0);
      // This builds lists for neighbor lists
      for (auto& atom : myPDB.getAtomNumbers()) {
        std::string atomName;
        unsigned pointIndex = 0;
        if (m_centerOfMass) {
          atomName = myPDB.getAtomName(atom);
          auto iRes = std::find (residueIndex.begin(), residueIndex.end(), myPDB.getResidueNumber(atom));
          if (iRes != residueIndex.end()) {
            pointIndex = iRes - residueIndex.begin();
          }
        } else {
          atomName = myPDB.getAtomName(atom);
          pointIndex = atom.index();
        }
        if (atomName == atomTypes[iType]) {
          if (indexList[pointIndex] == 0) {
            // adding the atomnumber to the atom/COM list for pairs
            pairList[iType].push_back(atom);
            indexList[pointIndex] = atom.index() + 1;
            outputID = pointIndex;
          }
          // adding the atomnumber to list of atoms for every COM/Atoms
          comAtoms[indexList[pointIndex] - 1].push_back(atom);
        }
      }
      // Output Lists
      log << "  Groups of type  " << iType << ": " << pairList[iType].size() << " \n";
      std::string groupName;
      unsigned groupSize;
      if (m_centerOfMass) {
        groupName = myPDB.getResidueName(comAtoms[indexList[outputID] - 1][0]);
        groupSize = comAtoms[indexList[outputID] - 1].size();
      } else {
        groupName = myPDB.getAtomName(comAtoms[indexList[outputID] - 1][0]);
        groupSize = 1;
      }
      log.printf("    %6s %3s %13s %10i %6s\n", "type  ", groupName.c_str(),
                 "   containing ", groupSize," atoms");
    }
 
    // This is to build the list with all the atoms
    std::vector<AtomNumber> listAllAtom;
    for (auto& atom : myPDB.getAtomNumbers()) {
      listAllAtom.push_back(atom);
    }
 
    if (m_doNeighbor) {
      log << "Creating Neighbor Lists \n";
      for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        log << "For Block " << iBlock + 1 
            << ", neighbor list cut-off=" << neighborListCut[iBlock]
            << ", stride=" << neighborListStride[iBlock]
            << ", shell=" << m_neighborListSkin[iBlock] << "\n";
      }
      // WARNING: is neighborListCut meaningful here?
      m_neighborListAll= new NeighborList(listAllAtom, m_pbc, getPbc(),
                               neighborListCut[0], neighborListStride[0]);
      //if (m_centerOfMass) {
      for (unsigned i = 0; i < m_positionCOM.size(); i++) {
        // WARNING: is neighborListCut meaningful here?
        m_neighborlistCOM[i] = new NeighborList(comAtoms[i], m_pbc, getPbc(),
                                     neighborListCut[0], neighborListStride[0]);
      }
      unsigned ncnt = 0;
      // Direct blocks AA, BB, CC, ...
      if (m_direct) {
        log << "Number of blocks: " << m_neighborList.size() << ", number of atom"
            << "types: " << m_nAtomTypes << "\n";
        for (unsigned iType = 0; iType < m_nAtomTypes; iType++) {
          m_neighborList[ncnt] = new NeighborList(pairList[iType], m_pbc, getPbc(),
                                       neighborListCut[iType], neighborListStride[iType]);
          ncnt += 1;
        }
      }
      // Cross blocks AB, AC, BC, ...
      if (m_cross) {
        log << "Number of blocks: " << m_neighborList.size() << ", number of atom"
            << "types: " << m_nAtomTypes << "\n";
        for (unsigned iType1 = 0; iType1 < m_nAtomTypes; iType1++) {
          for (unsigned iType2 = iType1 + 1; iType2 < m_nAtomTypes; iType2++) {
            m_neighborList[ncnt] = new NeighborList(pairList[iType1], pairList[iType2],
                                         false, m_pbc, getPbc(),
                                         neighborListCut[ncnt], neighborListStride[ncnt]);
            ncnt += 1;
          }
        }
      }
    } else {
      log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
      m_neighborListAll = new NeighborList(listAllAtom, m_pbc, getPbc());
      for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        m_neighborList[iBlock] = new NeighborList(pairList[iBlock], pairList[iBlock],
                                                  true, m_pbc, getPbc());
      }
    }
    // Output neighborlist
    log << "Total Nlists: " << m_numberBlocks << " \n";
    for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      log << "  list " << iBlock + 1 << "   size " << m_neighborList[iBlock]->size() << " \n";
    }
    // Calculate COM masses once and for all from lists
    if (m_centerOfMass) {
      //log << "Computing Center Of Mass masses  \n";
      for (unsigned iCOM = 0; iCOM < m_positionCOM.size(); iCOM++) {
        double massCOM = 0.;
        for (auto& atom : m_neighborlistCOM[iCOM]->getFullAtomList()) {
          massCOM += myPDB.getOccupancy()[atom.index()];
        }
        for(auto& atom : m_neighborlistCOM[iCOM]->getFullAtomList()) {
          if (massCOM > 0.) {
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
    if (m_doScaleVolume) {
      m_volumeFactor = cbrt(m_volume0 / myPbc.getBox().determinant());
      log << "Scaling atom distances by  " << m_volumeFactor << " \n";
    } else {
      log << "Using unscaled atom distances \n";
    }

    // build COMs from positions if requested
    if (m_centerOfMass) {
      for(unsigned iCOM = 0; iCOM < m_positionCOM.size(); iCOM++) {
        m_positionCOM[iCOM][0] = 0.;
        m_positionCOM[iCOM][1] = 0.;
        m_positionCOM[iCOM][2] = 0.; 
        for (auto& atom : m_neighborlistCOM[iCOM]->getFullAtomList()) {
          m_positionCOM[iCOM] += m_massFactor[atom.index()] * myPDB.getPositions()[atom.index()];
        }
      }
    }
    // build the rPIV2 distances (transformation and sorting is done afterwards)
    if (m_computeDerivatives) {
      log << "  PIV2  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    }
    for(unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      for(unsigned iNl = 0; iNl < m_neighborList[iBlock]->size(); iNl++) {
        unsigned i0 = (m_neighborList[iBlock]->getClosePairAtomNumber(iNl).first).index();
        unsigned i1 = (m_neighborList[iBlock]->getClosePairAtomNumber(iNl).second).index();
        //calculate/get COM position of centers i0 and i1
        Vector position0, position1;
        if (m_centerOfMass) {
          //if (m_pbc) makeWhole();
          position0 = m_positionCOM[i0];
          position1 = m_positionCOM[i1];
        } else {
          position0 = myPDB.getPositions()[i0];
          position1 = myPDB.getPositions()[i1];
        }
        Vector pairDist;
        if (m_pbc) {
          pairDist = myPbc.distance(position0, position1);
        } else {
          pairDist = delta(position0, position1);
        }
        double df = 0.;
        // Transformation and sorting done at the first timestep to solve the r0 definition issue
        if (pairDist.modulo2() < m_cutOff[iBlock] * m_cutOff[iBlock]) {
          if (m_computeDerivatives) {
            m_refPIV2[iRef][iBlock].push_back(m_switchFunc[iBlock].calculate(pairDist.modulo() * m_volumeFactor, df));
          } else {
            m_refPIV2[iRef][iBlock].push_back(pairDist.modulo() * m_volumeFactor);
          }
        }
      }
      log << "reference PIV block " << iBlock + 1 << " has size: " 
          << m_refPIV2[iRef][iBlock].size() << " over a total of "
           << m_neighborList[iBlock]->size() << " atoms-atoms pair\n";
      if (m_computeDerivatives) {
        if (m_doSort[iBlock]) {
          std::sort(m_refPIV2[iRef][iBlock].begin(), m_refPIV2[iRef][iBlock].end());
        }
        int lmt0 = 0;
        int lmt1 = 0;
        for (unsigned iAtm = 0; iAtm < m_refPIV2[iRef][iBlock].size(); iAtm++) {
          if (m_refPIV2[iRef][iBlock][iAtm] > 0.9) { lmt1++; }
          if (m_refPIV2[iRef][iBlock][iAtm] < 0.1) { lmt0++; }
        }
        log.printf("       |%10i|%15i|%15i|%15i|\n", iBlock, m_refPIV2[iRef][iBlock].size(), lmt0, lmt1);
      } // if we compute derivatives
    } // loop over the number of blocks
  } // loop over the number of references states
  log << "\n";

  checkRead();

  // add N distance for the N reference states
  for (unsigned iRef = 1; iRef <= m_numberReferences; iRef++) {
    addComponentWithDerivatives("d" + std::to_string(iRef));
    componentIsNotPeriodic("d" + std::to_string(iRef));
  }
  // add the lambda component for the path collective variables
  addComponent("lambda");
  componentIsNotPeriodic("lambda");
  requestAtoms(m_neighborListAll->getFullAtomList());

  // set size of derivative vector
  m_deriv.resize(getNumberOfAtoms()); 
}

//---------------------------------------------------------
// DESTRUCTOR
//---------------------------------------------------------

PIV2::~PIV2()
{
  for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
    delete m_neighborList[iBlock];
  }
  if (m_centerOfMass) {
    for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      delete m_neighborlistCOM[iBlock];
    }
  }
  delete m_neighborListAll;
}

//---------------------------------------------------------
// CUT-OFF
//---------------------------------------------------------

void PIV2::setCutOff(){
  // just resolves the equation swf(x) = 1 / nPrecision for all switching functions
  for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++){
    SwitchingFunction switchFunc;
    std::string errors;
    switchFunc.set(m_switchData[iBlock], errors);
    if (errors.length() != 0){
      error("problem reading switch" + std::to_string(iBlock + 1) + " keyword : " + errors );
    }

    log << "For block " << iBlock + 1;
    if (m_doCutOff) {
      unsigned maxIteration = 0;
      double f, g, df;
      double espilon = 1. / double(m_nPrecision);
      
      double x0 =  switchFunc.get_r0() / 2.;
      //log << switchFunc.description().c_str();
      //log << " x_0 =" << x0 ;
      f = switchFunc.calculate(x0, df) - epsilon;
      while (fabs(f) > (epsilon / 100.) && maxIteration < 10000)
      {
        f = switchFunc.calculate(x0, df) - espilon;
        g = (switchFunc.calculate(x0 + f, df) - epsilon) / f - 1;
        x0 = x0 - f / g;
        maxIteration++;
      }
      m_cutOff.push_back(x0);
      log << " epsilon=" << espilon << " this correspond to a cut-off=";
    } else {
      double dmax = switchFunc.get_dmax();
      m_cutOff.push_back(10.*dmax);
      log << ", cut-off=";
    }
    log << m_cutOff[iBlock] << "\n";
  }
}

//---------------------------------------------------------
// DISTANCE 
//---------------------------------------------------------

inline Vector PIV2::distanceAB(Vector A, Vector B)
{
  if (m_pbc) {
    return pbcDistance(A, B);
  } else {
    return delta(A, B);
  }
}

//---------------------------------------------------------
// CALCULATE
//---------------------------------------------------------

void PIV2::calculate()
{
  std::vector<std::vector<double>> cPIV2(m_numberBlocks);
  std::vector<std::vector<int>> atmI0(m_numberBlocks);
  std::vector<std::vector<int>> atmI1(m_numberBlocks);
  std::vector<std::vector<int>> atmPrecI0(m_nPrecision);
  std::vector<std::vector<int>> atmPrecI1(m_nPrecision);
  unsigned stride = 1;
  unsigned rank = 0;

  if (!m_serial) {
    stride = comm.Get_size();
    rank = comm.Get_rank();
  }

  // Calculate volume scaling factor
  if (m_doScaleVolume) {
    m_volumeFactor = cbrt(m_volume0 / getBox().determinant());
  }

  // Transform (and sort) the rPIV2 before starting the dynamics
  if (!m_refIsIni && !m_computeDerivatives) {
    // Calculate the volume scaling factor
    // compute the N reference PIV
    for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
      //Set switching function parameters
      log << "\n";
      log << "REFERENCE PDB # " << iRef + 1 << " \n";
      log << "Switching Function Parameters \n";
      m_switchFunc.resize(m_numberBlocks);
      std::string errors;
      // setting the PIV reference for each atom pair blocks
      for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        if (m_doScaleVolume) {
          auto data = Tools::getWords(m_switchData[iBlock]);
          data.erase(data.begin());
          // reading old r0
          double r0;
          if (!Tools::parse(data, "R_0", r0)) {
            log << "Error with the R_0 parameter of the switching function\n";
          }
          std::string sR0; 
          Tools::convert(r0, sR0);
          // computing new r0
          r0 *= m_volumeFactor;
          auto pos = m_switchData[iBlock].find("R_0");
          //m_switchData[iBlock].replace(pos + 4, std::to_string(r0).size(), std::to_string(r0));
          m_switchData[iBlock].replace(pos + 4, sR0.size(), std::to_string(r0));
        }
        m_switchFunc[iBlock].set(m_switchData[iBlock], errors);
        if (errors.length() != 0){
          error("problem reading SWITCH" + std::to_string(iBlock + 1) + " keyword : " + errors );
        }
        m_r00[iBlock] = m_switchFunc[iBlock].get_r0();
        log << "  Swf: " << iBlock + 1 << "  r0 = " << (m_switchFunc[iBlock].description()).c_str() << " \n";
      }
      
      //Transform and sort
      log << "Building Reference PIV2 Vector \n";
      log << "  PIV2  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
      double df = 0.;
      for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        for (unsigned i = 0; i < m_refPIV2[iRef][iBlock].size(); i++) {
          m_refPIV2[iRef][iBlock][i] = m_switchFunc[iBlock].calculate(m_refPIV2[iRef][iBlock][i], df);
        }
        if (m_doSort[iBlock]) {
          std::sort(m_refPIV2[iRef][iBlock].begin(), m_refPIV2[iRef][iBlock].end());
        }
        unsigned lmt0 = 0, lmt1 = 0;
        for (unsigned iAtm = 0; iAtm < m_refPIV2[iRef][iBlock].size(); iAtm++) {
          if (m_refPIV2[iRef][iBlock][iAtm] > 0.9) { lmt1++; }
          if (m_refPIV2[iRef][iBlock][iAtm] < 0.1) { lmt0++; }
        }
        log.printf("        |%10i|%15i|%15i|%15i|\n", iBlock, m_refPIV2[iRef][iBlock].size(), lmt0, lmt1);
      } // for each block
    } // for each reference file

    // we compute lambda
    double distance = 0.;
    for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      unsigned size;
      unsigned last = m_numberReferences - 1;
      if (m_refPIV2[0][iBlock].size() > m_refPIV2[last][iBlock].size()) {
        size = m_refPIV2[last][iBlock].size();
      } else {
        size = m_refPIV2[0][iBlock].size();
      }
      for (unsigned iAtm = 0; iAtm < size; iAtm++) {
        //double coord = m_refPIV2[m_numberReferences - 1][iBlock][i] - m_refPIV2[0][iBlock][i];
        double coord = m_refPIV2[last][iBlock][iAtm] 
                       - m_refPIV2[0][iBlock][m_refPIV2[0][iBlock].size() 
                                              - m_refPIV2[last][iBlock].size() + iAtm];
        distance += m_blockScaling[iBlock] * coord * coord;
      }
    }
    m_lambda = 2.3 / distance;
    log << "lambda=" << m_lambda << " d_1n=" << distance << "\n";

    m_refIsIni = true;
    log << "\n";
  } // building of the reference PIV2 

  // Do the sorting with a stride defined by updatePIV2 
  if (getStep() % m_updatePIV2 == 0 || m_computeDerivatives) {
    if (m_computeDerivatives) {
      log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV2 \n";
    }

    // build COMs from positions if requested
    if (m_centerOfMass) {
      if (m_pbc) {
        makeWhole();
      }
      for(unsigned iCOM = 0; iCOM < m_positionCOM.size(); iCOM++) {
        m_positionCOM[iCOM][0] = 0.;
        m_positionCOM[iCOM][1] = 0.;
        m_positionCOM[iCOM][2] = 0.;
        for(auto& atom : m_neighborlistCOM[iCOM]->getFullAtomList()) {
          m_positionCOM[iCOM] += m_massFactor[atom.index()] * getPosition(atom.index());
        }
      }
    }

    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (m_doNeighbor) {
      // for first step list = actual position
      if (!m_neighborIsIni) {
        for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
          for (unsigned iAtm = 0; iAtm < m_neighborList[iBlock]->getFullAtomList().size(); iAtm++) {
            Vector position;
            if (m_centerOfMass) {
              position = m_positionCOM[iAtm];
            } else {
              position = getPosition(m_neighborList[iBlock]->getFullAtomList()[iAtm].index());
            }
            m_prevPosition[iBlock].push_back(position);
          }
          m_neighborList[iBlock]->update(m_prevPosition[iBlock]);
        }
        m_neighborIsIni = true;
      }

      // Decide whether to update lists based on atom displacement, every stride
      if (getStep() % m_neighborListAll->getStride() == 0) {
        bool doUpdate = false;
        for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
          for (unsigned iAtm = 0; iAtm < m_neighborList[iBlock]->getFullAtomList().size(); iAtm++) {
            Vector position;
            if (m_centerOfMass) {
              position = m_positionCOM[iAtm];
            } else {
              position = getPosition(m_neighborList[iBlock]->getFullAtomList()[iAtm].index());
            }
            if (pbcDistance(position, m_prevPosition[iBlock][iAtm]).modulo2()
                >= m_neighborListSkin[iBlock] * m_neighborListSkin[iBlock]) {
              doUpdate = true;
              m_prevPosition[iBlock][iAtm] = position;
            }
          }
          // update positions if needed
          if (doUpdate == true) {
            m_neighborList[iBlock]->update(m_prevPosition[iBlock]);
            if (getStep() % 10000 == 0) {
              if (iBlock == 0) 
                log << "  Step " << getStep() << " neighbour lists updated for block ";
              log << iBlock + 1 << ": " << m_neighborList[iBlock]->size() << "; ";
            }
          }
        } // for each block
        if (getStep() % 10000 == 0) log << "\n";
      } // if step % nlStride == 0
    } // if do neighbor 

    Vector pairDist;
    // Build "neigborlist" PIV2 blocks
    for(unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
      if (m_doSort[iBlock]) {
        // Vectors collecting occupancies: orderVec one rank, orderVecAll all ranks
        std::vector<int> orderVec(m_nPrecision, 0);
        cPIV2[iBlock].resize(0);
        atmI0[iBlock].resize(0);
        atmI1[iBlock].resize(0);

        // Building distances for the PIV2 vector at time t
        if (m_timer) stopwatch.start("1 Build cPIV2");

        // If we have N cores and Na atoms, we need (Na/N + 1) process by cores
        unsigned procByCore = unsigned(m_neighborList[iBlock]->size() / stride + 1);
        for (unsigned iAtm = rank * procByCore; (iAtm < ((rank + 1) * procByCore)) && (iAtm < m_neighborList[iBlock]->size()); iAtm++) { /*
        for(unsigned iAtm = rank; iAtm < m_neighborList[iBlock]->size(); iAtm += stride) {*/
          unsigned i0 = (m_neighborList[iBlock]->getClosePairAtomNumber(iAtm).first).index();
          unsigned i1 = (m_neighborList[iBlock]->getClosePairAtomNumber(iAtm).second).index();
          Vector position0, position1;
          if (m_centerOfMass) {
            position0 = m_positionCOM[i0];
            position1 = m_positionCOM[i1];
          } else {
            position0 = getPosition(i0);
            position1 = getPosition(i1);
          }
          pairDist = distanceAB(position0, position1);
          double df = 0.;
          //Transforming distances with the Switching function + real to integer transformation
          if (pairDist.modulo2() < m_cutOff[iBlock]*m_cutOff[iBlock]) {
            int vecInt = int(m_switchFunc[iBlock].calculate(pairDist.modulo() * m_volumeFactor, df)
                             * double(m_nPrecision - 1) + 0.5);
            //Integer transformed distance values as index of the Ordering Vector orderVec
            orderVec[vecInt] += 1;
            //Keeps track of atom indices for force and virial calculations
            atmPrecI0[vecInt].push_back(i0);
            atmPrecI1[vecInt].push_back(i1);
          } else {
            int vecInt = 0;
            orderVec[vecInt] += 1;
            //Keeps track of atom indices for force and virial calculations
            atmPrecI0[vecInt].push_back(i0);
            atmPrecI1[vecInt].push_back(i1);
          }
        }
         
        if (m_timer) stopwatch.stop("1 Build cPIV2");
        if (m_timer) stopwatch.start("2 Sort cPIV2");

        if (!m_serial) {
          // Vectors keeping track of the dimension and the starting-position 
          // of the rank-specific pair vector in the big pair vector.
          std::vector<int> vecDimension(stride, 0);
          std::vector<int> vecPos(stride, 0);
          // Vectors collecting occupancies: orderVec one rank, orderVecAll all ranks
          std::vector<int> orderVecAll(stride * m_nPrecision);
          // Big vectors containing all Atom indexes for every occupancy 
          // (atmI0O(Nprec,n) and atmI1O(Nprec,n) matrices in one vector)
          std::vector<int> atmI0F;
          std::vector<int> atmI1F;
          // Vector used to reconstruct arrays
          std::vector<unsigned> counter(stride, 0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for(unsigned i = 1; i < m_nPrecision; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            atmI0F.insert(atmI0F.end(), atmPrecI0[i].begin(), atmPrecI0[i].end());
            atmI1F.insert(atmI1F.end(), atmPrecI1[i].begin(), atmPrecI1[i].end());
            atmPrecI0[i].resize(0);
            atmPrecI1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV2 block
          atmPrecI0[0].resize(0);
          atmPrecI1[0].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          orderVec[0] = 0;
          orderVec[m_nPrecision - 1] = 0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim = atmI0F.size();
          // vecDimension and vecPos keep track of the dimension and the starting-position 
          // of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim, 1, &vecDimension[0], 1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int finalDimension = 0;
          for(unsigned i = 1; i < stride; i++) {
            vecPos[i] = vecPos[i-1] + vecDimension[i-1];
            finalDimension += vecDimension[i];
          }
          finalDimension += vecDimension[0];
          
          if (getStep() == 0)
            log << "finalDimension=" << finalDimension << " " << dim << "\n";

          // build big vectors for atom pairs on all ranks for all ranks
          std::vector<int> atmI0FinalAll(finalDimension);
          std::vector<int> atmI1FinalAll(finalDimension);
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          // Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV2
          comm.Allgather(&orderVec[0], m_nPrecision, &orderVecAll[0], m_nPrecision);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&atmI0F[0], atmI0F.size(), &atmI0FinalAll[0], &vecDimension[0], &vecPos[0]);
          comm.Allgatherv(&atmI1F[0], atmI1F.size(), &atmI1FinalAll[0], &vecDimension[0], &vecPos[0]);

          if (m_timer) stopwatch.stop("2 Sort cPIV2");
          if (m_timer) stopwatch.start("3 Reconstruct cPIV2");

          // Reconstruct the full vectors from collections of Allgathered parts 
          // This is a tricky serial step, to assemble toghether PIV2 and 
          // atom-pair info from head-tail big vectors. Loop before on l and 
          // then on i would be better but the allgather should be modified
          for (unsigned i = 1; i < m_nPrecision; i++) {
            for (unsigned l = 0; l < stride; l++) {
              for (unsigned m = 0; m < orderVecAll[i + l * m_nPrecision]; m++) {
                cPIV2[iBlock].push_back( double(i) / double(m_nPrecision - 1) );
                atmI0[iBlock].push_back( atmI0FinalAll[counter[l] + vecPos[l]] );
                atmI1[iBlock].push_back( atmI1FinalAll[counter[l] + vecPos[l]] );
                counter[l] += 1;
              } // loop on the number of head-to-tail pieces
            } // loop on the ranks
          } // loop on the ordering vector excluding zero (i = 1)

          if (m_timer) stopwatch.stop("3 Reconstruct cPIV2");

        } else {
          for(unsigned i = 1; i < m_nPrecision; i++) {
            for(unsigned m = 0; m < orderVec[i]; m++) {
              cPIV2[iBlock].push_back( double(i) / double(m_nPrecision - 1) );
              atmI0[iBlock].push_back(atmPrecI0[i][m]);
              atmI1[iBlock].push_back(atmPrecI1[i][m]);
            }
          }
        } // if serial or parallel
      } // if we sort the PIV
    } // for each block
  } // if step % stride == 0

  Vector distance;
  double dfunc = 0.;

  // This test may be run by specifying the TEST keyword as input, it pritnts referencePIV2 and cPIV2 and quits
  if (m_doTest) {
    unsigned limit = 0;
    for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
      for (unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        if (m_doSort[iBlock]) {
          limit = cPIV2[iBlock].size();
        } else {
          limit = m_refPIV2[iRef][iBlock].size();
        }
        log.printf("PIV2 Block:  %6i %12s %6i \n", iBlock, "      Size:", limit);
        log.printf("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV2   ",
                   "    r-PIV2   ","   i-j distance vector       ");
        for(unsigned i = 0; i < limit; i++) {
          unsigned i0 = 0;
          unsigned i1 = 0;
          if (m_doSort[iBlock]) {
            i0 = atmI0[iBlock][i];
            i1 = atmI1[iBlock][i];
          } else {
            i0 = (m_neighborList[iBlock]->getClosePairAtomNumber(i).first).index();
            i1 = (m_neighborList[iBlock]->getClosePairAtomNumber(i).second).index();
          }
          Vector position0,position1;
          if (m_centerOfMass) {
            position0 = m_positionCOM[i0];
            position1 = m_positionCOM[i1];
          } else {
            position0 = getPosition(i0);
            position1 = getPosition(i1);
          }
          distance = distanceAB(position0, position1);
          dfunc = 0.;
          double cP, rP;
          if (m_doSort[iBlock]) {
            cP = cPIV2[iBlock][i];
            rP = m_refPIV2[iRef][iBlock][m_refPIV2[iRef][iBlock].size() - cPIV2[iBlock].size() + i];
          } else {
            cP = m_switchFunc[iBlock].calculate(distance.modulo() * m_volumeFactor, dfunc);
            rP = m_refPIV2[iRef][iBlock][i];
          }
          log.printf("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
        }
      }
      log.printf("This was a test, now exit \n");
      exit();
    }
  }
  
  if (m_timer) stopwatch.start("4 Build For Derivatives");

  if (getStep() % m_updatePIV2 == 0) {
    // set to zero PIV2distance, derivatives and virial when they are calculated
    for(unsigned j = 0; j < m_deriv.size(); j++) {
      for(unsigned k = 0; k < 3; k++) {
        m_deriv[j][k] = 0.;
        if (j < 3) m_virial[j][k] = 0.;
      }
    }
    // compute PIV-PIV distance and derivatives for each reference
    for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
      m_distancePIV2[iRef] = 0.;
      // Re-compute atomic distances for derivatives and compute PIV2-PIV2 distance
      for(unsigned iBlock = 0; iBlock < m_numberBlocks; iBlock++) {
        unsigned limit = 0;
        if (m_doSort[iBlock]) {
          limit = cPIV2[iBlock].size();
        } else {
          limit = m_refPIV2[iRef][iBlock].size();
        }
        unsigned procByCore = unsigned(limit / stride + 1);
        for (unsigned i = rank * procByCore; (i < ((rank + 1) * procByCore)) && (i < limit); i++) {
          unsigned i0 = 0;
          unsigned i1 = 0;
          if (m_doSort[iBlock]) {
            i0 = atmI0[iBlock][i];
            i1 = atmI1[iBlock][i];
          } else {
            i0 = (m_neighborList[iBlock]->getClosePairAtomNumber(i).first).index();
            i1 = (m_neighborList[iBlock]->getClosePairAtomNumber(i).second).index();
          }
          Vector position0, position1;
          if (m_centerOfMass) {
            position0 = m_positionCOM[i0];
            position1 = m_positionCOM[i1];
          } else {
            position0 = getPosition(i0);
            position1 = getPosition(i1);
          }
          distance = distanceAB(position0, position1);
          dfunc = 0.;
          // this is needed for dfunc and dervatives
          double dm = distance.modulo();
          double tPIV2 = m_switchFunc[iBlock].calculate(dm * m_volumeFactor, dfunc);
          
          // PIV2 distance
          double coord = 0.;
          if (!m_doSort[iBlock] || m_computeDerivatives) {
            coord = tPIV2 - m_refPIV2[iRef][iBlock][i];
          } else {
            coord = cPIV2[iBlock][i] 
                    - m_refPIV2[iRef][iBlock][m_refPIV2[iRef][iBlock].size() 
                    - cPIV2[iBlock].size() + i];
          }
          // Calculate derivatives, virial, and variable = sum_j (scaling[j] *(cPIV2-rPIV2)_j^2)
          // WARNING: dfunc=dswf/(m_volumeFactor * dm)  (this may change in future Plumed versions)
          double tmp = 2. * m_blockScaling[iBlock] * coord
                       * m_volumeFactor * m_volumeFactor * dfunc;
          Vector tempDeriv = tmp * distance;
          // 0.5 * (x_i - x_k) * f_ik         (force on atom k due to atom i)
          if (m_centerOfMass) {
            Vector dist;
            for(auto& atom0 : m_neighborlistCOM[i0]->getFullAtomList()) {
              unsigned x0 = atom0.index();
              m_deriv[x0] -= tempDeriv * m_massFactor[x0];
              for(unsigned l = 0; l < 3; l++) {
                dist[l] = 0.;
              }
              Vector position0 = getPosition(x0);
              for(auto& atom1 : m_neighborlistCOM[i0]->getFullAtomList()) {
                dist += distanceAB(position0, getPosition(atom1.index()) );
              }
              for(auto& atom1 : m_neighborlistCOM[i1]->getFullAtomList()) {
                dist += distanceAB(position0, getPosition(atom1.index()) );
              }
              m_virial -= 0.25 * m_massFactor[x0] * Tensor(dist, tempDeriv);
            } // loop on first atom of each pair
            for(auto& atom1 : m_neighborlistCOM[i1]->getFullAtomList()) {
              unsigned x1 = atom1.index();
              m_deriv[x1] += tempDeriv * m_massFactor[x1];
              for(unsigned l = 0; l < 3; l++) {
                dist[l] = 0.;
              }
              Vector position1 = getPosition(x1);
              for(auto& atom0 : m_neighborlistCOM[i1]->getFullAtomList()) {
                dist += distanceAB(position1, getPosition(atom0.index()) );
              }
              for(auto& atom0 : m_neighborlistCOM[i0]->getFullAtomList()) {
                dist += distanceAB(position1, getPosition(atom0.index()) );
              }
              m_virial += 0.25 * m_massFactor[x1] * Tensor(dist, tempDeriv);
            } // loop on second atom of each pair
          } else {
            m_deriv[i0] -= tempDeriv;
            m_deriv[i1] += tempDeriv;
            m_virial    -= tmp * Tensor(distance, distance);
          } // if do center of mass
          if (m_doScaleVolume) {
            m_virial += 1./3. * tmp * dm * dm * Tensor::identity();
          }
          m_distancePIV2[iRef] += m_blockScaling[iBlock] * coord * coord;
        } // loop on atoms-atoms pairs
      } // loop on block

      if (!m_serial) {
        comm.Barrier();
        comm.Sum(&m_distancePIV2[iRef], 1);
        if (!m_deriv.empty()) {
          comm.Sum(&m_deriv[0][0], 3 * m_deriv.size());
        }
        comm.Sum(&m_virial[0][0], 9);
      }
    } // loop on the number of references
  } // if update_piv
  
  //Timing
  if (m_timer) stopwatch.stop("4 Build For Derivatives");
  if (m_timer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log << stopwatch;
  }

  // Update derivatives, virial, and variable (PIV2-distance^2)
  for(unsigned i = 0; i < m_deriv.size(); ++i) {
    setAtomsDerivatives(i, m_deriv[i]);
  }
  setBoxDerivatives(m_virial);
  for (unsigned iRef = 0; iRef < m_numberReferences; iRef++) {
    Value* pValDistance = getPntrToComponent("d" + std::to_string(iRef + 1));
    pValDistance->set(m_distancePIV2[iRef]);
  }
  Value* pValLambda = getPntrToComponent("lambda");
  pValLambda->set(m_lambda);
} // end of calculate

} // close namespace piv
} // close namespace PLMD
