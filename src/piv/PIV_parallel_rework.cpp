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

//+PLUMEDOC COLVAR PIVbin
/*
Calculates the PIV-distance: the squared Cartesian distance between the PIV \cite gallet2013structural,pipolo2017navigating
associated to the configuration of the system during the dynamics and a reference configuration provided
as input (PDB file format).
PIV can be used together with \ref FUNCPATHMSD to define a path in the PIV space.
\par Examples

The following example calculates PIV-distances from three reference configurations in Ref1.pdb, Ref2.pdb and Ref3.pdb
and prints the results in a file named colvar.
Three atoms (PIVATOMS=3) with names (pdb file) A B and C are used to construct the PIV and all PIV blocks (AA, BB, CC, AB, AC, BC) are considered.
SFACTOR is a scaling factor that multiplies the contribution to the PIV-distance given by the single PIV block.
NLIST sets the use of neighbor lists for calculating atom-atom distances.
The SWITCH keyword specifies the perameters of the switching function that transforms atom-atom distances.
SORT=1 meand that the PIV block elements are sorted (SORT=0 no sorting.)
Values for SORT, SFACTOR and Neighborlist parameters have to be specified for each block.
The order is the following: AA,BB,CC,AB,AC,BC. If ONLYDIRECT (ONLYCROSS) is used the order is AA,BB,CC (AB,AC,BC).
The sorting operation within each PIV block is performed using the counting sort algorithm, PRECISION specifies the size of the counting array.
\plumedfile
PIVbin ...
LABEL=piv
REF_FILE1=Ref1.pdb
REF_FILE2=Ref2.pdb
REF_FILE3=Ref2.pdb
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
... PIV

PRINT ARG=piv.d1,piv.d2,piv.d3 FILE=colvar
\endplumedfile

WARNING:
Both the "CRYST" and "ATOM" lines of the PDB files must conform precisely to the official pdb format, including the width of each alphanumerical field:

\verbatim
CRYST1   31.028   36.957   23.143  89.93  92.31  89.99 P 1           1
ATOM      1  OW1 wate    1      15.630  19.750   1.520  1.00  0.00
\endverbatim

In each pdb frame, atoms must be numbered in the same order and with the same element symbol as in the input of the MD program.

The following example calculates the PIV-distances from two reference configurations Ref1.pdb and Ref2.pdb
and uses PIV-distances to define a Path Collective Variable (\ref FUNCPATHMSD) with only two references (Ref1.pdb and Ref2.pdb).
With the VOLUME keyword one scales the atom-atom distances by the cubic root of the ratio between the specified value and the box volume of the initial step of the trajectory file.

\plumedfile
PIVbin ...
LABEL=piv
VOLUME=12.15
REF_FILE1=Ref1.pdb
REF_FILE2=Ref2.pdb
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.5 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIVbin

p1: FUNCPATHMSD ARG=piv.d1,piv.d2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=piv.d1,piv.d2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIV please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class PIVbin : public Colvar
{
public:
  static void registerKeywords ( Keywords& keys );
  PIVbin (const ActionOptions&);
  ~PIVbin ();
  // active methods:
  virtual void calculate ();
  void checkFieldsAllowed () {}
private:
  Vector distanceAB (Vector A, Vector B);
  void setCutOff (bool cutOff);
private:
  ForwardDecl<Stopwatch> stopwatch_fwd;
  /// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch = *stopwatch_fwd;
  bool mPBC, mSerial, mTimer, mFirstStep; 
  bool mScaleVolume, mCross, mDirect, mDoNeighbor;
  bool mTest, mComputeDerivatives, mDoCom;
  unsigned mPrecision, mAtomTypes, mBlocks; //, mBlockAtomsSize, 
  unsigned mNumberReferences, mUpdateStride;
  double mVolumeFactor, mVolume0, mLambda;
  Tensor mVirial;
  NeighborList* mBlockAtomsAll;
  std::vector<bool> mDoSort;
  std::vector<double> mBlockFactor, m_r00;
  std::vector<double> mBlockAtomsSkin;
  std::vector<double> mCutOff;
  std::vector<double> mDistancePIV;
  std::vector<Vector> mDerivatives;
  std::vector<std::string> mSwitchData;
  std::vector<NeighborList*> mBlockAtoms;
  std::vector<SwitchingFunction> mSwitchFunc;
  std::vector<std::vector<Vector>> mPrevPosition;
  // first vector for reference, second for block, third for atoms
  std::vector<std::vector<std::vector<double>>> mRefPIV;
};

PLUMED_REGISTER_ACTION (PIVbin, "PIVbin")

void PIVbin::registerKeywords (Keywords& keys)
{
  Colvar::registerKeywords (keys);
  keys.add ("numbered", "SWITCH", "The switching functions parameter."
           "You should specify a Switching function for all PIV blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add ("numbered", "REF_FILE", "PDB file name that contains the i-th reference structure."
           "If you indicate n reference structure, you will get d_n PIV distances that can be used.");
  keys.add ("compulsory", "SORT", "Whether to sort or not the PIV block.");
  keys.add ("compulsory", "ATOMTYPES", "The atomtypes to use for PIV.");
  keys.add ("optional", "PRECISION", "the precision for approximating reals with integers in sorting.");
  keys.add ("optional", "SFACTOR", "Scale the PIV-distance by such block-specific factor");
  keys.add ("optional", "VOLUME", "Scale atom-atom distances by the cubic root of the cell volume."
           "The input volume is used to scale the R_0 value of the switching function. ");
  keys.add ("optional", "UPDATEPIV", "Frequency (timesteps) at which the PIV is updated.");
  keys.add ("optional", "NL_CUTOFF", "Neighbour lists cutoff.");
  keys.add ("optional", "NL_STRIDE", "Update neighbour lists every NL_STRIDE steps.");
  keys.add ("optional", "NL_SKIN", "The maximum atom displacement tolerated for the neighbor lists update.");
  keys.addFlag ("TEST", false, "Print the actual and reference PIV and exit");
  keys.addFlag ("ONLYCROSS", false, "Use only cross-terms (A-B, A-C, B-C, ...) in PIV");
  keys.addFlag ("ONLYDIRECT", false, "Use only direct-terms (A-A, B-B, C-C, ...) in PIV");
  keys.addFlag ("DERIVATIVES", false, "Activate the calculation of the PIV for every class"
               "(needed for numerical derivatives).");
  keys.addFlag ("NO_NLIST", false, "Don't use a neighbour list for distance calculations.");
  keys.addFlag ("SERIAL", false, "Perform the calculation in serial - for debug purpose");
  keys.addFlag ("TIMER", false, "Permorm timing analysis on heavy loops.");
  keys.addFlag ("NO_CUTOFF", false, "Don't use cut-off to reduce the computational cost of the PIV.");
  keys.reset_style ("SWITCH", "compulsory");
  componentsAreNotOptional (keys);
  // output
  keys.addOutputComponent ("lambda", "default", "Optimal lambda needed for the pathCV "
                           "It verifies that: lambda * D_{1-n} = 2.3");
  keys.addOutputComponent ("d1", "default", "PIV distance between the first "
                           "reference state and the current state");
  keys.addOutputComponent ("d2", "default", "PIV distance between the second "
                           "reference state and the current state");
  keys.addOutputComponent ("d3", "default", "PIV distance between the third "
                           "reference state and the current state");
  keys.addOutputComponent ("d4", "default", "PIV distance between the fourth "
                           "reference state and the current state");
  keys.addOutputComponent ("di", "default", "PIV distance between the i-th "
                           "reference state and the current state");
}

//---------------------------------------------------------
// CONSTRUCTOR
//---------------------------------------------------------

PIVbin::PIVbin (const ActionOptions&ao):
  PLUMED_COLVAR_INIT (ao),
  mPBC (true),
  mSerial (false),
  mTimer (false),
  mFirstStep (true),
  mScaleVolume (false),
  mCross (true),
  mDirect (true),
  mDoNeighbor (true),
  mTest (false),
  mComputeDerivatives (false),
  mDoCom (false),
  mPrecision (1000),
  mNumberReferences (0),
  mUpdateStride (1),
  mVolumeFactor (1.),
  mVolume0 (1.),
  mLambda (1.),
  mBlockAtoms (std::vector<NeighborList*> (1))
{
  log << "Starting PIV Constructor\n";
  bool onlyCross = false, onlyDirect = false, noPbc = !mPBC, noNeighbor = !mDoNeighbor, noCut = false;

  // parse all the mandatory inputs that are not vector
  parse ("VOLUME", mVolume0);
  // parse all the options
  parseFlag ("TEST", mTest);
  parseFlag ("NOPBC", noPbc);
  parseFlag ("TIMER", mTimer);
  parseFlag ("SERIAL", mSerial);
  parseFlag ("NO_NLIST", noNeighbor);
  parseFlag ("NO_CUTOFF", noCut);
  parseFlag ("ONLYCROSS", onlyCross);
  parseFlag ("ONLYDIRECT", onlyDirect);
  parseFlag ("DERIVATIVES", mComputeDerivatives);

  // parse the atom names 
  std::vector <std::string> atomTypes;
  parseVector ("ATOMTYPES", atomTypes);
  mAtomTypes = atomTypes.size();
  // Stride for which the PIV are computed
  if (keywords.exists ("UPDATEPIV")) {
    parse ("UPDATEPIV", mUpdateStride);
  }
  // Precision on the real-to-integer transformation for the sorting
  if (keywords.exists ("PRECISION")) {
    parse ("PRECISION", mPrecision);
  }
  log << "Precision N = " << mPrecision << "\n";
  if (mPrecision < 2) { 
    error ("Precision must be => 2");
  }
  // PBC
  mPBC = !noPbc;
  if (mPBC) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }
  // Serial or parallel
  if (mSerial) {
    log << "Serial PIV construction\n";
  } else     {
    log << "Parallel PIV construction\n";
  }
  // Derivatives
  if (mComputeDerivatives) {
    log << "Computing Derivatives\n";
  }
  // Timing
  if (mTimer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }
  // Volume Scaling
  if (mVolume0 > 0) {
    mScaleVolume = true;
  }
  // PIV direct and cross blocks
  if (onlyCross && onlyDirect) {
    error("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  if (onlyCross) {
    mDirect = false;
    log << "Using only CROSS-PIV blocks\n";
  } 
  if (onlyDirect) {
    mCross = false;
    log << "Using only DIRECT-PIV blocks\n";
  }
  mBlocks = 0;
  if (mCross) {
    mBlocks += unsigned(double(mAtomTypes * (mAtomTypes - 1)) / 2.);
  }
  if (mDirect) {
    mBlocks += unsigned(mAtomTypes);
  }

  // resizing all class vector according to mBlocks
  mDoSort.resize (mBlocks);
  mBlockFactor.resize (mBlocks);
  m_r00.resize (mBlocks);
  mSwitchData.resize (mBlocks); 
  mBlockAtomsSkin.resize (mBlocks);
  mPrevPosition.resize (mBlocks);
  mBlockAtoms.resize (mBlocks);

  // setting neighborlist parameters 
  std::vector<double> neighborListCut(mBlocks, 0.);
  std::vector<int> neighborListStride(mBlocks, 0);
  mDoNeighbor = !noNeighbor;
  if (mDoNeighbor) {
    parseVector ("NL_CUTOFF", neighborListCut);
    parseVector ("NL_STRIDE", neighborListStride);
    parseVector ("NL_SKIN", mBlockAtomsSkin);
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      if (neighborListCut[iBloc] <= 0.0) {
        error("NL_CUTOFF should be explicitly specified and positive");
      }
      if (neighborListStride[iBloc] <= 0) {
        error("NL_STRIDE should be explicitly specified and positive");
      }
      if (mBlockAtomsSkin[iBloc] <= 0.) {
        error("NL_SKIN should be explicitly specified and positive");
      }
      neighborListCut[iBloc] = neighborListCut[iBloc] + mBlockAtomsSkin[iBloc];
      log << "For Block " << iBloc + 1 
          << ", neighbor list cut-off=" << neighborListCut[iBloc]
          << ", stride=" << neighborListStride[iBloc]
          << ", shell=" << mBlockAtomsSkin[iBloc] << "\n";
    }
  }

  // Sorting
  std::vector<unsigned> yesNoSort(mBlocks);
  parseVector ("SORT", yesNoSort);
  for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
    if (yesNoSort[iBloc] == 0 || mComputeDerivatives) {
      mDoSort[iBloc] = false;
      log << "Not sorting block " << iBloc + 1 << ". ";
    } else {
      mDoSort[iBloc] = true;
      log << "Sorting block " << iBloc + 1 << ". ";
    }
  }
  log << "\n";

  // PIV scaled option
  for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
    mBlockFactor[iBloc] = 1.;
  }
  if (keywords.exists("SFACTOR")) {
    parseVector ("SFACTOR", mBlockFactor);
  }

  // read parameters and set-up switching functions here only if computing derivatives
  for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
    if ( !parseNumbered("SWITCH", iBloc + 1, mSwitchData[iBloc]) ){
      log << "Problem while reading the switching function parameters.\n";
      break;
    }
  }
  if (mComputeDerivatives) {
    log << "Switching Function Parameters \n";
    mSwitchFunc.resize(mBlocks);
    std::string errors;
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      // std::string num;
      // Tools::convert(iBloc + 1, num);
      mSwitchFunc[iBloc].set(mSwitchData[iBloc], errors);
      if (errors.length() != 0){
        error("problem reading switch" + std::to_string(iBloc + 1) + " keyword : " + errors );
      }
      m_r00[iBloc] = mSwitchFunc[iBloc].get_r0();
      log << "  Swf: " << iBloc + 1 << "  r0=" << (mSwitchFunc[iBloc].description()).c_str() << "\n";
    }
  }

  // set the cut-off
  setCutOff (!noCut);

  // Reference PDB file 
  for (unsigned iRef = 0; ; iRef++) {
    std::string referenceFile; 
    parseNumbered("REF_FILE", iRef + 1, referenceFile);
    if (referenceFile.empty()) {
      break;
    }
    mNumberReferences++;
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
 
    // Build Atom lists of AtomNumbers (this might be done in PBC.cpp)
    std::vector<std::vector<AtomNumber>> pairList(mAtomTypes);
    // number of atoms in the PDB cell
    unsigned neighborListSize = myPDB.getAtomNumbers().size();
    log << "Atoms numbers " << myPDB.getAtomNumbers().size() << "\n";

    // Build residueIndex
    std::vector<unsigned> residueIndex;
    for (unsigned iAtm = 0; iAtm < myPDB.getAtomNumbers().size(); iAtm++) {
      unsigned residueNumber = myPDB.getResidueNumber(myPDB.getAtomNumbers()[iAtm]);
      if (std::find(residueIndex.begin(), residueIndex.end(), residueNumber)
          == residueIndex.end()) {
        residueIndex.push_back (residueNumber);
      }
    }
    unsigned nResidues = residueIndex.size();
 
    // indexList is the index of atoms used in neighborlists 
    unsigned indexListSize = neighborListSize;
    std::vector<unsigned> indexList(indexListSize);
    log << "Total Atoms: " << mAtomTypes * nResidues << " \n";

    // Build lists of Atoms for neighbor lists
    for (unsigned iType = 0; iType < mAtomTypes; iType++) {
      //unsigned outputID = 0;
      std::fill (indexList.begin(), indexList.end(), 0);
      // This builds lists for neighbor lists
      for (auto& atom : myPDB.getAtomNumbers()) {
        std::string atomName;
        unsigned pointIndex = 0;
        atomName = myPDB.getAtomName(atom);
        pointIndex = atom.index();
        if (atomName == atomTypes[iType]) {
          if (indexList[pointIndex] == 0) {
            // adding the atomnumber to the atom list for pairs
            pairList[iType].push_back (atom);
            indexList[pointIndex] = atom.index() + 1;
            //outputID = pointIndex;
          }
        }
      }
      // Output Lists
      /*
      log << "  Groups of type  " << iType << ": " << pairList[iType].size() << " \n";
      std::string groupName;
      unsigned groupSize;
      groupName = myPDB.getAtomName(comAtoms[indexList[outputID] - 1][0]);
      groupSize = 1;
      log.printf ("    %6s %3s %13s %10i %6s\n", "type  ", groupName.c_str(),
                 "   containing ", groupSize," atoms");
      */
    }
 
    // This is to build the list with all the atoms
    std::vector<AtomNumber> listAllAtom;
    for (auto& atom : myPDB.getAtomNumbers()) {
      listAllAtom.push_back (atom);
    }
 
    if (mDoNeighbor) {
      log << "Creating Neighbor Lists \n";
      for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
        log << "For Block " << iBloc + 1 
            << ", neighbor list cut-off=" << neighborListCut[iBloc]
            << ", stride=" << neighborListStride[iBloc]
            << ", shell=" << mBlockAtomsSkin[iBloc] << "\n";
      }
      // WARNING: is neighborListCut meaningful here?
      mBlockAtomsAll= new NeighborList(listAllAtom, mPBC, getPbc(),
                             neighborListCut[0], neighborListStride[0]);
      unsigned ncnt = 0;
      // Direct blocks AA, BB, CC, ...
      if (mDirect) {
        log << "Number of blocks: " << mBlockAtoms.size() << ", number of atom"
            << "types: " << mAtomTypes << "\n";
        for (unsigned iType = 0; iType < mAtomTypes; iType++) {
          mBlockAtoms[ncnt] = new NeighborList(pairList[iType], mPBC, getPbc(),
                                     neighborListCut[iType], neighborListStride[iType]);
          ncnt += 1;
        }
      }
      // Cross blocks AB, AC, BC, ...
      if (mCross) {
        log << "Number of blocks: " << mBlockAtoms.size() << ", number of atom"
            << "types: " << mAtomTypes << "\n";
        for (unsigned iType1 = 0; iType1 < mAtomTypes; iType1++) {
          for (unsigned iType2 = iType1 + 1; iType2 < mAtomTypes; iType2++) {
            mBlockAtoms[ncnt] = new NeighborList(pairList[iType1], pairList[iType2],
                                       false, mPBC, getPbc(),
                                       neighborListCut[ncnt], neighborListStride[ncnt]);
            ncnt += 1;
          }
        }
      }
    } else {
      log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
      mBlockAtomsAll = new NeighborList(listAllAtom, mPBC, getPbc());
      for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
        mBlockAtoms[iBloc] = new NeighborList(pairList[iBloc], pairList[iBloc],
                                                  true, mPBC, getPbc());
      }
    }
    // Output neighborlist
    log << "Total Nlists: " << mBlocks << " \n";
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      log << "  list " << iBloc + 1 << "   size " << mBlockAtoms[iBloc]->size() << " \n";
    }

    //build box vectors and correct for pbc
    log << "Building the box from PDB data ... \n";
    Tensor box = myPDB.getBoxVec();
    log << "  Done! A,B,C vectors in Cartesian space:  \n";
    log.printf ("  A:  %12.6f%12.6f%12.6f\n", box[0][0], box[0][1], box[0][2]);
    log.printf ("  B:  %12.6f%12.6f%12.6f\n", box[1][0], box[1][1], box[1][2]);
    log.printf ("  C:  %12.6f%12.6f%12.6f\n", box[2][0], box[2][1], box[2][2]);
    log << "Changing the PBC according to the new box \n";
    Pbc myPbc;
    myPbc.setBox(box);
    log << "The box volume is " << myPbc.getBox().determinant() << " \n";

    //Compute scaling factor
    if (mScaleVolume) {
      mVolumeFactor = cbrt (mVolume0 / myPbc.getBox().determinant());
      log << "Scaling atom distances by  " << mVolumeFactor << " \n";
    } else {
      log << "Using unscaled atom distances \n";
    }
    
    // build the rPIV distances (transformation and sorting is done afterwards)
    std::vector<std::vector<double>> refPIV;
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      std::vector<double> blockRefPIV;
      for (unsigned iNl = 0; iNl < mBlockAtoms[iBloc]->size(); iNl++) {
        unsigned i0 = (mBlockAtoms[iBloc]->getClosePairAtomNumber(iNl).first).index();
        unsigned i1 = (mBlockAtoms[iBloc]->getClosePairAtomNumber(iNl).second).index();

        auto& positions = myPDB.getPositions();
        Vector pairDist;
        if (mPBC) {
          pairDist = myPbc.distance (positions[i0], positions[i1]);
        } else {
          pairDist = delta (positions[i0], positions[i1]);
        }
        double df = 0.;
        // Transformation and sorting done at the first timestep to solve the r0 definition issue
        //if (pairDist.modulo2() < mCutOff[iBloc] * mCutOff[iBloc]) {
        if (mComputeDerivatives) {
          blockRefPIV.push_back (mSwitchFunc[iBloc].calculate (pairDist.modulo() * mVolumeFactor, df) );
        } else {
          blockRefPIV.push_back (pairDist.modulo () * mVolumeFactor);
        }
        //}
      }
      refPIV.push_back (blockRefPIV);
    }
    mRefPIV.push_back (refPIV);
    
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      log << "reference PIV block " << iBloc + 1 << " has size: " 
          << mRefPIV[iRef][iBloc].size() << " over a total of "
           << mBlockAtoms[iBloc]->size() << " atoms-atoms pair\n";
      if (mComputeDerivatives) {
        if (mDoSort[iBloc]) {
          std::sort(mRefPIV[iRef][iBloc].begin(), mRefPIV[iRef][iBloc].end());
        }
        int lmt0 = 0;
        int lmt1 = 0;
        for (unsigned iAtm = 0; iAtm < mRefPIV[iRef][iBloc].size(); iAtm++) {
          if (mRefPIV[iRef][iBloc][iAtm] > 0.9) { lmt1++; }
          if (mRefPIV[iRef][iBloc][iAtm] < 0.1) { lmt0++; }
        }
        if (mComputeDerivatives && iBloc == 0) {
          log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
        }
        log.printf ("       |%10i|%15i|%15i|%15i|\n", iBloc, mRefPIV[iRef][iBloc].size(), lmt0, lmt1);
      } // if we compute derivatives
    } // loop over the number of blocks
  } // loop over the number of references states
  log << "\n";

  mDistancePIV.resize (mRefPIV.size());

  checkRead();

  // add N distance for the N reference states
  for (unsigned iRef = 1; iRef <= mRefPIV.size(); iRef++) {
    addComponentWithDerivatives("d" + std::to_string(iRef));
    componentIsNotPeriodic("d" + std::to_string(iRef));
  }
  // add the lambda component for the path collective variables
  addComponent("lambda");
  componentIsNotPeriodic("lambda");
  requestAtoms(mBlockAtomsAll->getFullAtomList());

  // set size of derivative vector
  mDerivatives.resize(getNumberOfAtoms()); 
}

//---------------------------------------------------------
// DESTRUCTOR
//---------------------------------------------------------

PIVbin::~PIVbin ()
{
  for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
    delete mBlockAtoms[iBloc];
  }
  delete mBlockAtomsAll;
}

//---------------------------------------------------------
// CUT-OFF
//---------------------------------------------------------

void PIVbin::setCutOff (bool doCutOff) 
{
  // just resolves the equation swf(x) = 1 / nPrecision for all switching functions
  for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++){
    SwitchingFunction switchFunc;
    std::string errors;
    switchFunc.set(mSwitchData[iBloc], errors);
    if (errors.length() != 0){
      error("problem reading switch" + std::to_string(iBloc + 1) + " keyword : " + errors );
    }

    log << "For block " << iBloc + 1;
    if (doCutOff) {
      unsigned maxIteration = 0;
      double f, g, df;
      double espilon = 2. / double(mPrecision);
      
      double x0 =  switchFunc.get_r0() / 2.;
      f = switchFunc.calculate(x0, df) - epsilon;
      while (fabs(f) > (epsilon / 100.) && maxIteration < 10000)
      {
        f = switchFunc.calculate(x0, df) - espilon;
        g = (switchFunc.calculate(x0 + f, df) - epsilon) / f - 1;
        x0 = x0 - f / g;
        maxIteration++;
      }
      mCutOff.push_back (x0);
      log << " epsilon=" << espilon << " this correspond to a cut-off=";
    } else {
      double dmax = switchFunc.get_dmax();
      mCutOff.push_back (dmax);
      log << ", cut-off=";
    }
    log << mCutOff[iBloc] << "\n";
  }
}

//---------------------------------------------------------
// DISTANCE 
//---------------------------------------------------------

inline Vector PIVbin::distanceAB (Vector A, Vector B)
{
  if (mPBC) {
    return pbcDistance (A, B);
  } else {
    return delta (A, B);
  }
}

//---------------------------------------------------------
// CALCULATE
//---------------------------------------------------------

void PIVbin::calculate()
{
  std::vector<std::vector<double>> currentPIV (mBlocks);
  /*
  std::vector<std::vector<int>> atmI0 (mBlocks);
  std::vector<std::vector<int>> atmI1 (mBlocks);
  std::vector<std::vector<int>> atmPrecI0 (mPrecision);
  std::vector<std::vector<int>> atmPrecI1 (mPrecision);
  */
  std::vector<std::vector<std::vector<int>>> atmI0 (mBlocks);
  std::vector<std::vector<std::vector<int>>> atmI1 (mBlocks);
  unsigned stride = 1;
  unsigned rank = 0;

  if (!mSerial) {
    stride = comm.Get_size();
    rank = comm.Get_rank();
  }

  // Compute volume scaling factor
  if (mScaleVolume) {
    mVolumeFactor = cbrt (mVolume0 / getBox().determinant());
  }

  ///////////////////////////////////////////////////////////////////
  // Transform (and sort) the rPIV before starting the dynamics
  if (mFirstStep && !mComputeDerivatives) {
    // Set switching function parameters
    log << "Switching Function Parameters \n";
    mSwitchFunc.resize(mBlocks);
    std::string errors;
    // setting the PIV reference for each atom pair blocks
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      if (mScaleVolume) {
        auto data = Tools::getWords(mSwitchData[iBloc]);
        data.erase(data.begin());
        // reading old r0
        double r0;
        if (!Tools::parse (data, "R_0", r0)) {
          log << "Error with the R_0 parameter of the switching function\n";
        }
        std::string sR0; 
        Tools::convert(r0, sR0);
        // computing new r0
        r0 *= mVolumeFactor;
        auto pos = mSwitchData[iBloc].find("R_0");
        mSwitchData[iBloc].replace(pos + 4, sR0.size(), std::to_string(r0));
      }
      mSwitchFunc[iBloc].set(mSwitchData[iBloc], errors);
      if (errors.length() != 0){
        error("problem reading SWITCH" + std::to_string(iBloc + 1) + " keyword : " + errors );
      }
      m_r00[iBloc] = mSwitchFunc[iBloc].get_r0();
      log << "  Swf: " << iBloc + 1 << "  r0 = " << (mSwitchFunc[iBloc].description()).c_str() << " \n";
    }

    // compute the N reference PIV
    for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++) {
      log << "\n";
      log << "REFERENCE PDB # " << iRef + 1 << " \n";
      //Transform and sort
      log << "Building Reference PIV Vector \n";
      log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
      double df = 0.;
      for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
        for (unsigned i = 0; i < mRefPIV[iRef][iBloc].size(); i++) {
          mRefPIV[iRef][iBloc][i] = mSwitchFunc[iBloc].calculate (mRefPIV[iRef][iBloc][i], df);
        }
        if (mDoSort[iBloc]) {
          std::sort (mRefPIV[iRef][iBloc].begin(), mRefPIV[iRef][iBloc].end());
        }
        unsigned lmt0 = 0, lmt1 = 0;
        for (unsigned iAtm = 0; iAtm < mRefPIV[iRef][iBloc].size(); iAtm++) {
          if (mRefPIV[iRef][iBloc][iAtm] > 0.9) { lmt1++; }
          if (mRefPIV[iRef][iBloc][iAtm] < 0.1) { lmt0++; }
        }
        log.printf ("        |%10i|%15i|%15i|%15i|\n", iBloc, mRefPIV[iRef][iBloc].size(), lmt0, lmt1);
      } // for each block
    } // for each reference file

    // we compute lambda
    double distance = 0.;
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      unsigned size;
      unsigned last = mRefPIV.size() - 1;
      if (mRefPIV[0][iBloc].size() > mRefPIV[last][iBloc].size()) {
        size = mRefPIV[last][iBloc].size();
      } else {
        size = mRefPIV[0][iBloc].size();
      }
      for (unsigned iAtm = 0; iAtm < size; iAtm++) {
        double coord = mRefPIV[last][iBloc][iAtm] 
                       - mRefPIV[0][iBloc][mRefPIV[0][iBloc].size() 
                                              - mRefPIV[last][iBloc].size() + iAtm];
        distance += mBlockFactor[iBloc] * coord * coord;
      }
    }
    double lambda = 2.3 / distance;
    log << "lambda=" << lambda << " d_1n=" << distance << "\n";
    Value* pValLambda = getPntrToComponent ("lambda");
    pValLambda->set (lambda);
    log << "\n";
  } // building of the reference PIV 

  // Vectors collecting occupancies: orderVec one rank, orderVecAll all ranks
  std::vector<std::vector<unsigned>> orderVec (mBlocks);
  ///////////////////////////////////////////////////////////////////
  // Do the sorting with a stride defined by updatePIV 
  if (getStep() % mUpdateStride == 0 || mComputeDerivatives || mFirstStep) {
    if (mComputeDerivatives) {
      log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV \n";
    }

    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (mDoNeighbor) {
      // for first step list = actual position
      if (mFirstStep) {
        for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
          for (unsigned iAtm = 0; iAtm < mBlockAtoms[iBloc]->getFullAtomList().size(); iAtm++) {
            mPrevPosition[iBloc].push_back ( 
              getPosition (mBlockAtoms[iBloc]->getFullAtomList ()[iAtm].index())); 
          }
          mBlockAtoms[iBloc]->update (mPrevPosition[iBloc]);
        }
      }

      // Decide whether to update lists based on atom displacement, every stride
      if (getStep() % mBlockAtomsAll->getStride() == 0) {
        bool doUpdate = false;
        std::vector< std::vector<Vector>> updatedPos (mBlocks);
        for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
          for (unsigned iAtm = 0; iAtm < mBlockAtoms[iBloc]->getFullAtomList().size(); iAtm++) {
            auto position = getPosition(mBlockAtoms[iBloc]->getFullAtomList()[iAtm].index());
            if (pbcDistance(position, mPrevPosition[iBloc][iAtm]).modulo2()
                >= mBlockAtomsSkin[iBloc] * mBlockAtomsSkin[iBloc]) {
              doUpdate = true;
              // mPrevPosition[iBloc][iAtm] = position;
            }
            updatedPos[iBloc].push_back (position);
          }
          // update positions if needed
          if (doUpdate == true) {
            mPrevPosition = updatedPos;
            mBlockAtoms[iBloc]->update (mPrevPosition[iBloc]);
            if (getStep() % 50000 == 0) {
              if (iBloc == 0) 
                log << "  Step " << getStep() << " neighbour lists updated for block ";
              log << iBloc + 1 << ": " << mBlockAtoms[iBloc]->size() << "; ";
            }
          }
        } // for each block
        if (getStep() % 50000 == 0) log << "\n";
      } // if step % nlStride == 0
    } // if do neighbor 

    Vector pairDist, deriv;
    // Build "neigborlist" PIV blocks
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      orderVec[iBloc].resize (mPrecision, 0);
      if (mDoSort[iBloc]) {
        currentPIV[iBloc].resize(0);
        atmI0[iBloc].resize(mPrecision);
        atmI1[iBloc].resize(mPrecision);

        // Building distances for the PIV vector at time t
        if (mTimer) stopwatch.start("1 Build currentPIV");

        // If we have N cores and Na atoms, we need (Na/N + 1) process by cores
        unsigned procByCore = unsigned(mBlockAtoms[iBloc]->size() / stride + 1);
        for (unsigned iAtm = rank * procByCore; (iAtm < ((rank + 1) * procByCore)) && (iAtm < mBlockAtoms[iBloc]->size()); iAtm++) { /*
        for (unsigned iAtm = rank; iAtm < mBlockAtoms[iBloc]->size(); iAtm += stride) {*/
          unsigned i0 = (mBlockAtoms[iBloc]->getClosePairAtomNumber(iAtm).first).index();
          unsigned i1 = (mBlockAtoms[iBloc]->getClosePairAtomNumber(iAtm).second).index();
          pairDist = distanceAB (getPosition (i0), getPosition (i1));

          double df = 0.;
          //Transforming distances with the Switching function + real to integer transformation
          if (pairDist.modulo2() < mCutOff[iBloc] * mCutOff[iBloc]) {
            int vecInt = int (mSwitchFunc[iBloc].calculate (pairDist.modulo() * mVolumeFactor, df)
                             * double(mPrecision - 1) + 0.5);
            //Integer transformed distance values as index of the Ordering Vector orderVec
            orderVec[iBloc][vecInt] += 1;
            //Keeps track of atom indices for force and virial calculations
            atmI0[iBloc][vecInt].push_back (i0);
            atmI1[iBloc][vecInt].push_back (i1);
          } else {
            int vecInt = 0;
            orderVec[iBloc][vecInt] += 1;
            //Keeps track of atom indices for force and virial calculations
            atmI0[iBloc][vecInt].push_back (i0);
            atmI1[iBloc][vecInt].push_back (i1);
          }
        }
         
        if (mTimer) stopwatch.stop("1 Build currentPIV");
        if (mTimer) stopwatch.start("2 Sort currentPIV");

        if (!mSerial) {
          // Vectors keeping track of the dimension and the starting-position 
          // of the rank-specific pair vector in the big pair vector.
          std::vector<int> vecDimension(stride, 0);
          std::vector<int> vecPos(stride, 0);

          // Avoid passing the zeros (i=1) for atom indices
          // orderVec[iBloc][0] = 0;
          // orderVec[iBloc][mPrecision - 1] = 0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          // Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV
          // comm.Allgather(&orderVec[0], mPrecision, &orderVecAll[0], mPrecision);
          comm.Sum (&orderVec[iBloc][0], mPrecision);
          if (mTimer) stopwatch.stop("2 Sort currentPIV");
        } /*
         else {
          for (unsigned i = 1; i < mPrecision; i++) {
            for (unsigned m = 0; m < orderVec[i]; m++) {
              currentPIV[iBloc].push_back ( double(i) / double(mPrecision - 1) );
              atmI0[iBloc].push_back (atmPrecI0[i][m]);
              atmI1[iBloc].push_back (atmPrecI1[i][m]);
            }
          }
        } // if serial or parallel
        */
      } // if we sort the PIV
    } // for each block
  } // if step % stride == 0

  Vector distance;
  double dfunc = 0.;

  ///////////////////////////////////////////////////////////////////
  // This test may be run by specifying the TEST keyword as input, 
  // it pritnts referencePIV and currentPIV and quits
  if (mTest) {
    unsigned limit = 0;
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      log.printf ("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV   ",
                 "    r-PIV   ","   i-j distance vector       ");
      std::vector<unsigned> counter (mRefPIV.size(), 0);
      for (unsigned iPrec = 1; iPrec < mPrecision; iPrec++) {
        if (atmI0[iBloc][iPrec].size() > 0) {
          unsigned i0 = atmI0[iBloc][iPrec][0];
          unsigned i1 = atmI1[iBloc][iPrec][0];
          distance = distanceAB (getPosition (i0), getPosition (i1));
          double cPIV = double (iPrec) / double (mPrecision - 1);
          for (unsigned i = 0; i < orderVec[iBloc][iPrec]; i++) {
            log.printf ("%6i%6i%12.6f%12.6f%12.6f%12.6f",
                        i0, i1, distance[0], distance[1], distance[2], cPIV);
            for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++){
              (counter[iRef])++;
              double rPIV = mRefPIV[iRef][iBloc][counter[iRef] - 1];
              log.printf ("%12.6f", rPIV);
            }
            log << "\n";
          }
        }
      }
    }
    log.printf ("This was a test, now exit \n");
    exit();
  }
  
  if (mTimer) stopwatch.start("4 Build For Derivatives");

  ///////////////////////////////////////////////////////////////////
  // compute derivatives with a given stride
  if (getStep() % mUpdateStride == 0 || mFirstStep) {
    // set to zero PIVdistance, derivatives and virial when they are calculated
    for (unsigned j = 0; j < mDerivatives.size(); j++) {
      for (unsigned k = 0; k < 3; k++) {
        mDerivatives[j][k] = 0.;
        if (j < 3) mVirial[j][k] = 0.;
      }
    }
    // compute PIV-PIV distance and derivatives for each reference
    unsigned i0 = 0;
    unsigned i1 = 0;
    double dm = 0;
    double tPIV = 0;
    for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++) {
      mDistancePIV[iRef] = 0.;
    }
    // Re-compute atomic distances for derivatives and compute PIV-PIV distance
    for (unsigned iBloc = 0; iBloc < mBlocks; iBloc++) {
      std::vector<unsigned> counter (mRefPIV.size(), 0);

      for (unsigned iPrec = 1; iPrec < mPrecision; iPrec++) {
        if (atmI0[iBloc][iPrec].size() > 0) {
          i0 = atmI0[iBloc][iPrec][0];
          i1 = atmI1[iBloc][iPrec][0];
          distance = distanceAB (getPosition (i0), getPosition (i1));
          // this is needed for dfunc and dervatives
          tPIV = 0.;
          dfunc = 0.;
          dm = distance.modulo();
          if (dm < mCutOff[iBloc]) {
            mSwitchFunc[iBloc].calculate (dm * mVolumeFactor, dfunc);
          }
          // PIV distance
          tPIV = double (iPrec) / double (mPrecision - 1);
          for (unsigned iRef = 0; i < mRefPIV.size(); iRef++) {
            for (unsigned i = 0; i < orderVec[iBloc][iPrec]; i++) {
              coord[iRef].push_back (tPIV - mRefPIV[iRef][iBloc][counter]);
              (counter[iRef])++;
            }
          }
          for (unsigned i = 0; i < atmI0[iBloc][iPrec].size(); i++) {
            for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++){
              double tmp = 2. * mBlockFactor[iBloc] * coord[iRef][i]
                           * mVolumeFactor*mVolumeFactor * dfunc;
              Vector tempDeriv = tmp * distance;
              // 0.5 * (x_i - x_k) * f_ik         (force on atom k due to atom i)
              unsigned index0 = atmI0[iBloc][iPrec][i];
              unsigned index1 = atmI1[iBloc][iPrec][i];
              mDerivatives[index0] -= tempDeriv;
              mDerivatives[index1] += tempDeriv;
              mVirial -= tmp * Tensor (distance, distance);         
              if (mScaleVolume) {
                mVirial += 1./3. * tmp * dm*dm * Tensor::identity();
              }
              mDistancePIV[iRef] += mBlockFactor[iBloc] * coord[iRef][i]*coord[iRef][i];
            }
          }
        } // if bin is not empty 
      } // loop on histogram
    } // loop on block

    if (!mSerial) {
      comm.Barrier();
      for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++) {
        comm.Sum(&mDistancePIV[iRef], 1);
      }
      if (!mDerivatives.empty()) {
        comm.Sum (&mDerivatives[0][0], 3 * mDerivatives.size());
      }
      comm.Sum (&mVirial[0][0], 9);
    }
  } // if update_piv
  
  //Timing
  if (mTimer) stopwatch.stop ("4 Build For Derivatives");
  if (mTimer) stopwatch.start ("5 Update Derivatives");

  // Update derivatives, virial
  setBoxDerivatives(mVirial);
  for (unsigned i = 0; i < mDerivatives.size(); ++i) {
    setAtomsDerivatives(i, mDerivatives[i]);
  }
  // update variables
  for (unsigned iRef = 0; iRef < mRefPIV.size(); iRef++) {
    Value* pValDistance = getPntrToComponent ("d" + std::to_string(iRef + 1));
    pValDistance->set (mDistancePIV[iRef]);
  }

  mFirstStep = false;
  
  if (mTimer) stopwatch.stop ("5 Update Derivatives");
  if (mTimer) {
    log.printf ("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log << stopwatch;
  }
} // end of calculate
} // close namespace piv
} // close namespace PLMD
