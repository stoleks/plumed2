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

#ifndef PIV_WIP_H
#define PIV_WIP_H

#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace PLMD
{
namespace piv
{

//+PLUMEDOC COLVAR PIV
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
PIV ...
LABEL=piv
REF_FILE1=Ref1.pdb
REF_FILE2=Ref2.pdb
REF_FILE3=Ref2.pdb
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8
NL_STRIDE=10
NL_SKIN=0.1
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
PIV ...
LABEL=piv
VOLUME=12.15
REF_FILE1=Ref1.pdb
REF_FILE2=Ref2.pdb
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.5 MM=10 NN=5}
NL_CUTOFF=1.2
NL_STRIDE=10
NL_SKIN=0.1
... PIV

p1: FUNCPATHMSD ARG=piv.d1,piv.d2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=piv.d1,piv.d2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIV please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

using PtrNeighborList = std::unique_ptr <NeighborList>;

class PIV : public Colvar
{
public:
  /**
   * parse parameters and set-up PIV calculation
   */
  PIV (const ActionOptions&);
  /**
   * compute PIV each step
   */
  virtual void calculate ();
  /**
   * plumed stuff
   */
  void checkFieldsAllowed () {}
  static void registerKeywords ( Keywords& keys );
private:
  /**
   * compute delta or pbc distance if set
   */
  Vector distanceAB (
         const Vector& A,
         const Vector& B);
  /**
   * called at first step to initialize switching
   * functions, references and neighbors list
   */
  void initializeSwitchingFunction (bool doScaleVolume = true);
  void initializeReferencePIV ();
  void initializeNeighborList ();
  /**
   * called each step to update neighbor list
   */
  void updateNeighborList ();
  /**
   * parse all options and parameters
   */
  void parseOptions (
         bool& cross,
         bool& direct,
         double& neighborCut, 
         unsigned& neighborStride,
         std::vector <std::string>& atomTypes);
  /**
   * parse reference pdb
   */
  void parseReferences (
         const bool cross,
         const bool direct,
         const double neighborCut, 
         const unsigned neighborStride,
         const std::vector <std::string>& atomTypes);
  /**
   * initialize neighbor list with pdb data
   */
  void makeNeighborLists (
         const bool cross,
         const bool direct,
         const double neighborCut, 
         const unsigned neighborStride,
         const std::vector <AtomNumber>& allAtomsList,
         const std::vector <std::vector <AtomNumber>>& pairList,
         const std::vector <std::vector <AtomNumber>>& comAtoms);
private:
  ForwardDecl<Stopwatch> stopwatch_fwd;
  /// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch = *stopwatch_fwd;
  bool mPBC;
  bool mSerial;
  bool mTimer;
  bool mFirstStep;
  bool mScaleVolume;
  bool mTest;
  bool mComputeDerivatives;
  bool mDoCom;
  unsigned mPrecision;
  unsigned mAtomTypes;
  unsigned mBlocks;
  unsigned mUpdateStride;
  double mVolumeFactor;
  double mVolume0;
  double mLambda;
  double mAtomsSkin;
  Tensor mVirial;
  PtrNeighborList mBlockAtomsAll;
  std::vector<bool> mDoSort;
  std::vector<double> mBlockFactor;
  std::vector<double> mMassFactor;
  std::vector<double> mDistancePIV;
  std::vector<Vector> mPosCOM;
  std::vector<Vector> mDerivatives;
  std::vector<std::string> mSwitchData;
  std::vector<PtrNeighborList> mBlockAtoms;
  std::vector<PtrNeighborList> mBlockAtomCOM;
  std::vector<SwitchingFunction> mSwitchFunc;
  // 2-dimensional vector (block, atoms)
  std::vector<std::vector<Vector>> mPrevPosition;
  // 3-dimensional vector (reference, block, atoms)
  std::vector<std::vector<std::vector<double>>> mRefPIV;
};

PLUMED_REGISTER_ACTION (PIV, "PIV")

void PIV::registerKeywords (Keywords& keys)
{
  Colvar::registerKeywords (keys);
  // input 
  keys.add (
    "numbered", "SWITCH",
    "Switching functions parameter. You should specify "
    "a Switching function for all PIV blocks. Details of "
    "the various functions are provided on \\ref switchingfunction."
  );
  keys.add (
    "numbered", "REF_FILE",
    "PDB file that contains the i-th reference structure. "
    "If you indicate n reference, you will get n PIV distances."
  );
  keys.add (
    "compulsory", "ATOMTYPES",
    "The atomtypes to use for PIV."
  );
  keys.add (
    "optional", "SORT",
    "Whether to sort or not PIV blocks."
  );
  keys.add (
    "optional", "PRECISION",
    "Precision for approximating reals with integers in sorting."
  );
  keys.add (
    "optional", "SFACTOR",
    "Scale the PIV-distance "
    "by such block-specific factor"
  );
  keys.add (
    "optional", "VOLUME",
    "Scale atom-atom distances by the cubic root of the cell volume. "
    "Input volume is used to scale the R_0 value of the switching function."
  );
  keys.add (
    "optional", "UPDATEPIV",
    "Update PIV every UPDATEPIV steps."
  );
  keys.add (
    "optional", "NL_CUTOFF",
    "Neighbour lists cutoff."
  );
  keys.add (
    "optional", "NL_STRIDE",
    "Update neighbour lists every NL_STRIDE steps."
  );
  keys.add (
    "optional", "NL_SKIN", 
    "Maximum atom displacement accepted for the neighbor lists update."
  );
  keys.addFlag (
    "TEST", false,
    "Print actual and reference PIV, exit"
  );
  keys.addFlag (
    "COM", false,
    "Use centers of mass of groups of atoms instead of atoms as specified in the .pdb file"
  );
  keys.addFlag (  
    "ONLYCROSS", false,
    "Use only cross-terms in adjancy matrix (A-B, A-C, B-C ...)"
  );
  keys.addFlag (
    "ONLYDIRECT", false,
    "Use only direct-terms in adjancy matrix (A-A, B-B, ...)"
  );
  keys.addFlag (
    "DERIVATIVES", false,
    "Activate computation of the PIV for every class"
    " (needed for numerical derivatives)."
  );
  keys.addFlag (
    "SERIAL", false,
    "Compute in serial mode (use it if plumed was not built with MPI)"
  );
  keys.addFlag (
    "TIMER", false,
    "Perform timing analysis on heavy loops to profile code - it's a debug option."
  );
  keys.reset_style ("SWITCH", "compulsory");
  componentsAreNotOptional (keys);
  // output we add d1,...,d5 to avoid warnings in typical use case
  keys.addOutputComponent (
    "lambda", "default",
    "Lambda for pathCV such that: lambda * D_12 = 2.3"
  );
  keys.addOutputComponent (
    "d1", "default",
    "PIV distance between 1st reference and current state."
  );
  keys.addOutputComponent (
    "d2", "default",
    "PIV distance between 2nd reference and current state."
  );
  keys.addOutputComponent (
    "d3", "default",
    "PIV distance between 3rd reference and current state."
  );
  keys.addOutputComponent (
    "d4", "default",
    "PIV distance between 4th reference and current state."
  );
  keys.addOutputComponent (
    "d5", "default",
    "PIV distance between 5th reference and current state."
  );
}

} // close namespace piv
} // close namespace PLMD

#endif // PIV_WIP_H
