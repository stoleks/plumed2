/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 of Alexandre Jedrecy, Pipolo Silvio and Fabio Pietrucci.

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

using namespace std;

namespace PLMD
{
namespace piv
{

//+PLUMEDOC COLVAR PIVg
/*
Calculates the PIVg-distance: the squared Cartesian distance between the PIVg \cite gallet2013structural,pipolo2017navigating
associated to the configuration of the system during the dynamics and a reference configuration provided
as input (PDB file format).
PIVg can be used together with \ref FUNCPATHMSD to define a path in the PIVg space.
\par Examples

The following example calculates PIVg-distances from three reference configurations in Ref1.pdb, Ref2.pdb and Ref3.pdb
and prints the results in a file named colvar.
Three atoms with names (pdb file) A B and C are used to construct the PIVg and all PIVg blocks (AA, BB, CC, AB, AC, BC) are considered.
SFACTOR is a scaling factor that multiplies the contribution to the PIVg-distance given by the single PIVg block.
NLIST sets the use of neighbor lists for calculating atom-atom distances.
The SWITCH keyword specifies the perameters of the switching function that transforms atom-atom distances.
SORT=1 meand that the PIVg block elements are sorted (SORT=0 no sorting.)
Values for SORT, SFACTOR and Neighborlist parameters have to be specified for each block.
The order is the following: AA,BB,CC,AB,AC,BC. If ONLYDIRECT (ONLYCROSS) is used the order is AA,BB,CC (AB,AC,BC).
The sorting operation within each PIVg block is performed using the counting sort algorithm, you can use PRECISION to specify the size of the counting array.
\plumedfile
PIVg ...
LABEL=piv
NLIST
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
... PIVg

PRINT ARG=piv.d1,piv.d2,piv.d3 FILE=colvar
\endplumedfile

WARNING:
Both the "CRYST" and "ATOM" lines of the PDB files must conform precisely to the official pdb format, including the width of each alphanumerical field:

\verbatim
CRYST1   31.028   36.957   23.143  89.93  92.31  89.99 P 1           1
ATOM      1  OW1 wate    1      15.630  19.750   1.520  1.00  0.00
\endverbatim

In each pdb frame, atoms must be numbered in the same order and with the same element symbol as in the input of the MD program.

The following example calculates the PIVg-distances from two reference configurations Ref1.pdb and Ref2.pdb
and uses PIVg-distances to define a Path Collective Variable (\ref FUNCPATHMSD) with only two references (Ref1.pdb and Ref2.pdb).
With the VOLUME keyword one scales the atom-atom distances by the cubic root of the ratio between the specified value and the box volume of the initial step of the trajectory file.

\plumedfile
PIVg ...
LABEL=piv
VOLUME=12.15
NLIST
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
... PIVg

p1: FUNCPATHMSD ARG=piv.d1,piv.d2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=piv.d1,piv.d2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIVg please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class PIVg : public Colvar
{
public:
  PIVg (const ActionOptions&);
  ~PIVg ();
  static void registerKeywords (Keywords& keys);
  void checkFieldsAllowed () {}
  // active methods:
  virtual void calculate ();
private:
  Vector distanceAB (const Vector& A, const Vector& B);
private:
  ForwardDecl<Stopwatch> stopwatch_fwd;
  /// The stopwatch that times the different parts of the calculation
  Stopwatch& stopwatch = *stopwatch_fwd;
  // parameters
  bool mPBC, mTest, mCross, mDirect, mDoCOM, mTimer;
  bool mSerial, mDoNeigh, mScaleVolume, mComputeDerivative;
  unsigned mPrecision, mUdpatePIV, mAtoms, mBlocks, NLsize;
  // mVolumeFactor: volume scaling factor for distances
  double mVolumeFactor, mVolume0, mPIVdistance;
  Tensor mVirial;
  NeighborList* mAllAtoms;
  std::vector<bool> dosort;
  std::vector<double> scaling, mR_0;
  std::vector<double> mNeighborListSkin;
  std::vector<double> mMassFactor;
  std::vector<Vector> mDeriv;
  std::vector<Vector> mPosCOM;
  std::vector<string> mSwitchParam;
  std::vector<NeighborList*> mBlockAtoms;
  std::vector<NeighborList*> mBlockAtomCOM;
  std::vector<SwitchingFunction> mSwitchFunction;
  std::vector<std::vector<double>> mRefPIV;
};

///////////////////////////////////////////////////////////
// REGISTER ACTION AND KEYWORDS
///////////////////////////////////////////////////////////

PLUMED_REGISTER_ACTION(PIVg,"PIVg")

void PIVg::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add ("numbered", "SWITCH", "The switching functions parameter."
            "You should specify a Switching function for all "
            "PIV blocks. Details of the various switching functions "
            "you can use are provided on \\ref switchingfunction.");
  keys.add ("compulsory", "REF_FILE", "PDB file name that contains "
            "the i-th reference structure.");
  keys.add ("compulsory", "SORT", "Whether to sort or not the PIV block.");
  keys.add ("compulsory", "ATOMTYPES", "The atomtypes to use for PIV.");
  keys.add ("optional", "PRECISION", "The precision for reals ap"
            "proximation with integers in sorting.");
  keys.add ("optional", "SFACTOR", "Scale the PIV-distance by such"
            "block-specific factor");
  keys.add ("optional", "VOLUME", "Scale atom-atom distances by the"
            "cubic root of the cell volume. The input volume is used"
            "to scale the R_0 value of the switching function.");
  keys.add ("optional", "NL_CUTOFF", "Neighbour lists cutoff.");
  keys.add ("optional", "NL_STRIDE", "Update neighbour lists every "
            "NL_STRIDE steps.");
  keys.add ("optional", "NL_SKIN", "The maximum atom displacement "
            "tolerated for the neighbor lists update.");
  keys.add ("optional", "UPDATEPIV", "Frequency (timesteps) at which "
            "the PIV is updated.");
  keys.addFlag ("TEST", false, "Print current+reference PIV and exit");
  keys.addFlag ("COM", false, "Use centers of mass of groups of atoms "
               "instead of atoms as specified in the Pdb file");
  keys.addFlag ("ONLYCROSS", false, "Use only cross-terms (A-B, A-C, "
                "B-C, ...) in PIV");
  keys.addFlag ("ONLYDIRECT", false, "Use only direct-terms (A-A, B-B,"
                "C-C, ...) in PIV");
  keys.addFlag ("DERIVATIVES", false, "Activate the calculation of "
                "the PIV for every class (needed for numerical "
                "derivatives).");
  keys.addFlag ("NO_NLIST", false, "Use a neighbour list for "
                "distance calculations.");
  keys.addFlag ("SERIAL", false,"Perform the calculation in serial "
                "- for debug purpose");
  keys.addFlag ("TIMER",false,"Permorm timing analysis on heavy loops.");
  keys.reset_style("SWITCH", "compulsory");
}

///////////////////////////////////////////////////////////
// CONSTRUCTOR
///////////////////////////////////////////////////////////

PIVg::PIVg(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  mPBC (true),
  mTest (false),
  mCross (true),
  mDirect (true),
  mDoCOM (false),
  mTimer (false),
  mSerial (false),
  mDoNeigh (false),
  mScaleVolume (false),
  mComputeDerivative (false),
  mPrecision (1000),
  mUdpatePIV (1),
  mVolumeFactor (1.),
  mVolume0 (0.),
  mAtoms (1),
  /*
  dosort (std::vector<bool> (mBlocks)),
  scaling (std::vector<double> (mBlocks)),
  mR_0 (std::vector<double> (mBlocks)),
  mNeighborListSkin (std::vector<double>(mBlocks)),
  mMassFactor (std::vector<double>(mBlocks)),
  mRefPIV (std::vector<std::vector<double> >(mBlocks)),
  mSwitchParam (std::vector<string>(mBlocks)),
  mDeriv (std::vector<Vector>(1)),
  mPosCOM (std::vector<Vector>(NLsize))
  */
  mBlockAtoms (std::vector<NeighborList *> (mBlocks)),
  mBlockAtomCOM (std::vector<NeighborList *> (1))
{
  log << "Starting PIV Constructor\n";
  bool noPBC = !mPBC, noNeighbor = !mDoNeigh;
  bool oc = false, od = false;

  // parse all mandatory inputs that are not vector and flags
  parse ("VOLUME", mVolume0);
  parseFlag ("COM", mDoCOM);
  parseFlag ("TEST", mTest);
  parseFlag ("NOPBC", noPBC);
  parseFlag ("TIMER", mTimer);
  parseFlag ("SERIAL", mSerial);
  parseFlag ("NO_NLIST", noNeighbor);
  parseFlag ("ONLYCROSS", oc);
  parseFlag ("ONLYDIRECT", od);
  parseFlag ("DERIVATIVES", mComputeDerivative);

  // Precision on the real-to-integer transformation for the sorting
  if (keywords.exists ("PRECISION")) {
    parse ("PRECISION", mPrecision);
    if (mPrecision < 2) error ("Precision must be => 2");
  }
  // Stride for the PIV computation
  if (keywords.exists ("UPDATEPIV")) {
    parse ("UPDATEPIV",mUdpatePIV);
  }
  // scaling factor for the PIV block
  if (keywords.exists ("SFACTOR")) {
    parseVector ("SFACTOR",scaling);
  }

  // PBC
  mPBC = !noPBC;
  if (mPBC) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }

  // serial/parallel computation
  if (mSerial) {
    log << "Serial PIV construction\n";
  } else     {
    log << "Parallel PIV construction\n";
  }

  // Derivatives calculation for test purposes
  if (mComputeDerivative) log << "Computing Derivatives\n";

  // Timing
  if (mTimer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }

  // computing COM or not
  if (mDoCOM) {
    log << "Building PIV using COMs\n";
  }

  // Volume Scaling
  if (mVolume0 > 0) {
    mScaleVolume = true;
  }

  // PIV direct and cross blocks
  if (oc && od) {
    error ("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  if (oc) {
    mDirect = false;
    log << "Using only CROSS-PIV blocks\n";
  }
  if (od) {
    mCross = false;
    log << "Using only DIRECT-PIV blocks\n";
  }

  // Atoms for PIV
  std::vector<string> atomTypes;
  parseVector ("ATOMTYPES", atomTypes);
  mAtoms = atomTypes.size();

  // Reference PDB file
  std::string referenceFile;
  parse ("REF_FILE", referenceFile);
  PDB mypdb;
  FILE* fp = fopen (referenceFile.c_str(),"r");
  if (fp != NULL) {
    log << "Opening PDB file with reference frame: "
        <<referenceFile.c_str() << "\n";
    mypdb.readFromFilepointer (fp, 
            plumed.getAtoms ().usingNaturalUnits(),
            0.1/atoms.getUnits ().getLength());
    fclose (fp);
  } else {
    error ("Error in reference PDB file");
  }

  // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
  // Atomlist or Plist used to build pair lists
  std::vector<std::vector<AtomNumber>> Plist (mAtoms);
  // Atomlist used to build list of atoms for each COM
  std::vector<std::vector<AtomNumber>> comatm (1);
  // NLsize is the number of atoms in the pdb cell
  NLsize = mypdb.getAtomNumbers ().size();
  // In the following P stands for Point (either an Atom or a COM)
  unsigned resnum = 0;
  // Presind (array size: number of residues) contains the contains the residue number
  //   this is because the residue numbers may not alwyas be ordered from 1 to resnum
  std::vector<unsigned> Presind;
  // Build Presind
  for (unsigned i = 0; i < mypdb.getAtomNumbers().size(); i++) {
    unsigned rind=mypdb.getResidueNumber(mypdb.getAtomNumbers()[i]);
    bool oldres=false;
    for (unsigned j = 0; j < Presind.size(); j++) {
      if (rind == Presind[j]) {
        oldres = true;
      }
    }
    if (!oldres) {
      Presind.push_back (rind);
    }
  }
  resnum=Presind.size();

  // Pind0 is the atom/COM used in Nlists (for COM Pind0 is the first atom in the pdb belonging to that COM)
  unsigned Pind0size;
  if (mDoCOM) {
    Pind0size=resnum;
  } else {
    Pind0size=NLsize;
  }
  std::vector<unsigned> Pind0(Pind0size);
  // If COM resize important arrays
  comatm.resize(NLsize);
  if (mDoCOM) {
    mBlockAtomCOM.resize(NLsize);
    mPosCOM.resize(NLsize);
    mMassFactor.resize(NLsize,0.);
  }
  log << "Total COM/Atoms: " << mAtoms*resnum << " \n";
  // Build lists of Atoms/COMs for NLists
  //   comatm filled also for non_COM calculation for analysis purposes
  for (unsigned j=0; j < mAtoms; j++) {
    unsigned oind = 0;
    for (unsigned i=0; i < Pind0.size(); i++) {
      Pind0[i]=0;
    }
    for (unsigned i=0; i < mypdb.getAtomNumbers().size(); i++) {
      // Residue/Atom AtomNumber: used to build NL for COMS/Atoms pairs.
      AtomNumber anum=mypdb.getAtomNumbers()[i];
      // ResidueName/Atomname associated to atom
      string rname=mypdb.getResidueName(anum);
      string aname=mypdb.getAtomName(anum);
      // Index associated to residue/atom: used to separate COM-lists
      unsigned rind=mypdb.getResidueNumber(anum);
      unsigned aind=anum.index();
      // This builds lists for NL
      string Pname;
      unsigned Pind = 0;
      if (mDoCOM) {
        Pname = rname;
        for (unsigned l = 0; l < resnum; l++) {
          if (rind == Presind[l]) {
            Pind = l;
          }
        }
      } else {
        Pname=aname;
        Pind=aind;
      }
      if (Pname == atomTypes[j]) {
        if (Pind0[Pind] == 0) {
          // adding the atomnumber to the atom/COM list for pairs
          Plist[j].push_back(anum);
          Pind0[Pind]=aind+1;
          oind=Pind;
        }
        // adding the atomnumber to list of atoms for every COM/Atoms
        comatm[Pind0[Pind]-1].push_back(anum);
      }
    }
    // Output Lists
    log << "  Groups of type  " << j << ": " << Plist[j].size() << " \n";
    string gname;
    unsigned gsize;
    if (mDoCOM) {
      gname=mypdb.getResidueName(comatm[Pind0[oind]-1][0]);
      gsize=comatm[Pind0[oind]-1].size();
    } else {
      gname=mypdb.getAtomName(comatm[Pind0[oind]-1][0]);
      gsize=1;
    }
    log.printf("    %6s %3s %13s %10i %6s\n", "type  ", gname.c_str(),"   containing ",gsize," atoms");
  }

  // This is to build the list with all the atoms
  std::vector<AtomNumber> listall;
  for (unsigned i=0; i < mypdb.getAtomNumbers().size(); i++) {
    listall.push_back(mypdb.getAtomNumbers()[i]);
  }

  // PIV blocks and Neighbour Lists
  mBlocks=0;
  // Direct adds the A-A ad B-B blocks (N)
  if (mDirect) {
    mBlocks += unsigned(mAtoms);
  }
  // Cross adds the A-B blocks (N*(N-1)/2)
  if (mCross) {
    mBlocks += unsigned(double(mAtoms*(mAtoms-1))/2.);
  }
  // Resize vectors according to mBlocks
  mRefPIV.resize (mBlocks);

  // PIV scaled option
  scaling.resize (mBlocks);
  for (unsigned j=0; j < mBlocks; j++) {
    scaling[j] = 1.;
  }

  // Neighbour Lists option
  mBlockAtoms.resize (mBlocks);
  mNeighborListSkin.resize (mBlocks);
  mDoNeigh = !noNeighbor;
  if (mDoNeigh) {
    std::vector<double> neighborListCutOff (mBlocks, 0.);
    std::vector<int> neighborListStride (mBlocks, 0);
    parseVector ("NL_CUTOFF",neighborListCutOff);
    //if (neighborListCutOff.size()!=getNumberOfArguments() && neighborListCutOff.size()!=0) error ("not enough values for NL_CUTOFF");
    parseVector ("NL_STRIDE",neighborListStride);
    //if (neighborListStride.size()!=getNumberOfArguments() && neighborListStride.size()!=0) error ("not enough values for NL_STRIDE");
    parseVector ("NL_SKIN",mNeighborListSkin);
    //if (mNeighborListSkin.size()!=getNumberOfArguments() && mNeighborListSkin.size()!=0) error ("not enough values for NL_SKIN");
    for (unsigned j = 0; j < mBlocks; j++) {
      if (neighborListCutOff[j] <= 0.0) error ("NL_CUTOFF should be explicitly specified and positive");
      if (neighborListStride[j] <= 0) error ("NL_STRIDE should be explicitly specified and positive");
      if (mNeighborListSkin[j] <= 0.) error ("NL_SKIN should be explicitly specified and positive");
      neighborListCutOff[j]=neighborListCutOff[j]+mNeighborListSkin[j];
    }
    log << "Creating Neighbor Lists \n";
    // WARNING: is neighborListCutOff meaningful here?
    mAllAtoms= new NeighborList (
                     listall, mPBC, getPbc(),
                     neighborListCutOff[0],
                     neighborListStride[0]);
    if (mDoCOM) {
      //Build lists of Atoms for every COM
      for (unsigned i=0; i<mPosCOM.size(); i++) {
        // WARNING: is neighborListCutOff meaningful here?
        mBlockAtomCOM[i]= new NeighborList (
                                comatm[i], mPBC, getPbc(),
                                neighborListCutOff[0],
                                neighborListStride[0]);
      }
    }
    unsigned ncnt=0;
    // Direct blocks AA, BB, CC, ...
    if (mDirect) {
      for (unsigned j=0; j < mAtoms; j++) {
        mBlockAtoms[ncnt]= new NeighborList (
                        Plist[j], mPBC, getPbc(),
                        neighborListCutOff[j],
                        neighborListStride[j]);
        ncnt+=1;
      }
    }
    // Cross blocks AB, AC, BC, ...
    if (mCross) {
      for (unsigned j=0; j < mAtoms; j++) {
        for (unsigned i=j+1; i < mAtoms; i++) {
          mBlockAtoms[ncnt]= new NeighborList(
                          Plist[i], Plist[j], 
                          false, mPBC, getPbc(),
                          neighborListCutOff[ncnt],
                          neighborListStride[ncnt]);
          ncnt+=1;
        }
      }
    }
  } else {
    log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
    mAllAtoms= new NeighborList(listall, mPBC, getPbc());
    for (unsigned j=0; j < mBlocks; j++) {
      mBlockAtoms[j]= new NeighborList(
                   Plist[j], Plist[j],
                   true, mPBC, getPbc());
    }
  }
  // Output mBlocks
  log << "Total Nlists: " << mBlocks << " \n";
  for (unsigned j=0; j < mBlocks; j++) {
    log << "  list " << j+1 << "   size " << mBlockAtoms[j]->size() << " \n";
  }
  // Calculate COM masses once and for all from lists
  if (mDoCOM) {
    //log << "Computing COM masses  \n";
    for (unsigned j=0; j<mPosCOM.size(); j++) {
      double commass=0.;
      for (unsigned i=0; i<mBlockAtomCOM[j]->getFullAtomList().size(); i++) {
        unsigned andx=mBlockAtomCOM[j]->getFullAtomList()[i].index();
        commass+=mypdb.getOccupancy()[andx];
      }
      for (unsigned i=0; i<mBlockAtomCOM[j]->getFullAtomList().size(); i++) {
        unsigned andx=mBlockAtomCOM[j]->getFullAtomList()[i].index();
        if (commass>0.) {
          mMassFactor[andx]=mypdb.getOccupancy()[andx]/commass;
        } else {
          mMassFactor[andx]=1.;
        }
      }
    }
  }

  // Sorting
  dosort.resize(mBlocks);
  std::vector<int> ynsort(mBlocks);
  parseVector ("SORT",ynsort);
  for (unsigned i=0; i<mBlocks; i++) {
    if (ynsort[i] == 0 || mComputeDerivative) {
      dosort[i]=false;
    } else {
      dosort[i]=true;
    }
  }

  //build box vectors and correct for pbc
  log << "Building the box from PDB data ... \n";
  Tensor Box=mypdb.getBoxVec();
  log << "  Done! A,B,C vectors in Cartesian space:  \n";
  log.printf("  A:  %12.6f%12.6f%12.6f\n", Box[0][0],Box[0][1],Box[0][2]);
  log.printf("  B:  %12.6f%12.6f%12.6f\n", Box[1][0],Box[1][1],Box[1][2]);
  log.printf("  C:  %12.6f%12.6f%12.6f\n", Box[2][0],Box[2][1],Box[2][2]);
  log << "Changing the PBC according to the new box \n";
  Pbc mypbc;
  mypbc.setBox(Box);
  log << "The box volume is " << mypbc.getBox().determinant() << " \n";

  //Compute scaling factor
  if (mScaleVolume) {
    mVolumeFactor=cbrt(mVolume0/mypbc.getBox().determinant());
    log << "Scaling atom distances by  " << mVolumeFactor << " \n";
  } else {
    log << "Using unscaled atom distances \n";
  }

  mR_0.resize (mBlocks);
  mSwitchParam.resize (mBlocks);
  for (unsigned iBlock = 0; iBlock < mBlocks; iBlock++) {
    if (!parseNumbered("SWITCH", iBlock+1, mSwitchParam[iBlock])) {
      log << "Problem reading switching function parameters.\n";
      break;
    }
  }
  if (mComputeDerivative) {
    // Set switching function parameters here if computing derivatives
    log << "Switching Function Parameters\n";
    mSwitchFunction.resize (mBlocks);
    std::string errors;
    for (unsigned iBlock = 0; iBlock < mBlocks; iBlock++) {
      /*
      if (mScaleVolume) {
        double r0;
        vector<string> data = Tools::getWords(mSwitchParam[iBlock]);
        data.erase(data.begin());
        bool tmp=Tools::parse (data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=mVolumeFactor;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = mSwitchParam[iBlock].find("R_0");
        mSwitchParam[iBlock].replace(pos+4,old_r0.size(),new_r0);
      }
      */
      mSwitchFunction[iBlock].set (mSwitchParam[iBlock], errors);
      std::string num;
      Tools::convert(iBlock + 1, num);
      if (errors.length() != 0) {
        error ("problem reading SWITCH" +num+ " keyword : " + errors);
      }
      mR_0[iBlock] = mSwitchFunction[iBlock].get_r0();
      log << "  Swf: " << iBlock + 1 << "  r0=" 
          << (mSwitchFunction[iBlock].description()).c_str() << " \n";
    }
  }

  // build COMs from positions if requested
  if (mDoCOM) {
    for (unsigned j=0; j<mPosCOM.size(); j++) {
      mPosCOM[j][0]=0.;
      mPosCOM[j][1]=0.;
      mPosCOM[j][2]=0.;
      for (unsigned i=0; i<mBlockAtomCOM[j]->getFullAtomList().size(); i++) {
        unsigned andx=mBlockAtomCOM[j]->getFullAtomList()[i].index();
        mPosCOM[j]+=mMassFactor[andx]*mypdb.getPositions()[andx];
      }
    }
  }
  // build the mRefPIV distances (transformation and sorting is done afterwards)
  if (mComputeDerivative) {
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
  }
  for (unsigned j=0; j<mBlocks; j++) {
    for (unsigned i=0; i<mBlockAtoms[j]->size(); i++) {
      unsigned i0=(mBlockAtoms[j]->getClosePairAtomNumber(i).first).index();
      unsigned i1=(mBlockAtoms[j]->getClosePairAtomNumber(i).second).index();
      //calculate/get COM position of centers i0 and i1
      Vector Pos0,Pos1;
      if (mDoCOM) {
        //if (pbc) makeWhole();
        Pos0=mPosCOM[i0];
        Pos1=mPosCOM[i1];
      } else {
        Pos0=mypdb.getPositions()[i0];
        Pos1=mypdb.getPositions()[i1];
      }
      Vector ddist;
      if (mPBC) {
        ddist=mypbc.distance(Pos0,Pos1);
      } else {
        ddist=delta(Pos0,Pos1);
      }
      double df=0.;
      // Transformation and sorting done at the first timestep to solve the r0 definition issue
      if (mComputeDerivative) {
        mRefPIV[j].push_back(mSwitchFunction[j].calculate(ddist.modulo()*mVolumeFactor, df));
      } else {
        mRefPIV[j].push_back(ddist.modulo()*mVolumeFactor);
      }
    }
    if (mComputeDerivative) {
      if (dosort[j]) {
        std::sort(mRefPIV[j].begin(),mRefPIV[j].end());
      }
      int lmt0=0;
      int lmt1=0;
      for (unsigned i=0; i<mRefPIV[j].size(); i++) {
        if (int(mRefPIV[j][i]*double(mPrecision-1)) == 0) {
          lmt0+=1;
        }
        if (int(mRefPIV[j][i]*double(mPrecision-1)) == 1) {
          lmt1+=1;
        }
      }
      log.printf("       |%10i|%15i|%15i|%15i|\n", j, mRefPIV[j].size(), lmt0, lmt1);
    }
  }

  checkRead();
  // From the plumed manual on how to build-up a new Colvar
  addValueWithDerivatives();
  requestAtoms(mAllAtoms->getFullAtomList());
  setNotPeriodic();
  // getValue()->setPeridodicity(false);
  // set size of derivative vector
  mDeriv.resize(getNumberOfAtoms());
}

///////////////////////////////////////////////////////////
// DESTRUCTOR
///////////////////////////////////////////////////////////

PIVg::~PIVg()
{
  for (unsigned j=0; j<mBlocks; j++) {
    delete mBlockAtoms[j];
  }
  if (mDoCOM) {
    for (unsigned j=0; j<NLsize; j++) {
      delete mBlockAtomCOM[j];
    }
  }
  delete mAllAtoms;
}

///////////////////////////////////////////////////////////
// DISTANCE
///////////////////////////////////////////////////////////

inline Vector PIVg::distanceAB (const Vector& A, const Vector& B)
{
  if (mPBC) {
    return pbcDistance (A, B);
  } else {
    return delta (A, B);
  }
}

///////////////////////////////////////////////////////////
// CALCULATE
///////////////////////////////////////////////////////////

void PIVg::calculate()
{
  // Local variables
  // The following are probably needed as static arrays
  static int prev_stp=-1;
  static int init_stp=1;
  static std::vector<std::vector<Vector> > prev_pos(mBlocks);
  static std::vector<std::vector<double> > cPIV(mBlocks);
  static std::vector<std::vector<int> > Atom0(mBlocks);
  static std::vector<std::vector<int> > Atom1(mBlocks);
  std::vector<std::vector<int> > A0(mPrecision);
  std::vector<std::vector<int> > A1(mPrecision);
  unsigned stride=1;
  unsigned rank=0;

  if (!mSerial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  } else {
    stride=1;
    rank=0;
  }

  // Transform (and sort) the mRefPIV before starting the dynamics
  if (((prev_stp == -1) || (init_stp == 1)) &&!mComputeDerivative) {
    if (prev_stp!=-1){init_stp=0;}
    // Calculate the volume scaling factor
    if (mScaleVolume) {
      mVolumeFactor = cbrt(mVolume0/getBox().determinant());
    }
    //Set switching function parameters
    log << "\n";
    log << "REFERENCE PDB # " << prev_stp+2 << " \n";
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    mSwitchFunction.resize (mBlocks);
    std::string errors;
    for (unsigned iBlock = 0; iBlock < mBlocks; iBlock++) {
      if (mScaleVolume) {
        auto data = Tools::getWords (mSwitchParam[iBlock]);
        data.erase (data.begin());
        // reading old r0
        double r0;
        if (!Tools::parse (data, "R_0", r0)) {
          log << "Error with R_0 parameter for switching function\n";
        }
        std::string sR0;
        Tools::convert (r0, sR0);
        // computing new r0
        r0 *= mVolumeFactor;
        auto pos = mSwitchParam[iBlock].find ("R_0");
        mSwitchParam[iBlock].replace (pos + 4, sR0.size(), std::to_string (r0));
      }
      mSwitchFunction[iBlock].set (mSwitchParam[iBlock], errors);
      std::string num;
      Tools::convert (iBlock + 1, num);
      if (errors.length() != 0) {
        error ("problem reading SWITCH" +num+ " keyword : " + errors);
        }
      mR_0[iBlock] = mSwitchFunction[iBlock].get_r0();
      log << "  Swf: " << iBlock << "  r0=" 
          << (mSwitchFunction[iBlock].description()).c_str() << " \n";
    }
    //Transform and sort
    log << "Building Reference PIV Vector \n";
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    double df=0.;
    for (unsigned j = 0; j < mBlocks; j++) {
      for (unsigned i = 0; i < mRefPIV[j].size(); i++) {
        mRefPIV[j][i] = mSwitchFunction[j].calculate (mRefPIV[j][i], df);
      }
      if (dosort[j]) {
        std::sort (mRefPIV[j].begin(), mRefPIV[j].end());
      }
      int lmt0 = 0;
      int lmt1 = 0;
      for (unsigned i=0; i<mRefPIV[j].size(); i++) {
        if (int(mRefPIV[j][i]*double(mPrecision-1)) == 0) {
          lmt0 += 1;
        }
        if (int(mRefPIV[j][i]*double(mPrecision-1)) == 1) {
          lmt1 += 1;
        }
      }
      log.printf("       |%10i|%15i|%15i|%15i|\n", j, mRefPIV[j].size(), lmt0, lmt1);
    }
    log << "\n";
  }
  // Do the sorting only once per timestep to avoid building the PIV N times for N mRefPIV PDB structures!
  if ((getStep()>prev_stp && getStep() % mUdpatePIV == 0) || mComputeDerivative) {
    if (mComputeDerivative) log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV \n";
    //
    // build COMs from positions if requested
    if (mDoCOM) {
      if (mPBC) makeWhole();
      for (unsigned j=0; j<mPosCOM.size(); j++) {
        mPosCOM[j][0]=0.;
        mPosCOM[j][1]=0.;
        mPosCOM[j][2]=0.;
        for (unsigned i=0; i<mBlockAtomCOM[j]->getFullAtomList().size(); i++) {
          unsigned andx=mBlockAtomCOM[j]->getFullAtomList()[i].index();
          mPosCOM[j]+=mMassFactor[andx]*getPosition(andx);
        }
      }
    }
    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (mDoNeigh) {
      bool doupdate=false;
      // For the first step build previous positions = actual positions
      if (prev_stp == -1) {
        bool docom=mDoCOM;
        for (unsigned j=0; j<mBlocks; j++) {
          for (unsigned i=0; i<mBlockAtoms[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if (docom) {
              Pos=mPosCOM[i];
            } else {
              Pos=getPosition(mBlockAtoms[j]->getFullAtomList()[i].index());
            }
            prev_pos[j].push_back(Pos);
          }
        }
        doupdate=true;
      }
      // Decide whether to update lists based on atom displacement, every stride
      std::vector<std::vector<Vector> > tmp_pos(mBlocks);
      if (getStep() % mAllAtoms->getStride() == 0) {
        bool docom=mDoCOM;
        for (unsigned j=0; j<mBlocks; j++) {
          for (unsigned i=0; i<mBlockAtoms[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if (docom) {
              Pos=mPosCOM[i];
            } else {
              Pos=getPosition(mBlockAtoms[j]->getFullAtomList()[i].index());
            }
            tmp_pos[j].push_back(Pos);
            if (pbcDistance (tmp_pos[j][i], prev_pos[j][i]).modulo () 
                >=mNeighborListSkin[j]) {
              doupdate=true;
            }
          }
        }
      }
      // Update Nlists if needed
      if (doupdate == true) {
        for (unsigned j=0; j<mBlocks; j++) {
          for (unsigned i=0; i<mBlockAtoms[j]->getFullAtomList().size(); i++) {
            prev_pos[j][i]=tmp_pos[j][i];
          }
          mBlockAtoms[j]->update(prev_pos[j]);
          log << " Step " << getStep() << "  Neighbour lists updated " << mBlockAtoms[j]->size() << " \n";
        }
      }
    }
    // Calculate the volume scaling factor
    if (mScaleVolume) {
      mVolumeFactor = cbrt(mVolume0/getBox().determinant());
    }
    Vector ddist;
    // Global to local variables
    bool doserial=mSerial;
    // Build "mBlocks" PIV blocks
    for (unsigned j=0; j<mBlocks; j++) {
      if (dosort[j]) {
        // from global to local variables to speedup the for loop with if statements
        bool docom=mDoCOM;
        bool dopbc=mPBC;
        // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
        std::vector<int> OrdVec(mPrecision,0);
        cPIV[j].resize(0);
        Atom0[j].resize(0);
        Atom1[j].resize(0);
        // Building distances for the PIV vector at time t
        if (mTimer) stopwatch.start("1 Build cPIV");
        for (unsigned i=rank; i<mBlockAtoms[j]->size(); i+=stride) {
          unsigned i0=(mBlockAtoms[j]->getClosePairAtomNumber(i).first).index();
          unsigned i1=(mBlockAtoms[j]->getClosePairAtomNumber(i).second).index();
          Vector Pos0,Pos1;
          if (docom) {
            Pos0=mPosCOM[i0];
            Pos1=mPosCOM[i1];
          } else {
            Pos0=getPosition(i0);
            Pos1=getPosition(i1);
          }
          if (dopbc) {
            ddist=pbcDistance(Pos0,Pos1);
          } else {
            ddist=delta(Pos0,Pos1);
          }
          double df=0.;
          //Integer sorting ... faster!
          //Transforming distances with the Switching function + real to integer transformation
          int Vint=int(mSwitchFunction[j].calculate(ddist.modulo()*mVolumeFactor, df)*double(mPrecision-1)+0.5);
          //Integer transformed distance values as index of the Ordering Vector OrdVec
          OrdVec[Vint]+=1;
          //Keeps track of atom indices for force and virial calculations
          A0[Vint].push_back(i0);
          A1[Vint].push_back(i1);
        }
        if (mTimer) stopwatch.stop("1 Build cPIV");
        if (mTimer) stopwatch.start("2 Sort cPIV");
        if (!doserial) {
          // Vectors keeping track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          std::vector<int> Vdim(stride,0);
          std::vector<int> Vpos(stride,0);
          // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
          std::vector<int> OrdVecAll(stride*mPrecision);
          // Big vectors contining all Atom indexes for every occupancy (Atom0O(mPrecision,n) and Atom1O(mPrecision,n) matrices in one vector)
          std::vector<int> Atom0F;
          std::vector<int> Atom1F;
          // Vector used to reconstruct arrays
          std::vector<unsigned> k(stride,0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for (unsigned i=1; i<mPrecision; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=mPrecision-1
            // Can this be avoided ???
            Atom0F.insert(Atom0F.end(),A0[i].begin(),A0[i].end());
            Atom1F.insert(Atom1F.end(),A1[i].begin(),A1[i].end());
            A0[i].resize(0);
            A1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV block
          A0[0].resize(0);
          A1[0].resize(0);
          A0[mPrecision-1].resize(0);
          A1[mPrecision-1].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          OrdVec[0]=0;
          OrdVec[mPrecision-1]=0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim=Atom0F.size();
          // Vdim and Vpos keep track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim,1,&Vdim[0],1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int Fdim=0;
          for (unsigned i=1; i<stride; i++) {
            Vpos[i]=Vpos[i-1]+Vdim[i-1];
            Fdim+=Vdim[i];
          }
          Fdim+=Vdim[0];
          // build big vectors for atom pairs on all ranks for all ranks
          std::vector<int> Atom0FAll(Fdim);
          std::vector<int> Atom1FAll(Fdim);
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          //   Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV
          comm.Allgather(&OrdVec[0],mPrecision,&OrdVecAll[0],mPrecision);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&Atom0F[0],Atom0F.size(),&Atom0FAll[0],&Vdim[0],&Vpos[0]);
          comm.Allgatherv(&Atom1F[0],Atom1F.size(),&Atom1FAll[0],&Vdim[0],&Vpos[0]);

          // Reconstruct the full vectors from collections of Allgathered parts (this is a serial step)
          // This is the tricky serial step, to assemble toghether PIV and atom-pair info from head-tail big vectors
          // Loop before on l and then on i would be better but the allgather should be modified
          // Loop on blocks
          // for (unsigned m=0;m<mBlocks;m++) {
          // Loop on Ordering Vector size excluding zeros (i=1)
          if (mTimer) stopwatch.stop("2 Sort cPIV");
          if (mTimer) stopwatch.start("3 Reconstruct cPIV");
          for (unsigned i=1; i<mPrecision; i++) {
            // Loop on the ranks
            for (unsigned l=0; l<stride; l++) {
              // Loop on the number of head-to-tail pieces
              for (unsigned m=0; m<OrdVecAll[i+l*mPrecision]; m++) {
                // cPIV is the current PIV at time t
                cPIV[j].push_back(double(i)/double(mPrecision-1));
                Atom0[j].push_back(Atom0FAll[k[l]+Vpos[l]]);
                Atom1[j].push_back(Atom1FAll[k[l]+Vpos[l]]);
                k[l]+=1;
              }
            }
          }
          if (mTimer) stopwatch.stop("3 Reconstruct cPIV");
        } else {
          for (unsigned i=1; i<mPrecision; i++) {
            for (unsigned m=0; m<OrdVec[i]; m++) {
              cPIV[j].push_back(double(i)/double(mPrecision-1));
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
  if (mScaleVolume) {
    mVolumeFactor = cbrt(mVolume0/getBox().determinant());
  }

  // This test may be run by specifying the TEST keyword as input, it pritnts mRefPIV and cPIV and quits
  if (mTest) {
    unsigned limit=0;
    for (unsigned j=0; j<mBlocks; j++) {
      if (dosort[j]) {
        limit = cPIV[j].size();
      } else {
        limit = mRefPIV[j].size();
      }
      log.printf("PIV Block:  %6i %12s %6i \n", j, "      Size:", limit);
      log.printf("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV   ","    r-PIV   ","   i-j distance vector       ");
      for (unsigned i=0; i<limit; i++) {
        unsigned i0=0;
        unsigned i1=0;
        if (dosort[j]) {
          i0=Atom0[j][i];
          i1=Atom1[j][i];
        } else {
          i0=(mBlockAtoms[j]->getClosePairAtomNumber(i).first).index();
          i1=(mBlockAtoms[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector Pos0,Pos1;
        if (mDoCOM) {
          Pos0=mPosCOM[i0];
          Pos1=mPosCOM[i1];
        } else {
          Pos0=getPosition(i0);
          Pos1=getPosition(i1);
        }
        if (mPBC) {
          distance=pbcDistance(Pos0,Pos1);
        } else {
          distance=delta(Pos0,Pos1);
        }
        dfunc=0.;
        double cP,rP;
        if (dosort[j]) {
          cP = cPIV[j][i];
          rP = mRefPIV[j][mRefPIV[j].size()-cPIV[j].size()+i];
        } else {
          double dm=distance.modulo();
          cP = mSwitchFunction[j].calculate(dm*mVolumeFactor, dfunc);
          rP = mRefPIV[j][i];
        }
        log.printf("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
      }
    }
    log.printf("This was a test, now exit \n");
    exit();
  }

  if (mTimer) stopwatch.start("4 Build For Derivatives");
  // non-global variables Nder and Scalevol defined to speedup if structures in cycles
  bool Nder=mComputeDerivative;
  bool Scalevol = mScaleVolume;
  if (getStep()%mUdpatePIV == 0) {
    // set to zero PIVdistance, derivatives and virial when they are calculated
    for (unsigned j=0; j<mDeriv.size(); j++) {
      for (unsigned k=0; k<3; k++) {mDeriv[j][k]=0.;}
    }
    for (unsigned j=0; j<3; j++) {
      for (unsigned k=0; k<3; k++) {
        mVirial[j][k]=0.;
      }
    }
    mPIVdistance=0.;
    // Re-compute atomic distances for derivatives and compute PIV-PIV distance
    for (unsigned j=0; j<mBlocks; j++) {
      unsigned limit=0;
      // dosorting definition is to speedup if structure in cycles with non-global variables
      bool dosorting=dosort[j];
      bool docom=mDoCOM;
      bool dopbc=mPBC;
      if (dosorting) {
        limit = cPIV[j].size();
      } else {
        limit = mRefPIV[j].size();
      }
      for (unsigned i = rank; i < limit; i += stride) {
        unsigned i0 = 0;
        unsigned i1 = 0;
        if (dosorting) {
          i0 = Atom0[j][i];
          i1 = Atom1[j][i];
        } else {
          i0 = (mBlockAtoms[j]->getClosePairAtomNumber(i).first).index();
          i1=(mBlockAtoms[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector Pos0,Pos1;
        if (docom) {
          Pos0=mPosCOM[i0];
          Pos1=mPosCOM[i1];
        } else {
          Pos0=getPosition(i0);
          Pos1=getPosition(i1);
        }
        if (dopbc) {
          distance=pbcDistance(Pos0,Pos1);
        } else {
          distance=delta(Pos0,Pos1);
        }
        dfunc=0.;
        // this is needed for dfunc and dervatives
        double dm=distance.modulo();
        double tPIV = mSwitchFunction[j].calculate(dm*mVolumeFactor, dfunc);
        // PIV distance
        double coord=0.;
        if (!dosorting || Nder) {
          coord = tPIV - mRefPIV[j][i];
        } else {
          coord = cPIV[j][i] - mRefPIV[j][mRefPIV[j].size()-cPIV[j].size()+i];
        }
        // Calculate derivatives, virial, and variable=sum_j (scaling[j] *(cPIV-mRefPIV)_j^2)
        // WARNING: dfunc=dswf/(mVolumeFactor*dm)  (this may change in future Plumed versions)
        double tmp = 2.*scaling[j]*coord*mVolumeFactor*mVolumeFactor*dfunc;
        Vector tmpder = tmp*distance;
        // 0.5*(x_i-x_k)*f_ik         (force on atom k due to atom i)
        if (docom) {
          Vector dist;
          for (unsigned k=0; k<mBlockAtomCOM[i0]->getFullAtomList().size(); k++) {
            unsigned x0=mBlockAtomCOM[i0]->getFullAtomList()[k].index();
            mDeriv[x0] -= tmpder*mMassFactor[x0];
            for (unsigned l=0; l<3; l++) {
              dist[l]=0.;
            }
            Vector P0=getPosition(x0);
            for (unsigned l=0; l<mBlockAtomCOM[i0]->getFullAtomList().size(); l++) {
              unsigned x1=mBlockAtomCOM[i0]->getFullAtomList()[l].index();
              Vector P1=getPosition(x1);
              if (dopbc) {
                dist+=pbcDistance(P0,P1);
              } else {
                dist+=delta(P0,P1);
              }
            }
            for (unsigned l=0; l<mBlockAtomCOM[i1]->getFullAtomList().size(); l++) {
              unsigned x1=mBlockAtomCOM[i1]->getFullAtomList()[l].index();
              Vector P1=getPosition(x1);
              if (dopbc) {
                dist+=pbcDistance(P0,P1);
              } else {
                dist+=delta(P0,P1);
              }
            }
            mVirial    -= 0.25*mMassFactor[x0]*Tensor(dist,tmpder);
          }
          for (unsigned k=0; k<mBlockAtomCOM[i1]->getFullAtomList().size(); k++) {
            unsigned x1=mBlockAtomCOM[i1]->getFullAtomList()[k].index();
            mDeriv[x1] += tmpder*mMassFactor[x1];
            for (unsigned l=0; l<3; l++) {
              dist[l]=0.;
            }
            Vector P1=getPosition(x1);
            for (unsigned l=0; l<mBlockAtomCOM[i1]->getFullAtomList().size(); l++) {
              unsigned x0=mBlockAtomCOM[i1]->getFullAtomList()[l].index();
              Vector P0=getPosition(x0);
              if (dopbc) {
                dist+=pbcDistance(P1,P0);
              } else {
                dist+=delta(P1,P0);
              }
            }
            for (unsigned l=0; l<mBlockAtomCOM[i0]->getFullAtomList().size(); l++) {
              unsigned x0=mBlockAtomCOM[i0]->getFullAtomList()[l].index();
              Vector P0=getPosition(x0);
              if (dopbc) {
                dist+=pbcDistance(P1,P0);
              } else {
                dist+=delta(P1,P0);
              }
            }
            mVirial    += 0.25*mMassFactor[x1]*Tensor(dist,tmpder);
          }
        } else {
          mDeriv[i0] -= tmpder;
          mDeriv[i1] += tmpder;
          mVirial    -= tmp*Tensor(distance,distance);
        }
        if (Scalevol) {
          mVirial+=1./3.*tmp*dm*dm*Tensor::identity();
        }
        mPIVdistance    += scaling[j]*coord*coord;
      }
    }

    if (!mSerial) {
      comm.Barrier();
      comm.Sum(&mPIVdistance,1);
      if (!mDeriv.empty()) comm.Sum(&mDeriv[0][0],3*mDeriv.size());
      comm.Sum(&mVirial[0][0],9);
    }
  }
  prev_stp=getStep();

  //Timing
  if (mTimer) stopwatch.stop("4 Build For Derivatives");
  if (mTimer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }

  // Update derivatives, virial, and variable (PIV-distance^2)
  for (unsigned i = 0; i < mDeriv.size(); ++i) setAtomsDerivatives(i,mDeriv[i]);
  setValue (mPIVdistance);
  setBoxDerivatives (mVirial);
}
//Close Namespaces at the very beginning
}
}

