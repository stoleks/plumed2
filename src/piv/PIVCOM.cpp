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

#include "PIVCOM.h"

namespace std {
  // note: this implementation does not disable 'this' overload for array types
  template <typename Type, typename... Args>
  std::unique_ptr<Type> make_unique (Args&&... args)
  {
    return std::unique_ptr<Type> (new Type (std::forward<Args>(args)...));
  }
}

namespace PLMD
{
namespace piv
{
//---------------------------------------------------------
// CONSTRUCTOR
//---------------------------------------------------------

PIVCOM::PIVCOM (const ActionOptions&ao):
  PLUMED_COLVAR_INIT (ao),
  mPBC (true),
  mSerial (false),
  mTimer (false),
  mFirstStep (true),
  mScaleVolume (false),
  mTest (false),
  mComputeDerivatives (false),
  mPrecision (1000),
  mUpdateStride (1),
  mVolumeFactor (1.),
  mVolume0 (1.),
  mLambda (1.)
{
  log << "Starting PIVCOM Constructor\n";
  // parse options and parameters
  auto cross = true;
  auto direct = true;
  auto neighborCut = double (0);
  auto neighborStride = unsigned (0);
  auto atomTypes = std::vector <std::string> ();
  parseOptions (
    cross,
    direct,
    neighborCut, 
    neighborStride,
    atomTypes
  );
  // parse reference files
  parseReferences (
    cross,
    direct,
    neighborCut,
    neighborStride,
    atomTypes
  );
  // Plumed set-up
  checkRead();
  // add N distance for the N reference states
  for (unsigned ref = 1; ref <= mRefPIV.size(); ref++) {
    addComponentWithDerivatives ("d" + std::to_string (ref));
    componentIsNotPeriodic ("d" + std::to_string (ref));
  }
  // add the lambda component for the path collective variables
  addComponent ("lambda");
  componentIsNotPeriodic ("lambda");
  requestAtoms(mBlockAtomsAll->getFullAtomList());
  // set size of derivative vector
  mDerivatives.resize (getNumberOfAtoms()); 
}

//---------------------------------------------------------
// PARSE-OPTIONS
// parse all mandatory and optional parameters
//---------------------------------------------------------
void PIVCOM::parseOptions (
  bool& cross,
  bool& direct,
  double& neighborCut,
  unsigned& neighborStride,
  std::vector <std::string>& atomTypes) 
{
  auto onlyCross = false;
  auto onlyDirect = false;
  auto noPBC = !mPBC;
  parse ("NL_CUTOFF", neighborCut);
  parse ("NL_STRIDE", neighborStride);
  parse ("NL_SKIN", mAtomsSkin);
  parseFlag ("VOLUME", mScaleVolume);
  parseFlag ("TEST", mTest);
  parseFlag ("NOPBC", noPBC);
  parseFlag ("TIMER", mTimer);
  parseFlag ("SERIAL", mSerial);
  parseFlag ("ONLYCROSS", onlyCross);
  parseFlag ("ONLYDIRECT", onlyDirect);
  parseFlag ("DERIVATIVES", mComputeDerivatives);

  // parse the atom names 
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
  // PBC option
  mPBC = !noPBC;
  if (mPBC) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }
  // Serial or parallel option
  if (mSerial) {
    log << "Serial PIV construction\n";
  } else     {
    log << "Parallel PIV construction\n";
  }
  // Derivatives option
  if (mComputeDerivatives) {
    log << "Computing Derivatives\n";
  }
  // Timing option
  if (mTimer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }
  // Center of Mass option
  log << "Building PIV using COMs\n";

  // PIV direct and cross blocks option
  if (onlyCross && onlyDirect) {
    error ("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  cross = !onlyDirect;
  direct = !onlyCross;
  mBlocks = 0;
  if (cross) {
    mBlocks += static_cast<unsigned> (
      static_cast<double> (mAtomTypes * (mAtomTypes - 1)) / 2.
    );
    log << "Using CROSS-PIV blocks. ";
  }
  if (direct) {
    mBlocks += static_cast<unsigned> (mAtomTypes);
    log << "Using DIRECT-PIV blocks. ";
  }
  log << "There are " << mBlocks << " PIV blocks.\n";

  // resize all vector according to mBlocks
  mSwitchData.resize (mBlocks); 
  mPrevPosition.resize (mBlocks);
    
  // set neighborlist parameters
  if (neighborCut <= 0.0) {
    error("NL_CUTOFF should be explicitly specified and positive");
  }
  if (neighborStride <= 0) {
    error("NL_STRIDE should be explicitly specified and positive");
  }
  if (mAtomsSkin <= 0.) {
    error("NL_SKIN should be explicitly specified and positive");
  }
  neighborCut = neighborCut + mAtomsSkin;
  log << "Neighbor list cut-off=" << neighborCut
      << ", stride=" << neighborStride
      << ", shell=" << mAtomsSkin << "\n";

  // By default we sort all blocks
  mDoSort.resize (mBlocks, true);
  if (keywords.exists ("SORT")) {
    auto doSort = std::vector <unsigned> (mBlocks, 1);
    parseVector ("SORT", doSort);
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      mDoSort[bloc] = (!(doSort[bloc] == 0 || mComputeDerivatives));
    }
  }
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    mDoSort[bloc] ? log << "Sort " : log << "Don't sort ";
    log << "block " << bloc + 1 << "\n";
  }

  // PIV scaled option
  mScalingBlockFactor = std::vector<double> (mBlocks, 1.);
  if (keywords.exists ("SFACTOR")) {
    parseVector ("SFACTOR", mScalingBlockFactor);
  }

  // read parameters of switching functions 
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    if ( !parseNumbered ("SWITCH", bloc + 1, mSwitchData[bloc]) ){
      log << "Problem while reading the switching function parameters.\n";
      break;
    }
  }
  // et-up it here only if computing derivatives option is set
  if (mComputeDerivatives) {
    initializeSwitchingFunction (); // (false);
  }
}

//---------------------------------------------------------
// PARSE-REFERENCES
// parse all pdb files
//---------------------------------------------------------
void PIVCOM::parseReferences (
  const bool cross,
  const bool direct,
  const double neighborCut, 
  const unsigned neighborStride,
  const std::vector <std::string>& atomTypes)
{
  // read all reference files and compute mVolume0
  bool readAllRef = false;
  auto referenceCount = unsigned (0);
  auto referencePBC = std::vector <Pbc> ();
  auto referencePDB = std::vector <PDB> ();
  mVolume0 = 0;
  while (!readAllRef) {
    // set-up opening of the reference file
    std::string referenceFile; 
    parseNumbered ("REF_FILE", referenceCount + 1, referenceFile);

    // if file is empty we stop the procedure
    readAllRef = referenceFile.empty();
    if (readAllRef) break;

    // opening of the reference PBD file
    log << "\n----- Reference " << referenceCount + 1 << " -----\n";
    auto myPDB = PDB ();
    auto file = fopen (referenceFile.c_str(), "r");
    if (file != NULL) {
      log << "Opening PDB file with reference frame: " 
          << referenceFile << "\n";
      myPDB.readFromFilepointer (
        file, 
        plumed.getAtoms ().usingNaturalUnits (), 
        0.1 / atoms.getUnits ().getLength ()
      );
      fclose (file);
    } else {
      error ("Error in reference PDB file " + referenceFile);
    }
 
    // number of atoms in the PDB cell
    auto pdbAtoms = myPDB.getAtomNumbers();
    log << "Atoms numbers " << pdbAtoms.size () << "\n";

    // Build list of residue index
    auto residueIndex = std::unordered_map <unsigned, unsigned> ();
    for (unsigned i = 0; i < pdbAtoms.size (); i++) {
      auto resNumber = myPDB.getResidueNumber (pdbAtoms[i]);
      if (residueIndex.find (resNumber) == std::end (residueIndex)) {
        residueIndex[resNumber] = resNumber;
      }
    }

    // set-up center of mass parameter
    auto comAtoms = std::vector <std::vector <AtomNumber>> (pdbAtoms.size ());
    mPosCOM.resize (pdbAtoms.size ());
    mMassFactor.resize (pdbAtoms.size (), 0.);
    log << "Total COM/Atoms: " << mAtomTypes * residueIndex.size () << "\n";
    
    // build list of pair and index
    auto indexList = std::vector<unsigned> ();
    auto pairList = std::vector <std::vector <AtomNumber>> (mAtomTypes);
    for (unsigned type = 0; type < mAtomTypes; type++) {
      for (const auto& atom : pdbAtoms) {
        auto atomName = myPDB.getAtomName (atom);
        if (atomName == atomTypes [type]) {
          pairList [type].push_back (atom);
          indexList.push_back (atom.index () + 1);
        }
        comAtoms [indexList.back () - 1].push_back (atom);
      }
      // Output Lists
      log << "  Groups of type  " << type << ": " << pairList [type].size() << " \n";
      auto groupName = myPDB.getResidueName (comAtoms [indexList.back () - 1][0]);
      auto groupSize = comAtoms [indexList.back () - 1].size ();
      log.printf ("    %6s %3s %13s %10i %6s\n", "type  ", groupName.c_str(),
                  "   containing ", groupSize," atoms");
    }
    // create neighbor lists
    makeNeighborLists (
      cross,
      direct,
      neighborCut,
      neighborStride,
      pdbAtoms,
      pairList,
      comAtoms
    );
    // Calculate COM masses once and for all from lists
    for (unsigned i = 0; i < mPosCOM.size(); i++) {
      auto massCOM = double (0.);
      for (const auto& atom : mBlockAtomCOM[i]->getFullAtomList()) {
        massCOM += myPDB.getOccupancy()[atom.index()];
      }
      for (const auto& atom : mBlockAtomCOM[i]->getFullAtomList()) {
        if (massCOM > 0.) {
          mMassFactor[atom.index()] = 
            myPDB.getOccupancy()[atom.index()] / massCOM;
        } else {
          mMassFactor[atom.index()] = 1.;
        }
      }
    }
    // build box vectors
    log << "Building the box from PDB data ... \n";
    auto box = myPDB.getBoxVec();
    log << "  Done! A,B,C vectors in Cartesian space:  \n";
    log.printf ("  A:  %12.6f%12.6f%12.6f\n", box[0][0], box[0][1], box[0][2]);
    log.printf ("  B:  %12.6f%12.6f%12.6f\n", box[1][0], box[1][1], box[1][2]);
    log.printf ("  C:  %12.6f%12.6f%12.6f\n", box[2][0], box[2][1], box[2][2]);
    log << "Changing the PBC according to the new box \n";
    // set PBC box and compute reference volume
    auto myPBC = Pbc();
    myPBC.setBox (box);
    mVolume0 += myPBC.getBox().determinant();
    log << "The box volume is " << myPBC.getBox().determinant() << " \n";
    // count number of references and save pdb and pbc
    referencePBC.push_back (myPBC);
    referencePDB.push_back (myPDB);
    referenceCount++;
  } // while reading reference files

  // compute volume factor as average of reference's volume
  mVolume0 /= static_cast <double> (referenceCount);
  log << "\nAverage scaling volume V_0 = " << mVolume0 << "\n";

  // compute reference PIV
  for (unsigned ref = 0; ref < referenceCount; ref++) {
    //Compute scaling factor
    log << "For reference " << ref;
    if (mScaleVolume) {
      mVolumeFactor = cbrt (mVolume0 / referencePBC [ref].getBox().determinant());
      log << ", scaling atom distances by  " << mVolumeFactor << " \n";
    } else {
      log << ", using unscaled atom distances \n";
    }

    // build COMs from positions if requested
    for (unsigned i = 0; i < mPosCOM.size(); i++) {
      mPosCOM[i][0] = 0.;
      mPosCOM[i][1] = 0.;
      mPosCOM[i][2] = 0.; 
      for (const auto& atom : mBlockAtomCOM[i]->getFullAtomList()) {
        auto index = atom.index ();
        mPosCOM[i] += mMassFactor[index] * referencePDB [ref].getPositions()[index];
      }
    }
    
    // build the rPIV distances (transformation and sorting is done afterwards)
    auto refPIV = std::vector <std::vector <double>> ();
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      auto blockRefPIV = std::vector <double> ();
      for (unsigned atom = 0; atom < mBlockAtoms[bloc]->size(); atom++) {
        auto i0 = (mBlockAtoms[bloc]->getClosePairAtomNumber (atom).first).index();
        auto i1 = (mBlockAtoms[bloc]->getClosePairAtomNumber (atom).second).index();
        // get COM position of centers i0 and i1
        auto position0 = mPosCOM[i0];
        auto  position1 = mPosCOM[i1];
        Vector pairDist;
        if (mPBC) {
          pairDist = referencePBC [ref].distance (position0, position1);
        } else {
          pairDist = delta (position0, position1);
        }
        // transform + sort are done at first timestep to solve the r0 def issue
        blockRefPIV.push_back (pairDist.modulo () * mVolumeFactor);
      }
      refPIV.push_back (blockRefPIV);
    }
    mRefPIV.push_back (refPIV);
  } // loop over the number of references
  
  // print reference
  for (unsigned ref = 0; ref < referenceCount; ref++) {
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      log << "Reference " << ref + 1 << ", PIV block " << bloc + 1 << " has size: " 
          << mRefPIV [ref] [bloc].size() << " over a total of "
           << mBlockAtoms [bloc]->size() << " atoms-atoms pair\n";
      if (mComputeDerivatives) {
        if (mDoSort[bloc]) {
          std::sort(
            std::begin (mRefPIV [ref][bloc]),
            std::end (mRefPIV [ref][bloc])
          );
        }
        unsigned lmt0 = 0;
        unsigned lmt1 = 0;
        for (unsigned atm = 0; atm < mRefPIV [ref][bloc].size(); atm++) {
          if (mRefPIV [ref][bloc][atm] > 0.9) { lmt1++; }
          if (mRefPIV [ref][bloc][atm] < 0.1) { lmt0++; }
        }
        if (mComputeDerivatives && bloc == 0) {
          log << "  PIV  |  block   |     Size      |     "
              << "Zeros     |     Ones      |" << " \n";
        }
        log.printf ("       |%10i|%15i|%15i|%15i|\n", bloc,
                    mRefPIV [ref][bloc].size(), lmt0, lmt1);
      } // if we compute derivatives
    } // loop over the number of blocks
  } // loop over the number of references
  log << "\n";
  mDistancePIV.resize (mRefPIV.size());
}
    
//---------------------------------------------------------
// MAKE NEIGHBOR LISTS
// make neighbor list for all atoms, blocks and COM
//---------------------------------------------------------
void PIVCOM::makeNeighborLists (
  const bool cross,
  const bool direct,
  const double neighborCut, 
  const unsigned neighborStride,
  const std::vector <AtomNumber>& allAtomsList,
  const std::vector <std::vector <AtomNumber>>& pairList,
  const std::vector <std::vector <AtomNumber>>& comAtoms)
{
  log << "Creating Neighbor Lists \n";
  // WARNING: is neighborCut meaningful here?
  mBlockAtomsAll =
    std::make_unique <NeighborList> (
      allAtomsList,
      mPBC, getPbc(),
      neighborCut,
      neighborStride
    );

  for (unsigned com = 0; com < mPosCOM.size(); com++) {
    // WARNING: is neighborCut meaningful here?
    mBlockAtomCOM.emplace_back (
      std::make_unique <NeighborList> (
        comAtoms[com],
        mPBC, getPbc(),
        neighborCut,
        neighborStride
      )
    );
  }

  // Direct blocks AA, BB, CC, ...
  if (direct) {
    log << "Number of blocks: " << mBlockAtoms.size()
        << ", number of atom types: " << mAtomTypes << "\n";
    for (unsigned type = 0; type < mAtomTypes; type++) {
      mBlockAtoms.emplace_back (
        std::make_unique <NeighborList> (
          pairList[type],
          mPBC, getPbc(),
          neighborCut,
          neighborStride
        )
      );
    }
  }

  // Cross blocks AB, AC, BC, ...
  if (cross) {
    log << "Number of blocks: " << mBlockAtoms.size()
        << ", number of atom types: " << mAtomTypes << "\n";
    for (unsigned type1 = 0; type1 < mAtomTypes; type1++) {
      for (unsigned type2 = type1 + 1; type2 < mAtomTypes; type2++) {
        mBlockAtoms.emplace_back (
          std::make_unique <NeighborList> (
            pairList[type1], pairList[type2],
            false, mPBC, getPbc(),
            neighborCut,
            neighborStride
          )
        );
      }
    }
  }

  // Output neighborlist
  log << "Total neigbor lists: " << mBlocks << " \n";
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    log << "  list " << bloc + 1 << "   size " 
        << mBlockAtoms[bloc]->size() << " \n";
  }
}

//---------------------------------------------------------
// INITIALIZE SWITCHING FUNCTION
//---------------------------------------------------------
void PIVCOM::initializeSwitchingFunction (bool doScaleVolume)
{
  // Set switching function parameters
  log << "Switching Function Parameters \n";
  mSwitchFunc.resize(mBlocks);
  std::string errors;
  // setting the PIV reference for each atom pair blocks
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    if (mScaleVolume && doScaleVolume) {
      auto data = Tools::getWords (mSwitchData[bloc]);
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
      auto pos = mSwitchData[bloc].find ("R_0");
      mSwitchData[bloc].replace(
        pos + 4,
        sR0.size(),
        std::to_string (r0)
      );
    }
    mSwitchFunc[bloc].set(mSwitchData[bloc], errors);
    if (errors.length() != 0){
      error ("problem reading SWITCH"
             + std::to_string (bloc + 1)
             + " keyword : " + errors );
    }
    log << "  Swf: " << bloc + 1 << "  r0 = "
        << (mSwitchFunc[bloc].description()).c_str() << " \n";
  }
}
  

//---------------------------------------------------------
// INITIALIZE REFERENCE PIV
//---------------------------------------------------------
void PIVCOM::initializeReferencePIV ()
{
  for (unsigned ref = 0; ref < mRefPIV.size(); ref++) {
    log << "\n";
    log << "REFERENCE PDB # " << ref + 1 << " \n";
    //Transform and sort
    log << "Building Reference PIV Vector \n";
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      for (unsigned i = 0; i < mRefPIV[ref][bloc].size(); i++) {
        double df = 0.;
        mRefPIV[ref][bloc][i] = mSwitchFunc[bloc].calculate (mRefPIV[ref][bloc][i], df);
      }
      if (mDoSort[bloc]) {
        std::sort (mRefPIV[ref][bloc].begin(), mRefPIV[ref][bloc].end());
      }
      unsigned lmt0 = 0, lmt1 = 0;
      for (unsigned atm = 0; atm < mRefPIV[ref][bloc].size(); atm++) {
        if (mRefPIV[ref][bloc][atm] > 0.9) { lmt1++; }
        if (mRefPIV[ref][bloc][atm] < 0.1) { lmt0++; }
      }
      log.printf ("        |%10i|%15i|%15i|%15i|\n", bloc, mRefPIV[ref][bloc].size(), lmt0, lmt1);
    } // for each block
  } // for each reference file

  // we compute lambda if they are more than one reference
  auto lambda = double (0.);
  if (mRefPIV.size () > 1) {
    auto distance = distanceAB (mRefPIV[0], mRefPIV[1]);
    lambda = 2.3 / distance;
    log << "lambda=" << lambda << " d_12=" << distance << "\n";
  } else {
    log << "There is only one reference so we cannot compute lambda.\n";
  }
  auto pValLambda = getPntrToComponent ("lambda");
  pValLambda->set (lambda);
  log << "\n";
}

//---------------------------------------------------------
// INITIALIZE NEIGHBOR LIST 
//---------------------------------------------------------
void PIVCOM::initializeNeighborList ()
{
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    const auto& atomsOfBlock = mBlockAtoms[bloc]->getFullAtomList ();
    for (unsigned atm = 0; atm < atomsOfBlock.size(); atm++) {
      auto position = mPosCOM[atm];
      mPrevPosition[bloc].push_back (position);
    }
    mBlockAtoms[bloc]->update (mPrevPosition[bloc]);
  }
}

//---------------------------------------------------------
// UPDATE NEIGHBOR LIST 
//---------------------------------------------------------
void PIVCOM::updateNeighborList () 
{
  bool doUpdate = false;
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    auto updatedPos = std::vector <Vector> ();
    // get new position
    const auto& atomsOfBlock = mBlockAtoms[bloc]->getFullAtomList ();
    for (unsigned atm = 0; atm < atomsOfBlock.size(); atm++) {
      auto position = mPosCOM[atm];
      doUpdate = doUpdate ||
        pbcDistance (position, mPrevPosition[bloc][atm]).modulo2()
        >= mAtomsSkin * mAtomsSkin;
      updatedPos.push_back (position);
    }

    // update positions if needed
    if (doUpdate) {
      mPrevPosition[bloc] = std::move (updatedPos);
      mBlockAtoms[bloc]->update (mPrevPosition[bloc]);
      if (getStep() % 50000 == 0 && bloc == 0) {
        log << " Step " << getStep() << "  nl updated: "
            << mBlockAtoms[bloc]->size() << "\n";
      }
    }
  } // for each block
}

//---------------------------------------------------------
// PRINT CURRENT AND REFERENCE PIV (for test purpose) 
//---------------------------------------------------------
void PIVCOM::printCurrentAndReferencePIV (
  const std::vector <std::vector <int>>& atomIndex0,
  const std::vector <std::vector <int>>& atomIndex1,
  const std::vector <std::vector <double>>& currentPIV)
{
  for (unsigned ref = 0; ref < mRefPIV.size(); ref++) {
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      auto limit = unsigned (0);
      if (mDoSort[bloc]) {
        limit = currentPIV[bloc].size();
      } else {
        limit = mRefPIV[ref][bloc].size();
      }
      log.printf ("PIV Block:  %6i %12s %6i \n", bloc, "      Size:", limit);
      log.printf ("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV   ",
                 "    r-PIV   ","   i-j distance vector       ");
      for (unsigned i = 0; i < limit; i++) {
        auto i0 = unsigned (0);
        auto i1 = unsigned (0);
        if (mDoSort[bloc]) {
          i0 = atomIndex0[bloc][i];
          i1 = atomIndex1[bloc][i];
        } else {
          i0 = (mBlockAtoms[bloc]->getClosePairAtomNumber(i).first).index();
          i1 = (mBlockAtoms[bloc]->getClosePairAtomNumber(i).second).index();
        }
        auto position0 = atomPosition (i0);
        auto position1 = atomPosition (i1);
        auto distance = distanceAB (position0, position1);
        auto dfunc = double (0.);
        double cP, rP;
        if (mDoSort[bloc]) {
          cP = currentPIV[bloc][i];
          rP = mRefPIV[ref][bloc][mRefPIV[ref][bloc].size() - currentPIV[bloc].size() + i];
        } else {
          cP = mSwitchFunc[bloc].calculate(distance.modulo() * mVolumeFactor, dfunc);
          rP = mRefPIV[ref][bloc][i];
        }
        log.printf ("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
      }
    }
    log.printf ("This was a test, now exit \n");
    exit();
  }
}

//---------------------------------------------------------
// DISTANCE 
//---------------------------------------------------------
Vector PIVCOM::distanceAB (
  const Vector& A,
  const Vector& B)
{
  if (mPBC) {
    return pbcDistance (A, B);
  } else {
    return delta (A, B);
  }
}

double PIVCOM::distanceAB (
  const std::vector<std::vector<double>>& PIVA,
  const std::vector<std::vector<double>>& PIVB)
{
  auto distance = double (0.);
  for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
    auto sizeRef1 = PIVA [bloc].size ();
    auto sizeRef2 = PIVB [bloc].size ();
    auto maxSize = std::max (sizeRef1, sizeRef2);
    auto minSize = std::min (sizeRef1, sizeRef2);
    for (unsigned atm = 0; atm < maxSize; atm++) {
      auto fstRefValue = PIVA [bloc] [maxSize - minSize + atm];
      auto sndRefValue = PIVB [bloc] [atm];
      auto coord = fstRefValue - sndRefValue;
      distance += mScalingBlockFactor[bloc] * coord * coord;
    }
  }
  return distance;
}

//---------------------------------------------------------
// ATOM POSITION 
// return position for a given atom index
//---------------------------------------------------------
Vector PIVCOM::atomPosition (const unsigned atom)
{
  return mPosCOM[atom];
}

//---------------------------------------------------------
// CALCULATE
//---------------------------------------------------------
void PIVCOM::calculate()
{
  auto currentPIV = std::vector <std::vector <double>> (mBlocks);
  auto atmI0 = std::vector <std::vector <int>> (mBlocks);
  auto atmI1 = std::vector <std::vector <int>> (mBlocks);
  auto atmPrecI0 = std::vector <std::vector <int>> (mPrecision);
  auto atmPrecI1 = std::vector <std::vector <int>> (mPrecision);
  auto stride = int (1);
  auto rank = int (0);

  ///////////////////////////////////////////////////////////////////
  // get number of processors and current processor number
  if (!mSerial && comm.initialized ()) {
    stride = comm.Get_size();
    rank = comm.Get_rank();
  }
  
  ///////////////////////////////////////////////////////////////////
  // compute volume factor
  if (mScaleVolume) {
    mVolumeFactor = cbrt (mVolume0 / getBox().determinant());
  }

  ///////////////////////////////////////////////////////////////////
  // Transform (and sort) the rPIV before starting the dynamics
  // for first step neighborlist = actual position and build ref PIV
  if (mFirstStep || mComputeDerivatives) {
    initializeSwitchingFunction ();
    initializeReferencePIV ();
    initializeNeighborList ();
  } 

  ///////////////////////////////////////////////////////////////////
  // Decide whether to update lists based on atom displacement, every stride
  if (getStep() % mBlockAtomsAll->getStride() == 0) {
    updateNeighborList ();
  }

  ///////////////////////////////////////////////////////////////////
  // Do the sorting with a stride defined by updatePIV 
  if (getStep() % mUpdateStride == 0 || mFirstStep || mComputeDerivatives) {
    if (mComputeDerivatives) {
      log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV \n";
    }

    // build COMs from positions if requested
    if (mPBC) {
      makeWhole();
    }
    for (unsigned i = 0; i < mPosCOM.size(); i++) {
      mPosCOM[i][0] = 0.;
      mPosCOM[i][1] = 0.;
      mPosCOM[i][2] = 0.;
      for (const auto& atom : mBlockAtomCOM[i]->getFullAtomList()) {
        mPosCOM[i] += mMassFactor[atom.index()] * getPosition (atom.index());
      }
    }

    // Build "neigborlist" PIV blocks
    for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
      if (mDoSort[bloc]) {
        // Vectors collecting occupancies: orderVec one rank, orderVecAll all ranks
        auto orderVec = std::vector<int> (mPrecision, 0);
        currentPIV[bloc].resize(0);
        atmI0[bloc].resize(0);
        atmI1[bloc].resize(0);

        // Building distances for the PIV vector at time t
        if (mTimer) stopwatch.start("1 Build currentPIV");

        // If we have N cores and Na atoms, we need (Na/N + 1) process by cores
        for (unsigned atm = rank; atm < mBlockAtoms[bloc]->size(); atm += stride) {
          auto i0 = (mBlockAtoms[bloc]->getClosePairAtomNumber(atm).first).index();
          auto i1 = (mBlockAtoms[bloc]->getClosePairAtomNumber(atm).second).index();
          auto position0 = atomPosition (i0);
          auto position1 = atomPosition (i1);
          auto pairDist = distanceAB (position0, position1);
          auto df = double (0.);
          // transform distances with Switching function 
           auto vecInt = static_cast<int> (
             mSwitchFunc[bloc].calculate (
               pairDist.modulo() * mVolumeFactor, df)
             * static_cast<double> (mPrecision - 1) + 0.5
           );
          // tranform distance into int that serves as index for ordering vec
          orderVec[vecInt] += 1;
          // keep track of atom indices for force and virial calculations
          atmPrecI0[vecInt].push_back (i0);
          atmPrecI1[vecInt].push_back (i1);
        }
         
        if (mTimer) stopwatch.stop("1 Build currentPIV");
        if (mTimer) stopwatch.start("2 Sort currentPIV");

        if (!mSerial && comm.initialized ()) {
          // Vectors keeping track of the dimension and the starting-position 
          // of the rank-specific pair vector in the big pair vector.
          std::vector<int> vecDimension(stride, 0);
          std::vector<int> vecPos(stride, 0);
          // Vectors collecting occupancies: orderVec one rank, orderVecAll all ranks
          std::vector<int> orderVecAll(stride * mPrecision);
          // Big vectors containing all Atom indexes for every occupancy 
          // (atmI0O(Nprec,n) and atmI1O(Nprec,n) matrices in one vector)
          std::vector<int> atmI0F;
          std::vector<int> atmI1F;
          // Vector used to reconstruct arrays
          std::vector<int> counter(stride, 0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for (unsigned i = 1; i < mPrecision; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            atmI0F.insert(atmI0F.end(), atmPrecI0[i].begin(), atmPrecI0[i].end());
            atmI1F.insert(atmI1F.end(), atmPrecI1[i].begin(), atmPrecI1[i].end());
            atmPrecI0[i].resize(0);
            atmPrecI1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV block
          atmPrecI0[0].resize(0);
          atmPrecI1[0].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          orderVec[0] = 0;
          orderVec[mPrecision - 1] = 0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          unsigned dim = atmI0F.size();
          // vecDimension and vecPos keep track of the dimension and the starting-position 
          // of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim, 1, &vecDimension[0], 1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          unsigned finalDimension = 0;
          for (unsigned i = 1; i < stride; i++) {
            vecPos[i] = vecPos[i-1] + vecDimension[i-1];
            finalDimension += vecDimension[i];
          }
          finalDimension += vecDimension[0];

          // build big vectors for atom pairs on all ranks for all ranks
          std::vector<int> atmI0FinalAll (finalDimension);
          std::vector<int> atmI1FinalAll (finalDimension);
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          // Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV
          comm.Allgather(&orderVec[0], mPrecision, &orderVecAll[0], mPrecision);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&atmI0F[0], atmI0F.size(), &atmI0FinalAll[0], &vecDimension[0], &vecPos[0]);
          comm.Allgatherv(&atmI1F[0], atmI1F.size(), &atmI1FinalAll[0], &vecDimension[0], &vecPos[0]);

          if (mTimer) stopwatch.stop("2 Sort currentPIV");
          if (mTimer) stopwatch.start("3 Reconstruct currentPIV");

          // Reconstruct the full vectors from collections of Allgathered parts 
          // This is a tricky serial step, to assemble toghether PIV and 
          // atom-pair info from head-tail big vectors. Loop before on l and 
          // then on i would be better but the allgather should be modified
          for (unsigned i = 1; i < mPrecision; i++) {
            for (unsigned l = 0; l < stride; l++) {
              for (unsigned m = 0; m < orderVecAll[i + l * mPrecision]; m++) {
                currentPIV[bloc].push_back ( double(i) / double(mPrecision - 1) );
                atmI0[bloc].push_back ( atmI0FinalAll[counter[l] + vecPos[l]] );
                atmI1[bloc].push_back ( atmI1FinalAll[counter[l] + vecPos[l]] );
                counter[l] += 1;
              } // loop on the number of head-to-tail pieces
            } // loop on the ranks
          } // loop on the ordering vector excluding zero (i = 1)

          if (mTimer) stopwatch.stop("3 Reconstruct currentPIV");

        } else {
          for (unsigned i = 1; i < mPrecision; i++) {
            for (unsigned m = 0; m < orderVec[i]; m++) {
              currentPIV[bloc].push_back ( double(i) / double(mPrecision - 1) );
              atmI0[bloc].push_back (atmPrecI0[i][m]);
              atmI1[bloc].push_back (atmPrecI1[i][m]);
            }
          }
        } // if serial or parallel
      } // if we sort the PIV
    } // for each block

    // compute PIV-PIV distance and derivatives for each reference
    for (unsigned ref = 0; ref < mRefPIV.size(); ref++) {
      mDistancePIV[ref] = 0.;
    }
  } // if step % stride == 0

  ///////////////////////////////////////////////////////////////////
  // This test may be run by specifying the TEST keyword as input, it pritnts referencePIV and currentPIV and quits
  if (mTest) {
    printCurrentAndReferencePIV (
      atmI0,
      atmI1,
      currentPIV
    );
  }
  
  if (mTimer) stopwatch.start("4 Build For Derivatives");

  ///////////////////////////////////////////////////////////////////
  // compute derivatives 
  if (getStep() % mUpdateStride == 0) {
    for (unsigned ref = 0; ref < mRefPIV.size(); ref++) {
      // set to zero PIVdistance, derivatives and virial when they are calculated
      for (unsigned j = 0; j < mDerivatives.size(); j++) {
        for (unsigned k = 0; k < 3; k++) {
          mDerivatives[j][k] = 0.;
        }
      }
      for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < 3; j++) {
          mVirial[i][j] = 0.;
        }
      }
      // compute PIV-PIV distance and derivatives for each reference
      auto dm = double (0.);
      auto tPIV = double (0.);
      mDistancePIV[ref] = 0.;

      // Re-compute atomic distances for derivatives and compute PIV-PIV distance
      for (unsigned bloc = 0; bloc < mBlocks; bloc++) {
        unsigned limit = 0;
        if (mDoSort[bloc]) {
          limit = currentPIV[bloc].size();
        } else {
          limit = mRefPIV[ref][bloc].size();
        }

        // for each atm
        for (unsigned atm = rank; atm < limit; atm += stride) {
          auto i0 = unsigned (0);
          auto i1 = unsigned (0);
          // recompute PIV only once
          if (mDoSort[bloc]) {
            i0 = atmI0[bloc][atm];
            i1 = atmI1[bloc][atm];
          } else {
            i0 = (mBlockAtoms[bloc]->getClosePairAtomNumber(atm).first).index();
            i1 = (mBlockAtoms[bloc]->getClosePairAtomNumber(atm).second).index();
          }
          auto position0 = atomPosition (i0);
          auto position1 = atomPosition (i1);
          auto distance = distanceAB (position0, position1);
          auto dfunc = double (0.);
          // this is needed for dfunc and dervatives
          dm = distance.modulo();
          tPIV = mSwitchFunc[bloc].calculate (dm * mVolumeFactor, dfunc);

          // PIV distance
          double coord = 0.;
          if (!mDoSort[bloc] || mComputeDerivatives) {
            coord = tPIV - mRefPIV [ref] [bloc] [atm];
          } else {
            auto refSize = mRefPIV [ref] [bloc].size ();
            auto currentSize = currentPIV [bloc].size ();
            auto sizeDifference = refSize - currentSize;
            auto reference = mRefPIV [ref] [bloc] [sizeDifference + atm];
            auto current = currentPIV [bloc] [atm];
            coord = current - reference;
          }

          // Calculate derivatives, virial, and variable = sum_j (scaling[j] *(currentPIV-rPIV)_j^2)
          // WARNING: dfunc=dswf/(mVolumeFactor * dm)  (this may change in future Plumed versions)
          auto tmp = 2. * mScalingBlockFactor[bloc] * coord
                       * mVolumeFactor*mVolumeFactor * dfunc;
          auto tempDeriv = tmp * distance;
          // 0.5 * (x_i - x_k) * f_ik         (force on atom k due to atom i)
          const auto& atomsBlock0 = mBlockAtomCOM [i0]->getFullAtomList ();
          const auto& atomsBlock1 = mBlockAtomCOM [i1]->getFullAtomList ();
          // loop on first atom of each pair
          for (const auto& atom0 : atomsBlock0) {
            auto x0 = atom0.index();
            auto dist = Vector (0, 0, 0);
            auto position0 = getPosition (x0);
            mDerivatives[x0] -= tempDeriv * mMassFactor[x0];
            for (const auto& atom1 : atomsBlock0) {
              dist += distanceAB (position0, getPosition (atom1.index()) );
            }
            for (const auto& atom1 : atomsBlock1) {
              dist += distanceAB (position0, getPosition (atom1.index()) );
            }
            mVirial -= 0.25 * mMassFactor[x0] * Tensor (dist, tempDeriv);
          }

          // loop on second atom of each pair
          for (const auto& atom1 : atomsBlock1) {
            auto x1 = atom1.index();
            auto dist = Vector (0, 0, 0);
            auto position1 = getPosition (x1);
            mDerivatives[x1] += tempDeriv * mMassFactor[x1];
            for (const auto& atom0 : atomsBlock1) {
              dist += distanceAB (position1, getPosition (atom0.index()) );
            }
            for (const auto& atom0 : atomsBlock0) {
              dist += distanceAB (position1, getPosition (atom0.index()) );
            }
            mVirial += 0.25 * mMassFactor[x1] * Tensor (dist, tempDeriv);
          }

          if (mScaleVolume) {
            mVirial += 1./3. * tmp * dm*dm * Tensor::identity();
          }
          mDistancePIV[ref] += mScalingBlockFactor[bloc] * coord*coord;
        } // loop on atoms-atoms pairs
      } // loop on block

      if (!mSerial && comm.initialized ()) {
        comm.Barrier();
        comm.Sum(&mDistancePIV[ref], 1);
        if (!mDerivatives.empty()) {
          comm.Sum (&mDerivatives[0][0], 3 * mDerivatives.size());
        }
        comm.Sum (&mVirial[0][0], 9);
      }
    } // loop on the number of references
  } // if update_piv
  
  //Timing
  if (mTimer) stopwatch.stop ("4 Build For Derivatives");
  if (mTimer) {
    log.printf ("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log << stopwatch;
  }

  // Update derivatives, virial
  setBoxDerivatives(mVirial);
  for (unsigned i = 0; i < mDerivatives.size(); ++i) {
    setAtomsDerivatives (i, mDerivatives[i]);
  }
  // update variables
  for (unsigned ref = 0; ref < mRefPIV.size(); ref++) {
    auto pValDistance = getPntrToComponent ("d" + std::to_string (ref + 1));
    pValDistance->set (mDistancePIV[ref]);
  }
  mFirstStep = false;
} // end of calculate
} // close namespace piv
} // close namespace PLMD
