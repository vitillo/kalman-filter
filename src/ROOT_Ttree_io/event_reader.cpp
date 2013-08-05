//==================================================================================================
//
// Name         : EventReader.cpp
// Description  : Class to parse the test ROOT file.
//
//==================================================================================================

#include "ROOT_Ttree_io/event_reader.h"
// #if USE_EIGEN_MATRIX_INVERSION
// #include "common_eigen.h"
// #endif
// Constructor
EventReader::EventReader(void)
{
  // Initialize pointers with NULL
  mLocX       = nullptr;
  mLocY       = nullptr;
  mErr_locX   = nullptr;
  mErr_locY   = nullptr;
  mCov_locXY  = nullptr;
#if HAS_MATERIAL_EFFECTS
  mSigmaDeltaTheta = 0;
  mSigmaDeltaPhi = 0;
  mSigmaDeltaE = 0;
#endif // #if HAS_MATERIAL_EFFECTS

  mDetType    = nullptr;
  mBec        = nullptr;
  mRefLocX    = nullptr;
  mRefLocY    = nullptr;
  mRefPhi     = nullptr;
  mRefTheta   = nullptr;
  mRefQoverP  = nullptr;
  mJacobi     = nullptr;
  mD0			= nullptr;
  mZ0			= nullptr;
  mPhi        = nullptr;
  mTheta      = nullptr;
  mQoverP     = nullptr;

  mMCI		= nullptr;
  mTTTI		= nullptr;
  mTTD0		= nullptr;
  mTTZ0		= nullptr;
  mTTQoverP	= nullptr;
  mTTPhi		= nullptr;
  mTTTheta	= nullptr;
#if HAS_COMPETING_ROTS
  mCompetingRotIds = 0;
  mRotLocX = nullptr;
  mRotLocY = nullptr;
  mRotProbs = nullptr;
  mRotIds = nullptr;

  mRotCov00 = 0;
  mRotCov01 = 0;
  mRotCov10 = 0;
  mRotCov11 = 0;
#endif

}

// Destructor
EventReader::~EventReader(void)
{
  if (mJacobi != nullptr)
    delete mJacobi;
}

// Reset class members representing different track attributes
  void
EventReader::resetAttributeVectors (void)
{
  mLocX       = nullptr;
  mLocY       = nullptr;
  mErr_locX   = nullptr;
  mErr_locY   = nullptr;
  mCov_locXY  = nullptr;
#if HAS_MATERIAL_EFFECTS
  mSigmaDeltaTheta = 0;
  mSigmaDeltaPhi = 0;
  mSigmaDeltaE = 0;
#endif // #if HAS_MATERIAL_EFFECTS
  mDetType    = nullptr;
  mBec        = nullptr;
  mRefLocX    = nullptr;
  mRefLocY    = nullptr;
  mRefPhi     = nullptr;
  mRefTheta   = nullptr;
  mRefQoverP  = nullptr;
  mD0			= nullptr;
  mZ0			= nullptr;
  mPhi        = nullptr;
  mTheta      = nullptr;
  mQoverP     = nullptr;

  if (mJacobi != nullptr)
    delete mJacobi;

  mJacobi = new std::vector<JacobiVal_t*>(ORDER * ORDER);

  for (size_t i = 0; i < mJacobi->size ( ); ++i)
    mJacobi->at (i) = nullptr;

  mMCI		= nullptr;
  mTTTI		= nullptr;
  mTTD0		= nullptr;
  mTTZ0		= nullptr;
  mTTQoverP	= nullptr;
  mTTPhi		= nullptr;
  mTTTheta	= nullptr;
#if HAS_COMPETING_ROTS
  mCompetingRotIds = 0;
  mRotLocX = nullptr;
  mRotLocY = nullptr;
  mRotIds = 0;
  mRotProbs = 0;
  mRotCov00 = 0;
  mRotCov01 = 0;
  mRotCov10 = 0;
  mRotCov11 = 0;
#endif

  return;
}

// Setup branch addresses for each track attribute
  int
EventReader::loadBranches (TTree *tree, const char* prefix)
{
  const std::string bNamePrefix = prefix;

  std::string bName;

  bName = bNamePrefix + kLocX;
  if (tree->SetBranchAddress (bName.c_str ( ), &mLocX) != 0) {
    PRINT_ERR ("SetBranchAddress failed for locX");
    return -1;
  }

  bName = bNamePrefix + kLocY;
  if (tree->SetBranchAddress (bName.c_str ( ), &mLocY) != 0) {
    PRINT_ERR ("SetBranchAddress failed for locY");
    return -1;
  }

  bName = bNamePrefix + kErrLocX;
  if (tree->SetBranchAddress (bName.c_str ( ), &mErr_locX) != 0) {
    PRINT_ERR ("SetBranchAddress failed for err_locX");
    return -1;
  }

  bName = bNamePrefix + kErrLocY;
  if (tree->SetBranchAddress (bName.c_str ( ), &mErr_locY) != 0) {
    PRINT_ERR ("SetBranchAddress failed for err_locY");
    return -1;
  }

  bName = bNamePrefix + kCovLocXY;
  if (tree->SetBranchAddress(bName.c_str(), &mCov_locXY) != 0) {
    PRINT_ERR("SetBranchAddress failed for cov_locXY");
    return -1;
  }
#if HAS_MATERIAL_EFFECTS
  bName = bNamePrefix + kSigmaDeltaPhiScat;
  if (tree->SetBranchAddress(bName.c_str(), &mSigmaDeltaPhi) != 0) {
    PRINT_ERR("SetBranchAddress failed for SigmaDeltaPhiScat");
    return -1;
  }
  bName = bNamePrefix + kSigmaDeltaThetaScat;
  if (tree->SetBranchAddress(bName.c_str(), &mSigmaDeltaTheta) != 0) {
    PRINT_ERR("SetBranchAddress failed for SigmaDeltaThetaScat");
    return -1;
  }
  bName = bNamePrefix + kSigmaDeltaEloss;
  if (tree->SetBranchAddress(bName.c_str(), &mSigmaDeltaE) != 0) {
    PRINT_ERR("SetBranchAddress failed for SigmaDeltaEloss");
    return -1;
  }
#endif // #if HAS_MATERIAL_EFFECTS

  bName = bNamePrefix + kDetType;
  if (tree->SetBranchAddress (bName.c_str ( ), &mDetType) != 0) {
    PRINT_ERR ("SetBranchAddress failed for detType");
    return -1;
  }

  bName = bNamePrefix + kBEC;
  if (tree->SetBranchAddress (bName.c_str ( ), &mBec) != 0) {
    PRINT_ERR ("SetBranchAddress failed for bec");
    return -1;
  }

  bName = bNamePrefix + kRefLocX;
  if (tree->SetBranchAddress (bName.c_str ( ), &mRefLocX) != 0) {
    PRINT_ERR ("SetBranchAddress failed for refLocX");
    return -1;
  }

  bName = bNamePrefix + kRefLocY;
  if (tree->SetBranchAddress (bName.c_str ( ), &mRefLocY) != 0) {
    PRINT_ERR ("SetBranchAddress failed for refLocY");
    return -1;
  }

  bName = bNamePrefix + kRefPhi;
  if (tree->SetBranchAddress (bName.c_str ( ), &mRefPhi) != 0) {
    PRINT_ERR ("SetBranchAddress failed for refPhi");
    return -1;
  }

  bName = bNamePrefix + kRefTheta;
  if (tree->SetBranchAddress (bName.c_str ( ), &mRefTheta) != 0) {
    PRINT_ERR ("SetBranchAddress failed for refTheta");
    return -1;
  }

  bName = bNamePrefix + kRefQoverP;
  if (tree->SetBranchAddress (bName.c_str ( ), &mRefQoverP)
      != 0) {
    PRINT_ERR ("SetBranchAddress failed for refQoverP");
    return -1;
  }

  bName = bNamePrefix + kQoverP;
  if (tree->SetBranchAddress (bName.c_str ( ), &mQoverP) != 0) {
    PRINT_ERR ("SetBranchAddress failed for qoverp");
    return -1;
  }

  bName = bNamePrefix + kD0;
  if (tree->SetBranchAddress (bName.c_str ( ), &mD0) != 0) {
    PRINT_ERR ("SetBranchAddress failed for d0");
    return -1;
  }

  bName = bNamePrefix + kZ0;
  if (tree->SetBranchAddress (bName.c_str ( ), &mZ0) != 0) {
    PRINT_ERR ("SetBranchAddress failed for z0");
    return -1;
  }

  bName = bNamePrefix + kPhi;
  if (tree->SetBranchAddress (bName.c_str ( ), &mPhi) != 0) {
    PRINT_ERR ("SetBranchAddress failed for phi");
    return -1;
  }

  bName = bNamePrefix + kTheta;
  if (tree->SetBranchAddress (bName.c_str ( ), &mTheta) != 0) {
    PRINT_ERR ("SetBranchAddress failed for theta");
    return -1;
  }

  // Jacobian Matrix elements row 0
  bName = bNamePrefix + kTransportJac00;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (0))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac00");
    return -1;
  }

  bName = bNamePrefix + kTransportJac01;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (1))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac01");
    return -1;
  }

  bName = bNamePrefix + kTransportJac02;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (2))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac02");
    return -1;
  }

  bName = bNamePrefix + kTransportJac03;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (3))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac03");
    return -1;
  }

  bName = bNamePrefix + kTransportJac04;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (4))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac04");
    return -1;
  }

  // Jacobian Matrix elements row 1
  bName = bNamePrefix + kTransportJac10;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (5))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac10");
    return -1;
  }

  bName = bNamePrefix + kTransportJac11;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (6))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac11");
    return -1;
  }

  bName = bNamePrefix + kTransportJac12;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (7))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac12");
    return -1;
  }

  bName = bNamePrefix + kTransportJac13;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (8))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac13");
    return -1;
  }

  bName = bNamePrefix + kTransportJac14;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (9))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac14");
    return -1;
  }

  // Jacobian Matrix elements row 2
  bName = bNamePrefix + kTransportJac20;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (10))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac20");
    return -1;
  }

  bName = bNamePrefix + kTransportJac21;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (11))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac21");
    return -1;
  }

  bName = bNamePrefix + kTransportJac22;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (12))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac22");
    return -1;
  }

  bName = bNamePrefix + kTransportJac23;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (13))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac23");
    return -1;
  }

  bName = bNamePrefix + kTransportJac24;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (14))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac24");
    return -1;
  }

  // Jacobian Matrix elements row 3
  bName = bNamePrefix + kTransportJac30;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (15))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac30");
    return -1;
  }

  bName = bNamePrefix + kTransportJac31;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (16))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac31");
    return -1;
  }

  bName = bNamePrefix + kTransportJac32;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (17))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac32");
    return -1;
  }

  bName = bNamePrefix + kTransportJac33;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (18))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac33");
    return -1;
  }

  bName = bNamePrefix + kTransportJac34;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (19))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac34");
    return -1;
  }

  // Jacobian Matrix elements row 4
  bName = bNamePrefix + kTransportJac40;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (20))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac40");
    return -1;
  }

  bName = bNamePrefix + kTransportJac41;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (21))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac41");
    return -1;
  }

  bName = bNamePrefix + kTransportJac42;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (22))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac42");
    return -1;
  }

  bName = bNamePrefix + kTransportJac43;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (23))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac43");
    return -1;
  }

  bName = bNamePrefix + kTransportJac44;
  if (tree->SetBranchAddress (bName.c_str ( ), &(mJacobi->at (24))) != 0) {
    PRINT_ERR ("SetBranchAddress failed for Jac44");
    return -1;
  }

  bName = bNamePrefix + kMci;
  if (tree->SetBranchAddress(bName.c_str(), &(mMCI)) != 0) {
    PRINT_ERR("SetBranchAddress failed for MCI");
    return -1;
  }

  bName = kMcPerigeeOk;
  if (tree->SetBranchAddress(bName.c_str(), &(mMCPerigeeOk)) != 0) {
    PRINT_ERR("SetBranchAddress failed for MCPerigeeOk");
    return -1;
  }

  bName = kTti;
  if (tree->SetBranchAddress(bName.c_str(), &(mTTTI)) != 0) {
    PRINT_ERR("SetBranchAddress failed for TTI");
    return -1;
  }

  bName = kTti;
  if (tree->SetBranchAddress(bName.c_str(), &(mTTTI)) != 0) {
    PRINT_ERR("SetBranchAddress failed for TTI");
    return -1;
  }

  //bName = "truthtrack_" + kQoverP;
  bName = "mc_perigee_" + kQoverP;
  if (tree->SetBranchAddress(bName.c_str(), &mTTQoverP) != 0) {
    PRINT_ERR("SetBranchAddress failed for qoverp");
    return -1;
  }

  //bName = "truthtrack_" + kD0;
  bName = "mc_perigee_" + kD0;
  if (tree->SetBranchAddress(bName.c_str(), &mTTD0) != 0) {
    PRINT_ERR("SetBranchAddress failed for d0");
    return -1;
  }

  //bName = "truthtrack_" + kZ0;
  bName = "mc_perigee_" + kZ0;
  if (tree->SetBranchAddress(bName.c_str(), &mTTZ0) != 0) {
    PRINT_ERR("SetBranchAddress failed for z0");
    return -1;
  }

  //bName = "truthtrack_" + kPhi;
  bName = "mc_perigee_" + kPhi;
  if (tree->SetBranchAddress(bName.c_str(), &mTTPhi) != 0) {
    PRINT_ERR("SetBranchAddress failed for phi");
    return -1;
  }

  //bName = "truthtrack_" + kTheta;
  bName = "mc_perigee_" + kTheta;
  if (tree->SetBranchAddress(bName.c_str(), &mTTTheta) != 0) {
    PRINT_ERR("SetBranchAddress failed for theta");
    return -1;
  }

  //rots
#if HAS_COMPETING_ROTS
  //
  bName = bNamePrefix + kCompetingRotIds;
  //
  if (tree->SetBranchAddress(bName.c_str(), &mCompetingRotIds) != 0) {
    PRINT_ERR("SetBranchAddress failed for mCompetingRotIds");
    return -1;
  }

  bName = bNamePrefix + kRotLocX;

  if (tree->SetBranchAddress(bName.c_str(), &mRotLocX) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotLocX");
    return -1;
  }
  //    std::cout << bName << " " << mRotLocX->size() << std::endl;
  bName = bNamePrefix + kRotLocY;
  if (tree->SetBranchAddress(bName.c_str(), &mRotLocY) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotLocY");
    return -1;
  }

  bName = bNamePrefix + kRotIds;
  if (tree->SetBranchAddress(bName.c_str(), &mRotIds) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotIds");
    return -1;
  }

  bName = bNamePrefix + kRotProbs;
  if (tree->SetBranchAddress(bName.c_str(), &mRotProbs) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotProbs");
    return -1;
  }
  bName = bNamePrefix + kRotCov00;
  if (tree->SetBranchAddress(bName.c_str(), &mRotCov00) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotCov00 = 0;");
    return -1;
  }
  bName = bNamePrefix + kRotCov01;
  if (tree->SetBranchAddress(bName.c_str(), &mRotCov01) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotCov01 = 0;");
    return -1;
  }
  bName = bNamePrefix + kRotCov10;
  if (tree->SetBranchAddress(bName.c_str(), &mRotCov10) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotCov10 = 0;");
    return -1;
  }
  bName = bNamePrefix + kRotCov11;
  if (tree->SetBranchAddress(bName.c_str(), &mRotCov11) != 0) {
    PRINT_ERR("SetBranchAddress failed for mRotCov11 = 0;");
    return -1;
  }
#endif
  //end rots

  return 0;
}

// Simple event sanity check. Test root dataset for matching event and track sizes
  int
EventReader::eventSanityCheck (void)
{
  size_t nTracks, nElements;

  if (!mLocX || !mLocY || !mErr_locX || !mErr_locY || !mCov_locXY || !mDetType
      || !mBec || !mRefLocX || !mRefLocY || !mRefPhi || !mRefTheta || !mRefQoverP
      || !mJacobi || !mD0 || !mZ0 ||!mPhi || !mTheta || !mQoverP) {

    PRINT_ERR ("Branch vector is null!");
    return -1;
  }
#if HAS_COMPETING_ROTS
  // !mCompetingRotIds ||
  if (!mRotLocX || !mRotLocY || !mRotIds || !mRotProbs) {

    PRINT_ERR("Branch vector of rots is null!");
    return -1;
  }
#endif
  for (size_t i = 0; i < mJacobi->size ( ); ++i) {
    if (!mJacobi->at (i)) {
      PRINT_ERR ("Jacobi matrix branch is null!");
      return -1;
    }
  }

  // Sanity check: If the vector size is not equal for all attributes, the data file
  // cannot be a valid data set.
  nTracks = mLocX->size ( );

  if (mLocY->size ( ) != nTracks || mErr_locX->size ( ) != nTracks
      || mErr_locY->size ( )  != nTracks || mCov_locXY->size ( ) != nTracks
      || mDetType->size ( )   != nTracks || mBec->size ( )       != nTracks
      || mRefLocX->size ( )   != nTracks || mRefLocY->size ( )   != nTracks
      || mRefPhi->size ( )    != nTracks || mRefTheta->size ( )  != nTracks
      || mRefQoverP->size ( ) != nTracks || mPhi->size ( )       != nTracks
      || mD0->size ( )        != nTracks || mZ0->size ( )        != nTracks
      || mTheta->size ( )     != nTracks || mQoverP->size ( )    != nTracks) {

    PRINT_ERR ("Number of tracks is not equal for all attributes!");
    return -1;
  }

  for (size_t i = 0; i < mJacobi->size ( ); ++i) {
    if (mJacobi->at (i)->size ( ) != nTracks) {
      PRINT_ERR ("Number of tracks is not equal for all attributes!");
      return -1;
    }
  }

  for (size_t track = 0; track < nTracks; ++track) {
    nElements = (*mLocX)[track].size ( );

    if ((*mLocY)[track].size ( ) != nElements
	|| (*mErr_locX)[track].size ( )  != nElements
	|| (*mErr_locY)[track].size ( )  != nElements
	|| (*mCov_locXY)[track].size ( ) != nElements
	|| (*mDetType)[track].size ( )   != nElements
	|| (*mBec)[track].size ( )       != nElements
	|| (*mRefLocX)[track].size ( )   != nElements
	|| (*mRefLocY)[track].size ( )   != nElements
	|| (*mRefPhi)[track].size ( )    != nElements
	|| (*mRefTheta)[track].size ( )  != nElements
	|| (*mRefQoverP)[track].size ( ) != nElements) {

      PRINT_ERR ("Hit count is not equal for all attributes!");
      return -1;
    }

    for (size_t i = 0; i < mJacobi->size ( ); ++i) {
      if (mJacobi->at (i)->at (track).size ( ) != nElements) {
	PRINT_ERR ("Number of Hits is not equal for all attributes!");
	return -1;
      }
    }
  }
#if HAS_COMPETING_ROTS
  for (size_t track = 0; track < nTracks; ++track) {
    nElements = (*mLocX)[track].size();
    //
    for (size_t cRot = 0; cRot < nElements; ++cRot) {
      size_t nRots = (*mRotLocX)[track][cRot].size();
      if (nRots != (*mRotIds)[track][cRot].size()
	 ) {
	//
	PRINT_ERR("Hit count is not equal for all attributes!");
	return -1;
      }

    }
  }
#endif
  return 0;
}

// Serialize the Jacobian matrix values into a simple float array matrix
  int
EventReader::setJacobiMatrix (scalar_t ( &matrix )[ORDER * ORDER],
    unsigned int trackIdx, unsigned int hitIdx)
{
  if (!matrix) {
    PRINT_ERR ("matrix reference is NULL!");
    return -1;
  }

  if (!mJacobi) {
    PRINT_ERR ("mJacobi is NULL!");
    return -1;
  }

  for (size_t i = 0; i < mJacobi->size ( ); ++i) {
    if (!mJacobi->at (i)) {
      PRINT_ERR ("mJacobi branch is NULL!");
      return -1;
    }

    if (trackIdx >= mJacobi->at (i)->size ( )) {
      PRINT_ERR ("trackIdx is out of bounds!");
      return -1;
    }

    if (hitIdx >= mJacobi->at (i)->at (trackIdx).size ( )) {
      PRINT_ERR ("hitIdx is out of bounds!");
      return -1;
    }

    matrix[i] = mJacobi->at (i)->at (trackIdx).at (hitIdx);
  }

  return 0;
}

// Parse the root file located at path and parse the Ttree treeName.
// prefix identifies the track which is to be used, e.g.
// trk_, trkSiSpSeeded_, pseudoTrk_
  int
EventReader::parse (const char* path, const char* treeName, const char* prefix)
{
  TFile *file;
  TTree *tree;

  // Open data file and read TTree
  file = new TFile (path, "READ", "TrackFitData Access");

  if (!file) {
    std::cerr << "Could not open " << path << std::endl;
    return -1;
  }

  tree = ( TTree* )file->Get (treeName);

  if (!tree) {
    std::cerr << "Could not get tree " << treeName << std::endl;
    delete file;
    return -1;
  }

  // Read events and their tracks. Each data attribute is represented
  // by a branch in the TTree.
  resetAttributeVectors ( );

  if (loadBranches (tree, prefix) == -1) {
    PRINT_ERR ("loadBranches failed!");
    resetAttributeVectors ( );
    delete file;
    return -1;
  }

  mEvents.clear ( );

  // Loop over all events
  for (Long64_t events = 0; events < tree->GetEntriesFast ( ); ++events) {
    Long64_t local_events = tree->LoadTree (events);

    if (local_events < 0) {
      PRINT_ERR ("LoadTree returned value < 0");
      resetAttributeVectors ( );
      delete file;
      return -1;
    }

    // Load branches into the matching vectors
    if (tree->GetEntry (events) <= 0) {
      PRINT_ERR ("tree->GetEntry failed!");
      resetAttributeVectors ( );
      delete file;
      return -1;
    }

    if (eventSanityCheck ( ) == -1) {
      PRINT_ERR ("Event sanity check failed!");
      resetAttributeVectors ( );
      delete file;
      return -1;
    }

    KF_Event_t newEvent;
    Track_t newTrack;
    TrackHit_t trackHit;
#if HAS_COMPETING_ROTS
    for (int i = 0; i < MAXROTS; ++i) {
      trackHit.rotIds[i] = 0;
    }
#endif
    newEvent.clear();

    // Loop over all tracks
    for (unsigned int track = 0; track < mLocX->size ( ); ++track) {
      newTrack.track.clear ( );

      // phi, theta, qoverp are values referring to a complete track
      newTrack.info.d0     = (*mD0)[track];
      newTrack.info.z0     = (*mZ0)[track];
      newTrack.info.phi    = (*mPhi)[track];
      newTrack.info.theta  = (*mTheta)[track];
      newTrack.info.qoverp = (*mQoverP)[track];

      // FIXME: has to check whether MC index is valid!
      //int tti = (*mTTTI)[ (*mMCI)[track] ];
      // FIXME: has to check whether truth track index is valid!
      //std::cout << "Track:" << track <<  " TruthTrackIndex: " << tti <<std::endl;
      newTrack.truthTrackInfo.d0 = -999.;
      newTrack.truthTrackInfo.z0 = -999.;
      newTrack.truthTrackInfo.phi = -999.;
      newTrack.truthTrackInfo.theta = -999.;
      newTrack.truthTrackInfo.qoverp = -999.;
      int mcindex = (*mMCI)[track];
      if (mcindex > -1) {
	if (mMCPerigeeOk->at(mcindex)) {
	  newTrack.truthTrackInfo.d0 = (*mTTD0)[mcindex];
	  newTrack.truthTrackInfo.z0 = (*mTTZ0)[mcindex];
	  newTrack.truthTrackInfo.phi = (*mTTPhi)[mcindex];
	  newTrack.truthTrackInfo.theta = (*mTTTheta)[mcindex];
	  newTrack.truthTrackInfo.qoverp = (*mTTQoverP)[mcindex];
	}
      }
      bool validTrack = true;
      //             std::cout<< "track " << track << std::endl;
      // Loop over all hits
      for (unsigned int hit = 0; hit < (*mLocX)[track].size ( ); ++hit) {
	// No NULL-pointer check, since eventSanityCheck passed
	trackHit.normal[0]  = (*mLocX)[track][hit];
	trackHit.normal[1]  = (*mLocY)[track][hit];
	trackHit.err_locX   = (*mErr_locX)[track][hit];
	trackHit.err_locY   = (*mErr_locY)[track][hit];
	trackHit.cov_locXY  = (*mCov_locXY)[track][hit];
	trackHit.detType    = (*mDetType)[track][hit];
	trackHit.bec        = (*mBec)[track][hit];
	trackHit.ref[0]     = (*mRefLocX)[track][hit];
	trackHit.ref[1]     = (*mRefLocY)[track][hit];
	trackHit.ref[2]     = (*mRefPhi)[track][hit];
	trackHit.ref[3]     = (*mRefTheta)[track][hit];
	trackHit.ref[4]     = (*mRefQoverP)[track][hit];
#if HAS_MATERIAL_EFFECTS
	trackHit.sigmaDeltaThetaSquared = pow(mSigmaDeltaTheta->at(track).at(hit), 2);
	trackHit.sigmaDeltaPhiSquared = pow(mSigmaDeltaPhi->at(track).at(hit), 2);
	//trackHit.sigmaDeltaQoverPSquared = (trackHit.ref[4]*trackHit.ref[4]>0.) ? ( mSigmaDeltaE->at(track).at(hit) * sqrt(PionMassSquared + 1./(trackHit.ref[4]*trackHit.ref[4])) * fabs(pow(trackHit.ref[4],3)) ) : 0.;
	const double qOverPSquared=trackHit.ref[4]*trackHit.ref[4];
	trackHit.sigmaDeltaQoverPSquared = (qOverPSquared>0.) ? ( pow(mSigmaDeltaE->at(track).at(hit),2) * (PionMassSquared + 1./qOverPSquared) * pow(qOverPSquared,3) ) : 0.;
#endif // #if HAS_MATERIAL_EFFECTS
#if HAS_COMPETING_ROTS
	trackHit.competingRotIds = (*mCompetingRotIds)[track][hit];
	//std::cout << "read compRot " << trackHit.competingRotIds << " track " << track << " hit " << hit << std::endl;
	for (unsigned int rot = 0; rot < (*mRotLocX)[track][hit].size(); ++rot) {
	  trackHit.rotLocX[rot] = (*mRotLocX)[track][hit][rot];
	  trackHit.rotLocY[rot] = (*mRotLocY)[track][hit][rot];
	  trackHit.rotIds[rot] = (*mRotIds)[track][hit][rot];
	  trackHit.rotProbs[rot] = (*mRotProbs)[track][hit][rot];

	  trackHit.rotCov[rot][0] = (*mRotCov00)[track][hit][rot];
	  trackHit.rotCov[rot][1] = (*mRotCov01)[track][hit][rot];
	  trackHit.rotCov[rot][2] = (*mRotCov10)[track][hit][rot];
	  trackHit.rotCov[rot][3] = (*mRotCov11)[track][hit][rot];
	  //                    if (DBG_LVL == 0) if ((*mRotLocY)[track][hit].size() > 1) std::cout << rot << "/" << (*mRotLocY)[track][hit].size() - 1 << " id " << trackHit.rotIds[rot] << " locY " << trackHit.rotLocY[rot] << " effY " << trackHit.measurement[1] << " cov " << trackHit.rotCov[rot][1] << std::endl;
	  //                    std::cout << "\trot " << trackHit.rotIds[rot] << std::endl;
	}
#endif
	if (abs(trackHit.normal[1] + 99999.) > 1.E-4) {
	  trackHit.is2Dim = 1;
	} else {
	  trackHit.is2Dim = 0;
	}
	// no measurement at all (pure material surface):
	if (!(fabs(trackHit.normal[0] + 99999.) > 1.E-4)) {
	  trackHit.is2Dim = 3;
	  //                     std::cout<< "mat lay: sigmaDeltaPhi=" << mSigmaDeltaPhi->at(track).at(hit)
	  //                              << " sigmaDeltaTheta=" << mSigmaDeltaTheta->at(track).at(hit)
	  //                              << " sigmaDeltaE=" << mSigmaDeltaE->at(track).at(hit) << std::endl;
	}// else {
	//                     std::cout<< "meas lay: sigmaDeltaPhi=" << mSigmaDeltaPhi->at(track).at(hit)
	//                              << " sigmaDeltaTheta=" << mSigmaDeltaTheta->at(track).at(hit)
	//                              << " sigmaDeltaE=" << mSigmaDeltaE->at(track).at(hit) << std::endl;
	//                 }
	//Check if this track has valid Data
	//if (abs(trackHit.normal[0] + 999.) < 1.E-4 || abs(trackHit.normal[1] + 999.) < 1.E-4 || abs(trackHit.ref[0] + 999.) < 1.E-4 || abs(trackHit.ref[1] + 999.) < 1.E-4) {
	if (abs(trackHit.ref[0] + 999.) < 1.E-4 || abs(trackHit.ref[1] + 999.) < 1.E-4) {
	  validTrack = false;
	  break;
	}            
	// Ignore Detector-Type 3
	if (trackHit.detType == 3) {
	  break;
	}

	if (setJacobiMatrix (trackHit.jacobi, track, hit) == -1) {
	  PRINT_WARN ("Could not set JacobiMatrix for hit!");
	  break;
	}
#if USE_EIGEN_MATRIX_INVERSION
	matrixInverseEigen(trackHit.jacobi, trackHit.jacobiInverse, ORDER, ORDER);
#else
	matrixInverse(trackHit.jacobi, trackHit.jacobiInverse, ORDER, ORDER);
#endif
	newTrack.track.push_back (trackHit);
      }
      //             if(validTrack)
      //             	newEvent.push_back (newTrack);
      if (!validTrack)
	newTrack.track.clear ( );
      newEvent.push_back (newTrack);
      }
      mEvents.push_back (newEvent);
    }

    resetAttributeVectors ( );

    delete file;

    return 0;
  }

  // Return the number of stored events
  int
    EventReader::getEventCount (void)
    {
      return mEvents.size ( );
    }

  // Return the number of tracks associated with the event eventIdx
  int
    EventReader::getTrackCount (unsigned int eventIdx)
    {
      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
	return 0;
      }

      return mEvents.at (eventIdx).size ( );
    }

  // Return the data for the track trackIdx associated with the event eventIdx
  Track_t
    EventReader::getTrack (unsigned int eventIdx, unsigned int trackIdx)
    {
      Track_t track;

      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
	return track;
      }

      if (trackIdx >= mEvents.at (eventIdx).size ( )) {
	PRINT_WARN ("trackIdx out of bounds!");
	return track;
      }

      return mEvents.at (eventIdx).at (trackIdx);
    }

  // Return the complete event eventIdx
  KF_Event_t
    EventReader::getEvent (unsigned int eventIdx)
    {
      KF_Event_t event;

      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
	return event;
      }

      return mEvents.at (eventIdx);
    }

  std::vector<KF_Event_t> EventReader::getEvents ( )
  {
    return mEvents;
  }

  void EventReader::print(unsigned int eventIdx, unsigned int trackIdx)
  {
    size_t i, j;

    if (eventIdx >= mEvents.size()) {
      PRINT_WARN("eventIdx out of bounds!");
      return;
    } else if (trackIdx >= mEvents.at(eventIdx).size()) {
      PRINT_WARN("trackIdx out of bounds!");
      return;
    }
    Track_t& track = mEvents.at(eventIdx).at(trackIdx);


    std::cout << "event " << eventIdx << " track " << trackIdx<< std::endl;
    std::cout << "\tNormalLocX\tRefLocX\t\tNormalLocY\t\tRefLocY\t\tRefQoverP\t\tdetType\t\tmeasType\t\tsigmaDeltaP\t\tsigmaDeltaTheta"
      << std::endl;

    TrackData_t& trackData = track.track;

    for (i = 0; i < trackData.size(); i++) {
      /*myfile  << trackData.at (i).measurement[0] << ";"
	<< trackData.at (i).ref[0] << ";" << XFit[i] << ";"
	<< trackData.at (i).measurement[1] << ";"
	<< trackData.at (i).ref[1] << ";" << YFit[i] << ";"
	<< trackData.at (i).detType << std::endl;
	*/
      printf("%lu.\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%4.3e\t%d\t%d\t%4.3e\t%4.3e\n", i,
	  trackData.at(i).normal[0], trackData.at(i).ref[0],
	  trackData.at(i).normal[1], trackData.at(i).ref[1], 
	  trackData.at(i).ref[4], 
	  trackData.at(i).detType, trackData.at(i).is2Dim,
#if HAS_MATERIAL_EFFECTS
	  trackData.at(i).sigmaDeltaQoverPSquared, trackData.at(i).sigmaDeltaThetaSquared);
#else 
      0., 0.);
#endif
    }

    printf("\n\n");
  }

  void
    EventReader::printLoc (const char*filepath, unsigned int eventIdx,
	unsigned int trackIdx, trackHitData *data)
    {
      //ofstream myfile;
      Track_t track;
      size_t i,j;

      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
      } else if (trackIdx >= mEvents.at (eventIdx).size ( )) {
	PRINT_WARN ("trackIdx out of bounds!");
      } else {
	track = mEvents.at (eventIdx).at (trackIdx);
      }

      /* myfile.open (filepath);

	 myfile  << "NormalLocX;RefLocX;XFit;NormalLocY;RefLocY;YFit;detType"
	 << std::endl;
	 */




      TrackData_t trackData;
      trackData = track.track;

      /*if( (data[trackData.size()*2].data[0] > 5) || (data[trackData.size()*2].data[0] < -5) ||
	(data[trackData.size()*2].data[1] > 5) || (data[trackData.size()*2].data[1] < -5))
	{*/
      printf("Event %d Track %d\n",eventIdx, trackIdx);
      std::cout << "\tNormalLocX\tRefLocX\t\tXFit\t\tNormalLocY\tRefLocY\t\tYFit\t\tdetType"
	<< std::endl;
      for (i = 0; i < trackData.size ( ); i++) {
	/*myfile  << trackData.at (i).normal[0] << ";"
	  << trackData.at (i).ref[0] << ";" << data[i].data[0] << ";"
	  << trackData.at (i).normal[1] << ";"
	  << trackData.at (i).ref[1] << ";" << data[i].data[0]<< ";"
	  << trackData.at (i).detType << std::endl;
	  */
	printf ("%lu.\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%d\n", i,
	    trackData.at (i).normal[0], trackData.at (i).ref[0],
	    data[i].data[0], trackData.at (i).normal[1],
	    trackData.at (i).ref[1], data[i].data[1], trackData.at (i).detType);
      }

      printf ("\n\n");


      for (i = 0, j=trackData.size(); i < trackData.size(); j++,i++) {
	/*myfile  << trackData.at (i).normal[0] << ";"
	  << trackData.at (i).ref[0] << ";" << data[j].data[0] << ";"
	  << trackData.at (i).normal[1] << ";"
	  << trackData.at (i).ref[1] << ";" << data[j].data[1] << ";"
	  << trackData.at (i).detType << std::endl;
	  */
	printf ("%lu.\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%d\n", i,
	    trackData.at (i).normal[0], trackData.at (i).ref[0],
	    data[j].data[0], trackData.at (i).normal[1],
	    trackData.at (i).ref[1], data[j].data[1], trackData.at (i).detType);
      }


      printf ("\n\n");


      for (i = 0, j=trackData.size()*2; i < trackData.size(); j++,i++) {
	/*myfile  << trackData.at (i).normal[0] << ";"
	  << trackData.at (i).ref[0] << ";" << data[j].data[0] << ";"
	  << trackData.at (i).normal[1] << ";"
	  << trackData.at (i).ref[1] << ";" << data[j].data[1] << ";"
	  << trackData.at (i).detType << std::endl;
	  */
	printf ("%lu.\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%d\n", i,
	    trackData.at (i).normal[0], trackData.at (i).ref[0],
	    data[j].data[0], trackData.at (i).normal[1],
	    trackData.at (i).ref[1], data[j].data[1], trackData.at (i).detType);
      }
      std::cout << std::endl;
      //}
      //myfile.close ( );
    }

  void
    EventReader::printC (unsigned int eventIdx,
	unsigned int trackIdx, scalar_t *C)
    {
      Track_t track;
      int i, j, k;
      int hitCount;


      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
      } else if (trackIdx >= mEvents.at (eventIdx).size ( )) {
	PRINT_WARN ("trackIdx out of bounds!");
      } else {
	track = mEvents.at (eventIdx).at (trackIdx);
      }
      TrackData_t trackData;
      trackData = track.track;
      hitCount = trackData.size();

      for (i = 0; i< hitCount; i++) {
	printf("%d\n",i);
	for(j=0; j<ORDER; j++){
	  for(k=0; k<ORDER; k++){
	    printf("%.10f\t",C[(j*ORDER+k)+(i*ORDER*ORDER)] );
	  }
	  printf("\n");
	}
	printf("\n");
      }

      printf("\n\n");

      for (i = 0; i< hitCount; i++) {
	printf("%d\n",i);
	for(j=0; j<ORDER; j++){
	  for(k=0; k<ORDER; k++){
	    printf("%.10f\t",C[(j*ORDER+k)+(i*ORDER*ORDER)+(ORDER*ORDER*hitCount)] );
	  }
	  printf("\n");
	}
	printf("\n");
      }

      printf("\n\n");
    }

  void
    EventReader::printData4Mathematica (unsigned int eventIdx, unsigned int trackIdx)
    {
      Track_t track;
      scalar_t Vxx, Vyy, Vxy, X, Y, Phi, Theta, QoverP;

      if (eventIdx >= mEvents.size ( )) {
	PRINT_WARN ("eventIdx out of bounds!");
      } else if (trackIdx >= mEvents.at (eventIdx).size ( )) {
	PRINT_WARN ("trackIdx out of bounds!");
      } else {
	track = mEvents.at (eventIdx).at (trackIdx);
      }

      TrackData_t trackData;
      for (size_t i = 0; i < trackData.size ( ); i++) {
	printf ("F%lu =\n", i);
	printf ("{");
	for (int j = 0; j < ORDER; j++) {
	  printf ("{");
	  for (int k = 0; k < ORDER; k++) {
	    if (k != 4)
	      printf ("%f,", trackData.at (i).jacobi[j * ORDER + k]);
	    else
	      printf ("%f", trackData.at (i).jacobi[j * ORDER + k]);
	  }
	  if (j != 4)
	    printf ("},");
	  else
	    printf ("}");
	}
	printf ("};\n");


	Vxx = trackData.at (i).err_locX * trackData.at (i).err_locX;
	Vyy = trackData.at (i).err_locY * trackData.at (i).err_locY;
	Vxy = trackData.at (i).cov_locXY;

	printf ("V%lu =\n", i);
	printf ("{{%f,%f,0,0,0},{%f,%f,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};\n",
	    Vxx, Vxy, Vxy, Vyy);


	X       = trackData.at (i).normal[0] - trackData.at (i).ref[0];
	Y       = trackData.at (i).normal[1] - trackData.at (i).ref[1];
	Phi     = trackData.at (i).normal[2] - trackData.at (i).ref[2];
	Theta   = trackData.at (i).normal[3] - trackData.at (i).ref[3];
	QoverP  = trackData.at (i).normal[4] - trackData.at (i).ref[4];

	printf ("m%lu =\n", i);
	printf ("{{%f}, {%f}, {%f}, {%f}, {%f}};\n", X, Y, Phi, Theta, QoverP);

	printf ("\n");
      }
    }
