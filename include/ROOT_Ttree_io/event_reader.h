//==================================================================================================
//
// Name         : EventReader.h
// Description  : Class to parse the test ROOT file.
//
//==================================================================================================

#ifndef H_GUARD_EVENTREADER
#define H_GUARD_EVENTREADER

#define PRINT_ERR(x) \
    std::cerr << "ERROR: " << __FILE__ << ", line: " << __LINE__ << " " << x << std::endl;

#define PRINT_WARN(x) \
    std::cerr << "WARNING: " << __FILE__ << ", line: " << __LINE__ << " " << x << std::endl;

    
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>

#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "TInterpreter.h"
#include "TSystem.h"

#include "serial/common.h"
#include "serial/common_gsl.h"

const std::string kLocX = "allhit_locX";
const std::string kLocY = "allhit_locY";
const std::string kErrLocX = "allhit_err_locX";
const std::string kErrLocY = "allhit_err_locY";
const std::string kCovLocXY = "allhit_cov_locXY";
const std::string kSigmaDeltaPhiScat = "allhit_sigmaDeltaPhiScat";
const std::string kSigmaDeltaThetaScat = "allhit_sigmaDeltaThetaScat";
const std::string kSigmaDeltaEloss = "allhit_sigmaDeltaEloss";

const std::string kDetType = "allhit_detType";
const std::string kBEC = "allhit_bec";
const std::string kRefLocX = "allhit_refLocX";
const std::string kRefLocY = "allhit_refLocY";
const std::string kRefPhi = "allhit_refPhi";
const std::string kRefTheta = "allhit_refTheta";
const std::string kRefQoverP = "allhit_refQoverP";
const std::string kQoverP = "qoverp";
const std::string kD0 = "d0";
const std::string kZ0 = "z0";
const std::string kPhi = "phi";
const std::string kTheta = "theta";
const std::string kMci = "mc_index";
const std::string kMcPerigeeOk = "mc_perigee_ok";
// not needed anymore... can directly use mc_perigee_d0, etc.
const std::string kTti = "mc_truthtrack_index";


const std::string kTransportJac00 = "allhit_transportJac00";
const std::string kTransportJac01 = "allhit_transportJac01";
const std::string kTransportJac02 = "allhit_transportJac02";
const std::string kTransportJac03 = "allhit_transportJac03";
const std::string kTransportJac04 = "allhit_transportJac04";
const std::string kTransportJac10 = "allhit_transportJac10";
const std::string kTransportJac11 = "allhit_transportJac11";
const std::string kTransportJac12 = "allhit_transportJac12";
const std::string kTransportJac13 = "allhit_transportJac13";
const std::string kTransportJac14 = "allhit_transportJac14";
const std::string kTransportJac20 = "allhit_transportJac20";
const std::string kTransportJac21 = "allhit_transportJac21";
const std::string kTransportJac22 = "allhit_transportJac22";
const std::string kTransportJac23 = "allhit_transportJac23";
const std::string kTransportJac24 = "allhit_transportJac24";
const std::string kTransportJac30 = "allhit_transportJac30";
const std::string kTransportJac31 = "allhit_transportJac31";
const std::string kTransportJac32 = "allhit_transportJac32";
const std::string kTransportJac33 = "allhit_transportJac33";
const std::string kTransportJac34 = "allhit_transportJac34";
const std::string kTransportJac40 = "allhit_transportJac40";
const std::string kTransportJac41 = "allhit_transportJac41";
const std::string kTransportJac42 = "allhit_transportJac42";
const std::string kTransportJac43 = "allhit_transportJac43";
const std::string kTransportJac44 = "allhit_transportJac44";

const std::string kCompetingRotIds = "allhit_id";
const std::string kRotLocX = "allhit_rotLocX";
const std::string kRotLocY = "allhit_rotLocY";
const std::string kRotIds = "allhit_CompetingRIO_ids";
const std::string kRotProbs = "allhit_CompetingRIO_assgnProbs";
const std::string kRotCov00 = "allhit_err_rotCov00";
const std::string kRotCov01 = "allhit_err_rotCov01";
const std::string kRotCov10 = "allhit_err_rotCov10";
const std::string kRotCov11 = "allhit_err_rotCov11";




typedef std::vector<TrackHit_t> TrackData_t;

struct TrackStruct {
	TrackData_t track;
	TrackInfo_t info;
	TrackInfo_t truthTrackInfo;
};

typedef TrackStruct Track_t;
typedef std::vector<Track_t> KF_Event_t;

typedef std::vector<std::vector<float> > JacobiVal_t;

class EventReader {
public:
	EventReader(void);
	~EventReader(void);

	int parse(const char* path, const char* treeName, const char* prefix);

	int getEventCount(void);

	int getTrackCount(unsigned int eventIdx);

	std::vector<KF_Event_t> getEvents();

	Track_t getTrack(unsigned int eventIdx, unsigned int trackIdx);
	KF_Event_t getEvent(unsigned int eventIdx);

	void printLoc(const char* filepath, unsigned int eventIdx,
			unsigned int trackIdx, trackHitData *data);
	void printC(unsigned int eventIdx, unsigned int trackIdx, scalar_t *C);
	void printData4Mathematica(unsigned int eventIdx, unsigned int trackIdx);
    void print(unsigned int eventIdx, unsigned int trackIdx);
private:
	void resetAttributeVectors(void);

	int eventSanityCheck(void);

	int loadBranches(TTree *tree, const char* prefix);

	int setJacobiMatrix(scalar_t (&matrix)[ORDER * ORDER],
			unsigned int trackIdx, unsigned int hitIdx);

	std::vector<KF_Event_t> mEvents;

	std::vector<std::vector<float> > * mLocX;
	std::vector<std::vector<float> > * mLocY;
	std::vector<std::vector<float> > * mErr_locX;
	std::vector<std::vector<float> > * mErr_locY;
	std::vector<std::vector<float> > * mCov_locXY;
	std::vector<std::vector<int> > * mDetType;
	std::vector<std::vector<int> > * mBec;
	std::vector<std::vector<float> > * mRefLocX;
	std::vector<std::vector<float> > * mRefLocY;
	std::vector<std::vector<float> > * mRefPhi;
	std::vector<std::vector<float> > * mRefTheta;
	std::vector<std::vector<float> > * mRefQoverP;
#if HAS_MATERIAL_EFFECTS
    std::vector<std::vector<float> > * mSigmaDeltaTheta;
    std::vector<std::vector<float> > * mSigmaDeltaPhi;
    std::vector<std::vector<float> > * mSigmaDeltaE;
#endif // #if HAS_MATERIAL_EFFECTS

	std::vector<float> * mD0;
	std::vector<float> * mZ0;
	std::vector<float> * mQoverP;
	std::vector<float> * mPhi;
	std::vector<float> * mTheta;

	std::vector<int> * mMCI;
	std::vector<int> * mTTTI;
    std::vector<int> * mMCPerigeeOk;
	std::vector<float> * mTTD0;
	std::vector<float> * mTTZ0;
	std::vector<float> * mTTQoverP;
	std::vector<float> * mTTPhi;
	std::vector<float> * mTTTheta;

    std::vector<JacobiVal_t*> * mJacobi;
#if HAS_COMPETING_ROTS
    std::vector<std::vector<unsigned int> > * mCompetingRotIds;
    std::vector< std::vector<std::vector<float> > >* mRotLocX;
    std::vector< std::vector<std::vector<float> > >* mRotLocY;
    std::vector< std::vector<std::vector<int> > >* mRotIds;
    std::vector< std::vector<std::vector<float> > >* mRotProbs;
    std::vector< std::vector<std::vector<float> > >* mRotCov00;
    std::vector< std::vector<std::vector<float> > >* mRotCov01;
    std::vector< std::vector<std::vector<float> > >* mRotCov10;
    std::vector< std::vector<std::vector<float> > >* mRotCov11;
#endif //#if HAS_COMPETING_ROTS
};

#endif
