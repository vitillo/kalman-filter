//==================================================================================================
//
// Name         : EventWriter.h
// Description  : Class to parse the test ROOT file.
//
//==================================================================================================

#ifndef H_GUARD_EVENTWRITER
#define H_GUARD_EVENTWRITER

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

const std::string kFitLocX = "allhit_fitLocX";
const std::string kFitLocY = "allhit_fitLocY";
const std::string kFitPhi = "allhit_fitPhi";
const std::string kFitTheta = "allhit_fitTheta";
const std::string kFitQoverP = "allhit_fitQoverP";
const std::string kFitPerigeeD0 = "fitD0";
const std::string kFitPerigeeZ0 = "fitZ0";
const std::string kFitPerigeePhi = "fitPhi";
const std::string kFitPerigeeTheta = "fitTheta";
const std::string kFitPerigeeQoverP = "fitQoverP";
const std::string kFitCov = "allhit_fitCov";
const std::string kFitPerigeeCov = "fitCov";



// typedef std::vector<TrackHit_t> TrackData_t;
// 
// struct TrackStruct {
// 	TrackData_t track;
// 	TrackInfo_t info;
// 	TrackInfo_t truthTrackInfo;
// };
// 
// typedef TrackStruct Track_t;
// typedef std::vector<Track_t> KF_Event_t;
// 
// typedef std::vector<std::vector<float> > JacobiVal_t;

class EventWriter {
public:
    EventWriter(void);
    ~EventWriter(void);

    int addTrackFitData(const KalmanFilterParameter_t& fittedParams);
    int fillEvent();
    int initTree(const char* path, const char* treeName, const char* prefix);
    int closeFile();
    
private:
    int setBranches(TTree *tree, const char* prefix);

    std::vector<std::vector<float> > * mFitLocX;
    std::vector<std::vector<float> > * mFitLocY;
    std::vector<std::vector<float> > * mFitPhi;
    std::vector<std::vector<float> > * mFitTheta;
    std::vector<std::vector<float> > * mFitQoverP;

    std::vector<float> * mFitPerigeeD0;
    std::vector<float> * mFitPerigeeZ0;
    std::vector<float> * mFitPerigeePhi;
    std::vector<float> * mFitPerigeeTheta;
    std::vector<float> * mFitPerigeeQoverP;

    //std::vector< std::vector<std::vector<float> >* >  mFitCovariance;
    std::vector<std::vector<float> >*  mFitCovariance[25];
    std::vector<float> *               mFitCovariancePerigee[25];
    
    TTree* mTree;
    TFile* mFile;
};

#endif
