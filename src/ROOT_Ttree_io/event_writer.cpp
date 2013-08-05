//==================================================================================================
//
// Name         : EventWriter.cpp
// Description  : Class to parse the test ROOT file.
//
//==================================================================================================

#include "ROOT_Ttree_io/event_writer.h"
#include "TString.h"

// Constructor
EventWriter::EventWriter(void) :
  // Initialize pointers with NULL
  mFitLocX(0),
  mFitLocY(0),
  mFitPhi(0),
  mFitTheta(0),
  mFitQoverP(0),
  mFitPerigeeD0(0),
  mFitPerigeeZ0(0),
  mFitPerigeePhi(0),
  mFitPerigeeTheta(0),
  mFitPerigeeQoverP(0),
  //mFitCovariance(25, nullptr),
  mTree(0),
  mFile(0)
{
  for (size_t i = 0; i < 25; i++) {
    mFitCovariance[i] = 0;
    mFitCovariancePerigee[i] = 0;
  }
}

// Destructor
EventWriter::~EventWriter(void)
{
  delete mTree;
  delete mFile;

  delete mFitLocX; mFitLocX=0;
  delete mFitLocY; mFitLocY=0;
  delete mFitPhi; mFitPhi=0;
  delete mFitTheta; mFitTheta=0;
  delete mFitQoverP; mFitQoverP=0;

  delete mFitPerigeeD0; mFitPerigeeD0=0;
  delete mFitPerigeeZ0; mFitPerigeeZ0=0;
  delete mFitPerigeePhi; mFitPerigeePhi=0;
  delete mFitPerigeeTheta; mFitPerigeeTheta=0;
  delete mFitPerigeeQoverP; mFitPerigeeQoverP=0;

  for (size_t i = 0; i < 25; i++) {
    delete mFitCovariance[i]; mFitCovariance[i] = 0;
    delete mFitCovariancePerigee[i]; mFitCovariancePerigee[i] = 0;
  }
}

int EventWriter::initTree(const char* path, const char* treeName, const char* prefix){
  if (mFile) {
    PRINT_ERR ("output file already set!");
    return -1;
  }

  // Open data file and create TTree
  mFile = new TFile (path, "UPDATE");

  if (!mFile || mFile->IsZombie()) {
    delete mFile; mFile = 0;
    std::cerr << "Could not open " << path << std::endl;
    return -1;
  }

  TTree* tree = new TTree (treeName, "fitter output");
  return setBranches (tree, prefix);

  //return 0;
}

int EventWriter::closeFile(){
  if (!mFile) return 0;
  if (!mTree) return 0;

  mFile->cd();
  mTree->Write();
  delete mTree; mTree=0;
  mFile->Close();    
  delete mFile; mFile=0;
  return 0;
}

// Setup branch addresses for each track attribute
int EventWriter::setBranches (TTree *tree, const char* prefix)
{

  if (mTree) {
    PRINT_ERR ("tree already set!");
    return -1;
  }
  mTree = tree;

  const std::string bNamePrefix = prefix;
  std::string bName;

  mFitLocX = new std::vector<std::vector<float> >;
  bName = bNamePrefix + kFitLocX;
  if (!tree->Branch(bName.c_str(), &mFitLocX)) {
    PRINT_ERR ("creating branch failed for fitLocX");
    return -1;
  }

  mFitLocY = new std::vector<std::vector<float> >;
  bName = bNamePrefix + kFitLocY;
  if (!tree->Branch(bName.c_str(), &mFitLocY)) {
    PRINT_ERR ("creating branch failed for fitLocY");
    return -1;
  }

  mFitPhi = new std::vector<std::vector<float> >;
  bName = bNamePrefix + kFitPhi;
  if (!tree->Branch(bName.c_str(), &mFitPhi)) {
    PRINT_ERR ("creating branch failed for fitPhi");
    return -1;
  }

  mFitTheta = new std::vector<std::vector<float> >;
  bName = bNamePrefix + kFitTheta;
  if (!tree->Branch(bName.c_str(), &mFitTheta)) {
    PRINT_ERR ("creating branch failed for fitTheta");
    return -1;
  }

  mFitQoverP = new std::vector<std::vector<float> >;
  bName = bNamePrefix + kFitQoverP;
  if (!tree->Branch(bName.c_str(), &mFitQoverP)) {
    PRINT_ERR ("creating branch failed for fitQoverP");
    return -1;
  }

  mFitPerigeeD0 = new std::vector<float>;
  bName = bNamePrefix + kFitPerigeeD0;
  if (!tree->Branch(bName.c_str(), &mFitPerigeeD0)) {
    PRINT_ERR ("creating branch failed for fitD0");
    return -1;
  }

  mFitPerigeeZ0 = new std::vector<float>;
  bName = bNamePrefix + kFitPerigeeZ0;
  if (!tree->Branch(bName.c_str(), &mFitPerigeeZ0)) {
    PRINT_ERR ("creating branch failed for fitZ0");
    return -1;
  }

  mFitPerigeePhi = new std::vector<float>;
  bName = bNamePrefix + kFitPerigeePhi;
  if (!tree->Branch(bName.c_str(), &mFitPerigeePhi)) {
    PRINT_ERR ("creating branch failed for fitPhi");
    return -1;
  }

  mFitPerigeeTheta = new std::vector<float>;
  bName = bNamePrefix + kFitPerigeeTheta;
  if (!tree->Branch(bName.c_str(), &mFitPerigeeTheta)) {
    PRINT_ERR ("creating branch failed for fitTheta");
    return -1;
  }

  mFitPerigeeQoverP = new std::vector<float>;
  bName = bNamePrefix + kFitPerigeeQoverP;
  if (!tree->Branch(bName.c_str(), &mFitPerigeeQoverP)) {
    PRINT_ERR ("creating branch failed for fitQoverP");
    return -1;
  }

  TString bName2(bNamePrefix + kFitCov);
  for (unsigned  i=0; i < 5; i++) {
    for (unsigned  j=0; j < 5; j++) {
      mFitCovariance[i+5*j] = new std::vector<std::vector<float> >;
      bName2 = bNamePrefix + kFitCov;
      bName2 += i;
      bName2 += j;
      if (!tree->Branch(bName2.Data(), &(mFitCovariance[i+5*j]))) {
	PRINT_ERR ("creating branch failed for " << bName2);
	return -1;
      }
      mFitCovariancePerigee[i+5*j] = new std::vector<float>;
      bName2 = bNamePrefix + kFitPerigeeCov;
      bName2 += i;
      bName2 += j;
      if (!tree->Branch(bName2.Data(), &(mFitCovariancePerigee[i+5*j]))) {
	PRINT_ERR ("creating branch failed for " << bName2);
	return -1;
      }            
    }
  }

  return 0;
}

int EventWriter::addTrackFitData(const KalmanFilterParameter_t& fittedParams){
  if (!mTree) {
    PRINT_ERR ("tree not initialised!");
    return -1;
  }

  std::vector<float> fitLocX;
  std::vector<float> fitLocY;
  std::vector<float> fitPhi;
  std::vector<float> fitTheta;
  std::vector<float> fitQoverP;
  std::vector< std::vector<float> > fitCov(25);

  if (!(fittedParams.fitsResult && fittedParams.hitCount && fittedParams.hits)) {
    mFitLocX->push_back(fitLocX);
    mFitLocY->push_back(fitLocY);
    mFitPhi->push_back(fitPhi);
    mFitTheta->push_back(fitTheta);
    mFitQoverP->push_back(fitQoverP);
    mFitPerigeeD0->push_back(-999.);
    mFitPerigeeZ0->push_back(-999.);
    mFitPerigeePhi->push_back(-999.);
    mFitPerigeeTheta->push_back(-999.);
    mFitPerigeeQoverP->push_back(-999.);
    for (unsigned int i=0; i<25; ++i) {
      mFitCovariance[i]->push_back(fitCov[i]);
      mFitCovariancePerigee[i]->push_back(-999.);
    }
    return -1;
  }
  if (*fittedParams.hitCount < 1) {
    mFitLocX->push_back(fitLocX);
    mFitLocY->push_back(fitLocY);
    mFitPhi->push_back(fitPhi);
    mFitTheta->push_back(fitTheta);
    mFitQoverP->push_back(fitQoverP);
    mFitPerigeeD0->push_back(-999.);
    mFitPerigeeZ0->push_back(-999.);
    mFitPerigeePhi->push_back(-999.);
    mFitPerigeeTheta->push_back(-999.);
    mFitPerigeeQoverP->push_back(-999.);
    for (unsigned int i=0; i<25; ++i) {
      mFitCovariance[i]->push_back(fitCov[i]);
      mFitCovariancePerigee[i]->push_back(-999.);
    }
    return 0;
  }
  size_t hitCount = *fittedParams.hitCount;
  fitLocX.reserve(hitCount);
  fitLocY.reserve(hitCount);
  fitPhi.reserve(hitCount);
  fitTheta.reserve(hitCount);
  fitQoverP.reserve(hitCount);
  for (size_t i = 0; i < hitCount; ++i) {
    fitLocX.push_back(fittedParams.fitsResult[i].data[0] + fittedParams.hits[i].ref[0]);
    fitLocY.push_back(fittedParams.fitsResult[i].data[1] + fittedParams.hits[i].ref[1]);
    fitPhi.push_back(fittedParams.fitsResult[i].data[2] + fittedParams.hits[i].ref[2]);
    fitTheta.push_back(fittedParams.fitsResult[i].data[3] + fittedParams.hits[i].ref[3]);
    fitQoverP.push_back(fittedParams.fitsResult[i].data[4] + fittedParams.hits[i].ref[4]);
    for (unsigned int c=0; c<5; ++c) {
      for (unsigned int k=0; k<5; ++k) {
	fitCov[c+5*k].push_back(fittedParams.C_Inverse[c*5+k+i*25]);
      }
    }

  }

  if (fittedParams.fittedPerigee && fittedParams.C_fittedPerigee) {
    mFitPerigeeD0->push_back(fittedParams.fittedPerigee->d0);
    mFitPerigeeZ0->push_back(fittedParams.fittedPerigee->z0);
    mFitPerigeePhi->push_back(fittedParams.fittedPerigee->phi);
    mFitPerigeeTheta->push_back(fittedParams.fittedPerigee->theta);
    mFitPerigeeQoverP->push_back(fittedParams.fittedPerigee->qoverp);        
    for (unsigned int c=0; c<5; ++c) {
      for (unsigned int k=0; k<5; ++k) {
	mFitCovariancePerigee[c+5*k]->push_back(fittedParams.C_fittedPerigee[c*5+k]);
      }
    }
  } else {
    mFitPerigeeD0->push_back(-999.);
    mFitPerigeeZ0->push_back(-999.);
    mFitPerigeePhi->push_back(-999.);
    mFitPerigeeTheta->push_back(-999.);
    mFitPerigeeQoverP->push_back(-999.);        
    for (unsigned int i=0; i<25; ++i) {
      mFitCovariancePerigee[i]->push_back(-999.);
    }
  }

  mFitLocX->push_back(fitLocX);
  mFitLocY->push_back(fitLocY);
  mFitPhi->push_back(fitPhi);
  mFitTheta->push_back(fitTheta);
  mFitQoverP->push_back(fitQoverP);
  for (unsigned int i=0; i<25; ++i) {
    mFitCovariance[i]->push_back(fitCov[i]);
  }

  return 0;
}

int EventWriter::fillEvent(){
  if (!mTree) {
    PRINT_ERR ("tree not initialised!");
    return -1;
  }

  // fill data into tree
  mTree->Fill();
  // clear vectors for next event
  mFitLocX->clear();
  mFitLocY->clear();
  mFitPhi->clear();
  mFitTheta->clear();
  mFitQoverP->clear();
  mFitPerigeeD0->clear();
  mFitPerigeeZ0->clear();
  mFitPerigeePhi->clear();
  mFitPerigeeTheta->clear();
  mFitPerigeeQoverP->clear();
  for (unsigned int i=0; i<25; ++i) {
    mFitCovariance[i]->clear();
    mFitCovariancePerigee[i]->clear();
  }
  return 0;
}

