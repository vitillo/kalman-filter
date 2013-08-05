#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>

#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"

#include "TInterpreter.h"
#include "TSystem.h"

#include "ROOT_Ttree_io/event_reader.h"
#include "ROOT_Ttree_io/event_writer.h"
#include "serial/KalmanFilterSerial.h"
#include "ispc/KalmanFilter.h"
#include "clock.h"

const char * const kTreeName = "InDetTrackTree";

static EventReader eventReader;
static EventWriter eventWriter;

using namespace std;

void runISPC(){
  int eventCount = eventReader.getEventCount();
  vector<vector<Track_t> > tracksPerEvent;
  vector<int> maxHitsPerEvent;

  for(int i = 0; i < eventCount; i++){
    KF_Event_t event = eventReader.getEvent(i);
    vector<Track_t> tracks;
    int trackCounts = event.size();
    int maxHits = 0;

    for(int j = 0; j < trackCounts; j++){
      Track_t track = eventReader.getTrack(i, j);

      if(track.track.empty())
	continue;

      maxHits = track.track.size() > maxHits ? track.track.size() : maxHits;
      tracks.push_back(track);
    }

    sort(tracks.begin(), tracks.end(), [](const Track_t &a, const Track_t &b){return a.track.size() > b.track.size();});
    tracksPerEvent.push_back(tracks);
    maxHitsPerEvent.push_back(maxHits);
  }

  ispc::KalmanFilter filter;
  ispc::KalmanFilter_initialize(&filter, 200);

  ispc::KalmanFilterParameter param;
  KalmanFilterParameter_initialize(&param, 200, 20);

  for(int i = 0; i < tracksPerEvent.size(); i++){
    vector<Track_t> &tracks = tracksPerEvent[i];

    ispc::KalmanFilter_reinitialize(&filter, tracks.size());
    KalmanFilterParameter_reinitialize(&param, tracks, maxHitsPerEvent[i]);
    ispc::startFilter(&filter, &param);

  }

  ispc::KalmanFilter_deallocate(&filter);
  ispc::KalmanFilterParameter_deallocate(&param);
}

void runSerial(){
  int eventCount = eventReader.getEventCount();

  for(int i = 0; i < eventCount; i++){
    KF_Event_t event = eventReader.getEvent(i);
    int trackCount = event.size();

    KalmanFilterSerial kalmanFilter;

    for(int j = 0; j < trackCount; j++){
      Track_t track = eventReader.getTrack(i, j);
      kalmanFilter.setTrackData(track);
      kalmanFilter.startFilter();
    }
  }
}

void benchmark(){
  float tp, t1;
  start_clock();
  runSerial();
  t1 = stop_clock();

  start_clock();
  runISPC();
  tp = stop_clock();

  cout << "ISPC speedup over serial: " << t1/tp << endl;
}

void benchmarkPreload(){
  int eventCount = eventReader.getEventCount();
  vector<vector<Track_t> > tracksPerEvent;
  vector<ispc::KalmanFilterParameter> eventParams;

  for(int i = 0; i < eventCount; i++){
    KF_Event_t event = eventReader.getEvent(i);
    int trackCounts = event.size();
    vector<Track_t> tracks;

    for(int j = 0; j < trackCounts; j++){
      Track_t track = eventReader.getTrack(i, j);

      if(track.track.empty())
	continue;

      tracks.push_back(track);
    }

    sort(tracks.begin(), tracks.end(), [](const Track_t &a, const Track_t &b){return a.track.size() > b.track.size();});

    ispc::KalmanFilterParameter param;
    KalmanFilterParameter_initialize(&param, tracks, 20);

    eventParams.push_back(param);
    tracksPerEvent.push_back(tracks);
  }

  double t1, tp;

  start_clock();
  for(int i = 0; i < tracksPerEvent.size(); i++){
    vector<Track_t> &tracks = tracksPerEvent[i];
    KalmanFilterSerial kalmanFilter;

    for(int j = 0; j < tracks.size(); j++){
      kalmanFilter.setTrackData(tracks[j]);
      kalmanFilter.startFilter();
    }
  }
  t1 = stop_clock();

  ispc::KalmanFilter filter;
  KalmanFilter_initialize(&filter, 200);

  start_clock();
  for(int i = 0; i < eventParams.size(); i++){
    vector<Track_t> &tracks = tracksPerEvent[i];

    ispc::KalmanFilter_reinitialize(&filter, tracks.size());
    ispc::startFilter(&filter, &eventParams[i]);
  }
  tp = stop_clock();
  ispc::KalmanFilter_deallocate(&filter);

  cout << "ISPC speedup over serial (with preload): " << t1/tp << endl;

  for(int i = 0; i < eventParams.size(); i++)
    ispc::KalmanFilterParameter_deallocate(&eventParams[i]);
}

int main (int argc, char *argv[]){
  KF_Event_t event;
  Track_t track;
  Int_t ret, eventCount;

  char* kOutFileLocation = new char[100];
  strncpy (kOutFileLocation, "dump.root",40);

  char* kFileLocation = new char[100];
  strncpy (kFileLocation, "testTrackD3PD-noExtraMatLayer.root", 50);

  char* kOutTreeName = new char[100];
  strncpy (kOutTreeName, "FitResults",20);

  char* kPrefix = new char[100]; strncpy (kPrefix, "pseudoTrk_", 20);

  for (int i = 1; i < argc; i++) {
    if (!strncmp(argv[i], "-i=", 3)) {
      kFileLocation = &(argv[i][3]);
    }
    if (!strncmp(argv[i], "-o=", 3)) {
      kOutFileLocation =  &(argv[i][3]);
    }
    if (!strncmp(argv[i], "-t=", 3)) {
      kOutTreeName = &(argv[i][3]);
    }
    if (!strncmp(argv[i], "-p=", 3)) {
      kPrefix = &(argv[i][3]);
    }
  }

  cout << "Input file: " << kFileLocation << endl;
  cout << "Prefix: " << kPrefix << endl;
  cout << "Output file: " << kOutFileLocation << endl;
  cout << "Tree name: " << kOutTreeName << endl;

  if ((ret = eventReader.parse (kFileLocation, kTreeName, kPrefix)) == -1) {
    cerr << "Parsing file " << kFileLocation << " failed" << endl;
    return (-1);
  }

  if (eventWriter.initTree(kOutFileLocation, kOutTreeName, kPrefix)) {
    cerr << "Creation of output file " << kOutFileLocation << " failed" << endl;
    return (-1);   
  }

  runSerial(); //warmup

  benchmarkPreload();
  benchmark();

  eventWriter.closeFile();
  return 0;
}
