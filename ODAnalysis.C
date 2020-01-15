// C++ includes
#include  <iostream>
#include  <stdlib.h>

//ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "math.h"
#include "TSpectrum2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH2.h"
#include "TF2.h"
#include "TMath.h"
#include "TROOT.h"


//WCSim Includes
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootLinkDef.hh"

//Definitions
#define PI 3.141592654
#define RadiusID 3.55e+03
#define HeightID 5.50e+03
#define RadiusOD 3.60e+03
#define HeightOD 5.60e+03
/*
*How to run:
*
* enter in the terminal root -l llib.C 'ODAnalysis.C("WCSim.root","outputFile.root",false)' to run the code
* where you replace WCSim.root with your file name and outputfile with the name you wish to save it under
*/

// a structure to hold a pmt id along with its charge
typedef struct {
  double charge;
  double id;
} pmt;

typedef struct {
  double charge;
  int xbin;
  int ybin;
  double x;
  double y;
  double z;
  double totCharge;
} section;

// a function to be used to sort pmt charges
bool CompPMT(pmt a, pmt b)
{
  return a.charge > b.charge;
}
// a function to be used to sort pmt charges
bool CompSection(section a, section b)
{
  return a.charge > b.charge;
}
// Degrees to Radians conversions
double RadToDeg(double x){
  return x*180/PI;
}

//Radians to Degrees conversions
double DegToRad(double x){
  return x*PI/180;
}
void PrintPMT(int pmt_ID1, WCSimRootGeom *geo){

  if (pmt_ID1 <0) exit(1);
    WCSimRootPMT pmt1 = geo->GetPMT(pmt_ID1);

    double pmt1_x = pmt1.GetPosition(0);
    double pmt1_y = pmt1.GetPosition(1);
    double pmt1_z = pmt1.GetPosition(2);

    std::cout << "PMT Position:" << std::endl;
    std::cout << "x: " << pmt1_x << std::endl;
    std::cout << "y: " << pmt1_y << std::endl;
    std::cout << "z: " << pmt1_z << std::endl;

}
double GetDistance(int pmt_ID1, int pmt_ID2, WCSimRootGeom *geo){

  double distance;

  if ( pmt_ID2 == -1){
    distance = -1;
  }

  else{
    WCSimRootPMT pmt1 = geo->GetPMT(pmt_ID1);
    WCSimRootPMT pmt2 = geo->GetPMT(pmt_ID2);

    double pmt1_x = pmt1.GetPosition(0);
    double pmt1_y = pmt1.GetPosition(1);
    double pmt1_z = pmt1.GetPosition(2);

    double pmt2_x = pmt2.GetPosition(0);
    double pmt2_y = pmt2.GetPosition(1);
    double pmt2_z = pmt2.GetPosition(2);

    distance = sqrt( (pow((pmt2_x - pmt1_x), 2)) + (pow((pmt2_y - pmt1_y), 2)) + (pow((pmt2_z - pmt1_z), 2)) );
    distance = distance / 100. ;
  }

  return distance;
}
// finds angles using dot product
double Angles(double a[2], double b[2]){

  double angle = 0;
  double dot = a[0]*b[0] + a[1]*b[1];
  double aval = sqrt( pow(a[0],2) + pow(a[1],2) );
  double bval = sqrt( pow(b[0],2) + pow(b[1],2) );

  angle = acos( dot/(aval*bval) );

  return angle;

}

void Reorder(double charges[3], double ids[3]){
  std::vector<pmt> vec;
  pmt a;
  for (int i = 0; i < 3; i++){
    a.charge = charges[i];
    a.id = ids[i];
    vec.push_back(a);
  }

  std::sort(vec.begin(), vec.end(), CompPMT); // sorts highest to lowest

  for (int i = 0; i < 3; i++){
    charges[i] = vec[i].charge;
    ids[i] = vec[i].id;
  }


}

// Clustering algorithms
void CylinderToSquare(double square[2], int tubeID, double charge, WCSimRootGeom *geo, double radius, double height){

  // Find out where the PMT is in the tank
  double tube[3];
  int cylLoc = geo->GetPMT(tubeID).GetCylLoc();
  tube[0] = geo->GetPMT(tubeID).GetPosition(0);
  tube[1] = geo->GetPMT(tubeID).GetPosition(1);
  tube[2] = geo->GetPMT(tubeID).GetPosition(2);

  //Top (OD || ID)
  if ( cylLoc == 5 || cylLoc == 0){
          square[0] = tube[0];
          square[1] = tube[1] + radius + height/2. + 1.;
  }
  //Bot (OD || ID)
  else if ( cylLoc == 3 || cylLoc == 2){
          square[0] = tube[0];
          square[1] = -(height/2. +radius +tube[1] + 1. );
  }
  //Barrel OD
  else {

          double zero[2] = {0,-1};
          double values[2] = {tube[0],tube[1]};
          double angle = Angles(zero, values);

          double length = angle*radius ;
          if (tube[0]<0) length *= -1;
          square[0] = length;
          square[1] = tube[2];
;
  }

}
void SquareToCylinder(double square[2], double cylinder[3], double radius, double height){

  // Find out where the PMT is in the tank

  //Top (OD || ID)
  if ( square[1] >= (height/2. + 1.) ){
          cylinder[0] = square[0];
          cylinder[1] = square[1] - (radius + height/2. + 1.);
          cylinder[2] = height/2.;
  }
  //Bot (OD || ID)
  else if ( square[1] <= - (height/2. + 1.) ){
          cylinder[0] = square[0];
          cylinder[1] = -(square[1] + height/2. +radius + 1.);
          cylinder[2] = -height/2.;
  }
  //Barrel OD
  else {

          double angle = square[0]/radius;
          cylinder[0] = radius*cos( angle - PI/2.);
          cylinder[1] = radius*sin( abs(angle) - PI/2.);
          cylinder[2] = square[1];

  }

}

int FindBin(int numBins, double low, double high, double val){
    double maxlength = high - low;
    double step = maxlength/ (double)numBins;

    double binNum = (val - low)/step;
    return (int) binNum;
}

void GetNeighbours(int a[8][2], int i, int j){

  a[0][0] = i+1; a[0][1] =  j+1;
  a[1][0] = i+1; a[1][1] =  j-1;
  a[2][0] = i-1; a[2][1] =  j+1;
  a[3][0] = i-1; a[3][1] =  j-1;

  a[4][0] = i; a[4][1] = j+1;
  a[5][0] = i; a[5][1] = j-1;
  a[6][0] = i+1; a[6][1] = j;
  a[7][0] = i-1; a[7][1] = j;

}

std::vector<section> FindPeaks(double array[40][20]){

std::vector<section> peaks;

  for (int i = 0; i < 40; i++){
    for (int j = 0; j < 20; j++){
      int larger = 0;
      double totcharge = array[i][j];
      int neighbours[8][2];
      GetNeighbours( neighbours, i, j );
      for (int n = 0; n < 8; n++){

        if ( neighbours[n][0] < 0
            || neighbours[n][1] < 0
            || neighbours[n][0] > 39
            || neighbours[n][1] > 19
          ) continue;
          double val = array[ neighbours[n][0] ][ neighbours[n][1] ];
        if (val < array[i][j]) {
          totcharge += val;
          larger++;
        }

      } // end of loop over n
      if (larger == 8){
        section a;
        a.xbin = i;
        a.ybin = j;
        a.charge = array[i][j];
        a.totCharge = totcharge;
        peaks.push_back(a);
      }

    } // end loop over y (j index)
  }// end loop over x (i index)

  // sort out the clusters from highest to lowest charge
  std::sort(peaks.begin(), peaks.end(), CompSection);

  // find out the distance between clusters
  int i = 0;
  int j = (peaks.size()-1);
  std::cout << "PeakSize: " <<  peaks.size() << std::endl;
  if(peaks.size() > 0){

    while (i < (peaks.size()-1) ){
      //std::cout << "Mr hello: " << i << " "<< j <<  std::endl;
      double Val1[2]; double Val2[2];

      Val1[0] = -20350 + ( ((double)peaks[i].xbin + 0.5))*2*20350/40;
      Val2[0] = -20350 + ( ((double)peaks[i+1].xbin + 0.5))*2*20350/40;
      Val1[1] = -10000 + ( ((double)peaks[i].ybin + 0.5))*2*20350/40;
      Val2[1] = -10000 + ( ((double)peaks[i+1].ybin + 0.5))*2*20350/40;
      double cluster1[3];
      double cluster2[3];
      SquareToCylinder(Val1, cluster1, RadiusOD, HeightOD);
      SquareToCylinder(Val2, cluster2, RadiusOD, HeightOD);

      double difflength = sqrt( pow( (cluster2[0] - cluster1[0]  ) ,2 )  + pow( (cluster2[1] - cluster1[1]  ) ,2 ) + pow( (cluster2[2] - cluster1[2]  ) ,2 )  );

      peaks[i].x = cluster1[0];
      peaks[i].y = cluster1[1];
      peaks[i].z = cluster1[2];

      if (difflength < 1400 ) { // less than the distance of one box corner to corner
        // merge clusters
        peaks[i].charge += peaks[i+1].charge;
        peaks[i].totCharge += peaks[i+1].totCharge; // this double counts some charges but we can remove the double counted charge later
        peaks[i].x = (cluster1[0] + cluster2[0])/2;
        peaks[i].y = (cluster1[1] + cluster2[1])/2;
        peaks[i].z = (cluster1[2] + cluster2[2])/2;

        peaks.erase( peaks.begin() +i + 1 ); // removes the merged element

      }else if(peaks[i+1].totCharge/peaks[0].totCharge < 0.2 ){
        peaks.erase( peaks.begin() +i + 1 ); // removes the small element
      }
      else{i++;}

    }

  }




  return peaks;

}


void loadlibs(){

  char *wcsimdirenv;
  wcsimdirenv = getenv("WCSIMDIR");

  if (wcsimdirenv != NULL) {
    gSystem->Load("$WCSIMDIR/libWCSimRoot.so");
    gSystem->Load("$WCSIMDIR/libWCSimRoot.rootmap");
    gSystem->Load("$WCSIMDIR/libWCSimRootDict_rdict.pcm");
  } else {
    std::cout << "ERROR: WCSIMDIR environment variable has not been set." << std::endl;
  }


}  //end of loadlibs function

void ODAnalysis( const char *inFileName = "wcsim.root", const char *outFileName = "ODHits.root", bool verbosity = 0){
  // loadlibs(); // loads the required libraries

  bool createCanvases = false;

  // Some nice formatting for text options
  std::cout << std::scientific; //all numbers appear in scientific notation
  std::cout << std::setprecision(2); //sets all numbers to output to no more than 2 D.P.
  std::cout << std::left; // Sets the text justification to the left
  const int txtW = 20; // Width of box 'holding' text
  const int numW = 10; // Width of box 'holding' numbers


  // Open the WCSim file

    TFile *inFile = new TFile(inFileName, "READ");
    if ( !inFile->IsOpen() ){
      std::cout << "Error: Could not open file \"" << inFileName << "\"." <<std::endl;

    } else if (verbosity) {
      std::cout << "Input file: " << inFileName << std::endl;

    }

	// Get a pointer to the tree from the input file
	TTree *wcsimTree = (TTree*) inFile->Get("wcsimT");

	// Get the number of events in the tree
	long int nEvent = wcsimTree->GetEntries();
	if (verbosity) { std::cout << "Number of events: "<< nEvent << std::endl;}

	// Create a WCSimRootEvent to put stuff from the tree in
	WCSimRootEvent *wcsimRoot = new WCSimRootEvent();
	WCSimRootEvent *wcsimRootID = new WCSimRootEvent();

	//ID Event readings

	TBranch *IDbranch = wcsimTree->GetBranch("wcsimrootevent");
	IDbranch->SetAddress(&wcsimRootID);

	//Force Deletion to prevent memory leak
	wcsimTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	//OD Event readings and branch creation

	// Set the branch address for reading from the tree
	TBranch *branch = wcsimTree->GetBranch("wcsimrootevent_OD");
	branch->SetAddress(&wcsimRoot);

	// Force deletion to prevent memory leak
	wcsimTree->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);

	// Load the geometry tree (only 1 "event")
	TTree* geoTree = (TTree*) inFile->Get("wcsimGeoT");
	WCSimRootGeom *geo = 0;
	geoTree->SetBranchAddress("wcsimrootgeom", &geo);
	if (verbosity) {std::cout << "Geotree has " << geoTree->GetEntries() << " entries." << std::endl;}
	geoTree->GetEntry(0);

	// Start with the main trigger as it always exists and contains most of the info
	WCSimRootTrigger *wcsimTriggerID;
	WCSimRootTrigger *wcsimTriggerOD;

	// Create an output file
	TFile *outFile = new TFile(outFileName, "RECREATE");
	TTree *outBranch = new TTree("simulation", "simulation");

  TH1D *clusterHist = new TH1D("clusterHist", "clusterHist", 21, 0, 20);

	// Detector Geometry Details

	int MAXPMT = geo->GetWCNumPMT(); //Get the maximum number of PMTs in the ID
	int MAXPMTA = geo->GetODWCNumPMT(); //Get the maximum number of PMTs in the OD


	// OD Event Analysis
  for (int ev = 0; ev < nEvent; ev++){ // Loop over events
  //for (int ev = 21; ev < 22; ev++){ // Loop over events
	  wcsimTree->GetEntry(ev);
	  wcsimTriggerOD = wcsimRoot->GetTrigger(0);
	  int numTriggers = wcsimRoot->GetNumberOfEvents();
	  int numSubTriggers = wcsimRoot->GetNumberOfSubEvents();

	  for (int nTrig = 0; nTrig < numTriggers; nTrig++){



	    wcsimTriggerOD = wcsimRoot->GetTrigger(nTrig);
	    int numTracks = wcsimTriggerOD->GetNtrack();

	    if ( numTracks != 0){
        WCSimRootTrack * trackOD = (WCSimRootTrack*) wcsimTriggerOD->GetTracks()->At(0);

        double tankArray[40][20] = {}; // split tank into 40 by 20

	      double vtxX = wcsimTriggerOD->GetVtx(0);
	      double vtxY = wcsimTriggerOD->GetVtx(1);
	      double vtxZ = wcsimTriggerOD->GetVtx(2);
	      double dirX = trackOD->GetDir(0);
	      double dirY = trackOD->GetDir(1);
	      double dirZ = trackOD->GetDir(2);

	      double energy = trackOD->GetE();
	      double rawHits = 0;
	      double digiHits = 0;
        double largestHit = -1;
        double currentDist = 9999999;
        int largestHitID = -1;
        int largestHitID2 = -1;

	      int numPMTsHit = wcsimTriggerOD->GetNcherenkovhits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
	      int numPMTsDigiHit = wcsimTriggerOD->GetNcherenkovdigihits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
        if (verbosity) {std::cout << "Number of OD Hits: " << numPMTsDigiHit << std::endl;}
        if (verbosity) {std::cout << "Number triggers: " << nTrig << std::endl;}

    		for (int i = 0; i < numPMTsDigiHit; i++){

    		  WCSimRootCherenkovDigiHit *cherenkovDigiHit = (WCSimRootCherenkovDigiHit*) wcsimTriggerOD->GetCherenkovDigiHits()->At(i);
    		  double tmpDigiHits = cherenkovDigiHit->GetQ();
          int tubeID = cherenkovDigiHit->GetTubeId();
          int ODtubeID = tubeID + MAXPMT - 1;

          double tube[3];
          int cylLoc = geo->GetPMT(ODtubeID).GetCylLoc();
          tube[0] = geo->GetPMT(ODtubeID).GetPosition(0);
          tube[1] = geo->GetPMT(ODtubeID).GetPosition(1);
          tube[2] = geo->GetPMT(ODtubeID).GetPosition(2);

          // testing
          //CylinderToSquare(square[2], tubeID, charge, geo, radius, height);
          //SquareToCylinder(square[2], Cylinder[3], radius, height);
          double square[2] = {};
          double cylinder[3] = {};
          CylinderToSquare(square, ODtubeID, tmpDigiHits, geo, RadiusOD, HeightOD);
          tankArray[FindBin(40, -20350, 20350, square[0])][FindBin(20, -10000, 10000, square[1])] += tmpDigiHits;
          //SquareToCylinder(square, cylinder, RadiusOD, HeightOD);

    		} // End of loop over digi hits


        std::vector<section> peaks = FindPeaks(tankArray);
        if (verbosity) {std::cout << "Number of Peaks found: " << peaks.size() << std::endl;}
        clusterHist->Fill(peaks.size());
        /*
        for (int chkpk = 0; chkpk < peaks.size(); chkpk++ ){
          std::cout << "Peak " << chkpk << " : " << peaks[chkpk].x << ", " << peaks[chkpk].y << ", " << peaks[chkpk].z << ", " << peaks[chkpk].totCharge << std::endl;
        }
        */

	    }

	  } // End of loop over triggers

	}


	if (createCanvases) {
    TCanvas *c1 = new TCanvas("c1", "c1");
    clusterHist->Draw();
	}


  clusterHist->Write();
	outFile->Write(); // Write all of objects to the output file.
	outFile->Close(); // Close the output file.

}
