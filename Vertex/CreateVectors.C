/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Created on 20 Jul 2020				       *
 * 							       *
 * (Modified from Mahdi's CreateVectors.C, which can be found  *
 * at: https://github.com/mahditaani/ODScripts)		       *
 * 							       *
 * This file generates nuance format vector files as the input *
 * to WCSim						       *
 * To run, do 	root -l CreateVectors.C			       *
 * 							       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//Include C++
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <string>

//Include ROOT
#include <TFile.h>
#include <TTree.h>


//WCIDDiameter = 65.8*m; // = 69m - 2*(60cm ID wall + 1m OD)
//WCIDHeight   = 67.8*m; // = 73m - 2*(60cm ID wall + 2m OD)
#define nParticles 1000

#define PI 3.141592654
#define SEED 12345

#define IDRadius 3290 // in cm
#define IDHalfHeight  3390 // in cm
#define ODRadius 3450 // in cm
#define ODHalfHeight 3650 //in cm
#define Deadspace 600 // in cm
#define LargeCircleRadius 5000 // in cm
#define LargeCircleOffset 200 // in cm
#define FiducialCut 200 // in cm

#define RockThickness 1000 // 10 m

// G4double ekin = energy - mass;
#define muMass 105.6583745// mu- mass (in MeV)
#define muEnergyCerenkovThreshold 159.7404473 //minimum energy to produce Cerenkov radiation (including mass)
#define muRadiativeLoss 2.3 // 2MeV/cm loss for muon (in MeV) 0.3 MeV per cm Cerenkov loss
//#define FixedKE 1e6 // 10**6 MeV kinetic energy
#define FixedKE 2e6 // 2000 GeV kinetic energy
//#define FixedKE 1e3 

// Enumerator for the various event types
enum EventType {fc, pc, upmus, upmut, cosmic, randomloc};

std::default_random_engine generator(SEED);

// This function generates a random double between a and b. (uniform distribution)
double Random(double a, double b)
{
  std::uniform_real_distribution<double> quick(a,b);
  return quick(generator);
}

// This function generates an energy/cm lower than that required to traverse the distance
double LessEnergy()
{
  return Random(0.01, muRadiativeLoss-0.1);
}

// This function generates an energy/cm higher than that required to traverse the distance
double MoreEnergy()
{
  return Random(muRadiativeLoss+0.1, muRadiativeLoss+3);
}

// function to generate a point inside a cylinder. height is the half height and pos is the output vertex
void GeneratePointInCylinder(double rad, double height, double Pos[3])
{
  double genR = Random(0., pow(height,2));
  double genTheta = Random(0., 2*PI);
  Pos[0] = sqrt(genR)*cos(genTheta);
  Pos[1] = sqrt(genR)*sin(genTheta);
  Pos[2] = Random(-height, height);

}

// function to generate a point on the surface of a cylinder. height is the half height, pos is the output vertex
void GeneratePointOnCylinder(double rad, double height, double Pos[3])
{
  // First work out the area of the cylinder caps and barrel
  double areaCap = PI*rad*rad;
  double areaBarrel = height*PI*4*rad;
  double totalArea = areaBarrel + 2*areaCap;
  areaCap /= totalArea;
  areaBarrel /= totalArea;
  
  // Using probablility to ensure a proportional number of points on the surface of the cylinder
  double prob = Random(0.,1.);
  double genR = 0;
  double genTheta = Random(0., 2*PI);
  
  if (prob >= 0 && prob < areaCap){ // generate on top cap
    genR = Random(0., pow(rad,2));
    Pos[2] = height;
  }
  else if (prob >= areaCap && prob < 2*areaCap){ // generate on bottom cap
    genR = Random(0., pow(rad,2));
    Pos[2] = -height;
  }
  else if (prob >= 2*areaCap){ // generate on barrel
    genR = rad*rad;
    Pos[2] = Random(-height, height);
  }
  
  Pos[0] = sqrt(genR)*cos(genTheta);
  Pos[1] = sqrt(genR)*sin(genTheta);
  
}

// Function to generate the name of the output file
std::string GetOutFileName(EventType ev)
{
  std::string evTypeName;
  
  switch (ev){ 
    case fc: evTypeName = "FC"; break;
    case pc: evTypeName = "PC"; break;
    case upmus: evTypeName = "UPMUS"; break;
    case upmut: evTypeName = "UPMUT"; break;
    case cosmic: evTypeName = "COSMIC"; break;
    case randomloc: evTypeName = "RANDOM"; break;
    default: std::cerr << "Unknown event type. Please select from: randomloc, fc, pc, upmus, upmut, cosmic." 
		       << std::endl; exit(1);
  }
  
  std::string out = "VectorFile" + evTypeName + std::to_string(int(FixedKE/1000)) + "GeV"
		    + std::to_string(int(nParticles)) + "Evnts.dat";
  return out;
}


// Function to get the length between two vectors and output a normalised 
// direction from the first point to the second
double GetL(double aPos[3], double bPos[3], double dir[3] )
{
  
  double l = sqrt( pow(bPos[0] - aPos[0], 2) + pow(bPos[1] - aPos[1], 2) + pow(bPos[2] - aPos[2], 2)  );
  dir[0] = (bPos[0] - aPos[0])/l;
  dir[1] = (bPos[1] - aPos[1])/l;
  dir[2] = (bPos[2] - aPos[2])/l;
  
  return l;
}
  
// Function to generate a random position for each event
void GenerateVertex(EventType ev, double aPos[3], double dir[3], double &energy )
{  
  // The idea is to select two point and draw a line between them in order to 
  // generate a specific kind of event. e.g. pick a point above and below the 
  // detector to generate a cosmic muon.
  
  
  double bPos[3]; // coordinates of second point
  double length; // distance between two points
  
  if (ev == randomloc){
    // this is to generate muon with 1TeV KE scattered randomly on the cave wall with
    // random directions pointing towards the ID region
     
    // generate point a on the surface of the cave walls (OD)
    GeneratePointOnCylinder(ODRadius+RockThickness, ODHalfHeight+RockThickness, aPos); 

    // generate point b inside the fiducial volume
    GeneratePointInCylinder(IDRadius-FiducialCut, IDHalfHeight-FiducialCut, bPos);

    // get distance and direction between the two points
    length = GetL(aPos, bPos, dir); // dir now filled

    energy = FixedKE + muMass;

  }
  
  else if (ev == fc || ev == pc){
    // generate point a inside the fiducial volume
    GeneratePointInCylinder(IDRadius-FiducialCut, IDHalfHeight-FiducialCut, aPos);
    
    // generate point b on the edge of the ID
    GeneratePointOnCylinder(IDRadius, IDHalfHeight, bPos);
    
    // get distance and direction between the two points
    length = GetL(aPos, bPos, dir); // dir now filled
    
    if (ev == fc){
      energy = muEnergyCerenkovThreshold + length*LessEnergy(); // not enough energy to leave
    }
    else if (ev == pc){
      energy = muEnergyCerenkovThreshold + length*MoreEnergy(); // more than enough energy to leave
    }
  
  }
  else if (ev == upmus || ev == upmut){
    // generate point on the lower half surface of the cave walls
    GeneratePointOnCylinder(ODRadius, ODHalfHeight, aPos);
    aPos[2] = -1*fabs(aPos[2]); // ensures it is on the lower surface
    
    // generate point b on the upper half surface of the ID 
    GeneratePointOnCylinder(IDRadius, IDHalfHeight, bPos);
    bPos[2] = fabs(bPos[2]); // ensures it is on the upper surface
    
    // get distance and direction between the two points
    length = GetL(aPos, bPos, dir); // dir now filled
    
    if (ev == upmus){
      energy = muEnergyCerenkovThreshold + length*LessEnergy(); // not enough energy to leave
    }
    else if (ev == upmut){
      energy = muEnergyCerenkovThreshold + length*MoreEnergy(); // more than enough energy to leave
    }
    
    }
  else if (ev == cosmic){
    // generate point a on the upper half surface of the cave walls
    GeneratePointOnCylinder(ODRadius, ODHalfHeight, aPos);
    aPos[2] = fabs(aPos[2]); // ensures it is on the upper surface
    
    // generate point b on the lower half surface of the ID
    GeneratePointOnCylinder(IDRadius, IDHalfHeight, bPos);
    bPos[2] = -1*fabs(bPos[2]); // ensures it is on the lower surface
    
    // get distance and direction between the two points
    length = GetL(aPos, bPos, dir); // dir now filled
    
    energy = muEnergyCerenkovThreshold + length*MoreEnergy();
  
  }

}


//Function to create a nuance format vector file from the root input
void CreateVector(EventType ev = fc)
{

  //PDG code of the particle you want to generate
  //int particlePDG = 11; // e-
  //int particlePDG = -11; // e+
  int particlePDG = 13; // mu-
  //int particlePDG = -13; // mu+
  //int particlePDG = 22; // gamma
  
  std::string outputFileName = GetOutFileName(ev);
   
  //Variables used to store the infomation for the output
  double energyVal;
  double dirVal[3];
  double vtxVal[3];
  double timeVal;
  
  
  //Create the output file
  ofstream outputFile;
  outputFile.open(outputFileName.c_str());
  std::cout << "Output file " << outputFileName << " has been created." << std::endl;
  
  //Check to make sure the file is open
  if (!outputFile)
  {
    std::cerr << "Error: " << outputFileName << " could not be opened." << std::endl;
    exit(1);
  }
  
  //For each entry in the file, create a nuance formatted entry to the output
  for (int i = 0; i< nParticles; i++)
  {
  
    GenerateVertex(ev, vtxVal, dirVal, energyVal );
    
    //Dummy time value (until I can figure out where the time info is)
    timeVal = 0;
    
    //Some text that must be here for WCSim to read the file correctly
    outputFile << "$ begin" << std::endl;
    outputFile << "$ nuance mode 1" << std::endl;
    //The creation point and time of the particle being written
    outputFile << "$ vertex " << vtxVal[0] << " " << vtxVal[1] << " " << vtxVal[2] <<" " << timeVal << std::endl;
    //Dummy values that must be written for WCSim to read the file correctly
    outputFile << "$ track 14 " << energyVal << " 0.00 0.00 0.00 -1" << std::endl;
    outputFile << "$ track 2212 " << energyVal << " 0.00 0.00 0.00 -1" << std::endl;
    outputFile << "$ info 2 949000 0.000E+00" << std::endl;
    //The particle being created its PDG code, energy (MeV), direction followed by a 0 (For WCSim to create this particle)
    outputFile << "$ track " << particlePDG << " " << energyVal <<" " << dirVal[0] << " " << dirVal[1] << " " << dirVal[2] << " 0" << std::endl;
    //Ends the creation of the particle
    outputFile << "$ end" << std::endl;

  }
  
  std::cout << "Number of particles in data file: " << nParticles << std::endl;
  //Close the output file
  outputFile.close();
  std::cout << "Finished." << std::endl;
}
//End of the CreateVector() function.

int main(){
CreateVector(randomloc);
}

/*
void CreateVectors(){
  CreateVector(fc);
  //CreateVector(pc);
  //CreateVector(upmus);
  //CreateVector(upmut);
  //CreateVector(cosmic);
}
*/
