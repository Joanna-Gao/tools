//Include C++
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

//Include ROOT
#include <TFile.h>
#include <TTree.h>


//WCIDDiameter = 70.8*m; // = 74m - 2*(60cm ID wall + 1m OD)
//WCIDHeight = 54.8*m; // = 60m - 2*(60cm ID wall + 2m OD)

#define PI 3.141592654
#define SEED 12345

#define IDRadius 3540 // in cm
#define IDHalfHeight  2740 // in cm
#define ODRadius 3700 // in cm
#define ODHalfHeight 3000 //in cm
#define Deadspace 600 // in cm
#define LargeCircleRadius 5000 // in cm
#define LargeCircleOffset 200 // in cm
#define FiducialCut 200 // in cm

// G4double ekin = energy - mass;
#define muMass 105.6583745// mu- mass (in MeV)
#define muEnergyCerenkovThreshold 159.7404473 //minimum energy to produce Cerenkov radiation (including mass)
#define muRadiativeLoss 2.3 // 2MeV/cm loss for muon (in MeV) 0.3 MeV per cm Cerenkov loss



enum EventType {fc, pc, upmus, upmut, cosmic};

std::default_random_engine generator(SEED);
/*
// generates a point within the fiducial cut from ID barrel walls^2
std::uniform_real_distribution<double> distributionPosRIDFiducial(0,pow(IDRadius-FiducialCut,2));
// generates a point within the ID radius^2
std::uniform_real_distribution<double> distributionPosRID(0,pow(IDRadius,2));
// generates a point within the fiducial cut from ID endcaps
std::uniform_real_distribution<double> distributionPosZIDFiducial(-(IDHalfHeight-FiducialCut), (IDHalfHeight-FiducialCut) );
// generates a point within the ID height
std::uniform_real_distribution<double> distributionPosZID(-IDHalfHeight , IDHalfHeight );
// generates an angle from 0 to 2*pi
std::uniform_real_distribution<double>distributionTheta(0.0,2*PI);
// generates a point within the OD barrel walls^2
std::uniform_real_distribution<double> distributionPosROD(pow(IDRadius+Deadspace,2),pow(ODRadius,2));
// generates a point on a large circle (for cosmics and upmu)
std::uniform_real_distribution<double> distributionPosRLarge(0,pow(LargeCircleRadius,2));
// generates a random number for probabilities
std::uniform_real_distribution<double> distributionProbability(0,1);

// generates energy for FC and stopping muon events
std::uniform_real_distribution<double> lessEnergy(0.01, muRadiativeLoss-0.1);
// generates energy for PC and through-going muon and cosmic events
std::uniform_real_distribution<double> moreEnergy(muRadiativeLoss+0.1, muRadiativeLoss+3);
*/
double Random(double a, double b){
	std::uniform_real_distribution<double> quick(a,b);
	return quick(generator);
}
double LessEnergy(){
	return Random(0.01, muRadiativeLoss-0.1);
}
double MoreEnergy(){
	return Random(muRadiativeLoss+0.1, muRadiativeLoss+3);
}
// function to generate a point inside a cylinder. height is the half height and pos is the output vertex
void GeneratePointInCylinder(double rad, double height, double Pos[3]){
	double genR = Random(0., pow(height,2));
	double genTheta = Random(0., 2*PI);
	Pos[0] = sqrt(genR)*cos(genTheta);
	Pos[1] = sqrt(genR)*sin(genTheta);
	Pos[2] = Random(-height, height);

}
// function to generate a point on the surface of a cylinder. height is the half height, pos is the output vertex
void GeneratePointOnCylinder(double rad, double height, double Pos[3]){
	// generate point b on the edge of the ID
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
std::string GetOutFileName(EventType ev){

	std::string evTypeName;

	switch (ev){
		case fc: evTypeName = "FC"; break;
		case pc: evTypeName = "PC"; break;
		case upmus: evTypeName = "UPMUS"; break;
		case upmut: evTypeName = "UPMUT"; break;
		case cosmic: evTypeName = "COSMIC"; break;
		default: std::cerr << "Unknown event type. Please select from: fc, pc, upmus, upmut, cosmic." << std::endl; exit(1);
	}

	std::string out = "VectorFile" + evTypeName + ".dat";
	return out;
}


// Function to get the length between two vectors and output a normalised direction from the first point to the second
double GetL(double aPos[3], double bPos[3], double dir[3] ){

	double l = sqrt( pow(bPos[0] - aPos[0], 2) + pow(bPos[1] - aPos[1], 2) + pow(bPos[2] - aPos[2], 2)  );
	dir[0] = (bPos[0] - aPos[0])/l;
	dir[1] = (bPos[1] - aPos[1])/l;
	dir[2] = (bPos[2] - aPos[2])/l;

	return l;
}

// Function to generate a random position for each event
void GenerateVertex(EventType ev, double aPos[3], double dir[3], double &energy ){

	// The idea is to select two point and draw a line between them in order to generate a specific kind of event. e.g. pick a point above and below the detector to generate a cosmic muon.

	//aPos[3]; // coordinates of first point
	double bPos[3]; // coordinates of second point
	//dir[3]; // normalised direction between two points
	double length; // distance between two points

	//double genThetaA = distributionTheta(generator);
	//double genThetaB = distributionTheta(generator);

	//double genR;

	if (ev == fc || ev == pc){
		// generate point a inside the fiducial volume
		GeneratePointInCylinder(IDRadius-FiducialCut, IDHalfHeight-FiducialCut, aPos);
		
		// generate point b on the edge of the ID
		GeneratePointOnCylinder(IDRadius, IDHalfHeight, bPos);

		// get distance and direction between the two points
		length = GetL(aPos, bPos, dir); // dir now filled

		if (ev == fc){
			energy = muEnergyCerenkovThreshold + length*LessEnergy();
		}
		else if (ev == pc){
			energy = muEnergyCerenkovThreshold + length*MoreEnergy();
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
			energy = muEnergyCerenkovThreshold + length*LessEnergy();
		}
		else if (ev == upmut){
			energy = muEnergyCerenkovThreshold + length*MoreEnergy();
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

	int nParticles = 100;


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

	//Output the number of particles processed, for every multiple of 1000
		if (i % 1000 == 0)
			{
			std::cout << "Number of particles in data file: " << i << std::endl;
			}


		}

	//Close the output file
	outputFile.close();
	std::cout << "Finished." << std::endl;
}
//End of the CreateVector() function.


void CreateVectors(){
	CreateVector(fc);
	CreateVector(pc);
	CreateVector(upmus);
	CreateVector(upmut);
	CreateVector(cosmic);
}
