#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

#include "Constants.h"
#include "Sidebands.h"


using namespace std;

/* 

 takes as argument:
 - number of random generated events
 - number of dimensions
 - number of dimensions for which the dimension is set
 - number of boxes thrown
 - lenght of the set side of the box
 
*/

bool searchset (int arr[], int size, int x) {
	int is;
	for (is=0; is < size; is++) if (arr[is] == x) return 1;
	return 0; 
}

void setvolumeroot( int NEvents = 10000, int Nvar = 10, int setdim = 4, int ratiocycles = 10000, double set = 1e-3) {
 double Sidemin[maxNvar] = {0.};
 double Sidemax[maxNvar] = {0.};
 double secondsidemin[maxNvar] = {0.};
 double secondsidemax[maxNvar] = {0.};
 double Blockmin[maxNvar];
 double Blockmax[maxNvar];
 double VolumeOrig      = 1.;
 double SidebandsVolume = 1.;
 double SecondsidebandsVolume = 1.;
 int dimensionofset;


 // make histograms for the distribution of the volumes
 double EV [NEvents] [Nvar];
 double event;
 double bbratio, sbbratio, secsbbratio;
 double lbound,rbound,adjust;

 default_random_engine generator;
 uniform_real_distribution<double> distribution(0.0,1.0);		// we want volume of 10^-2 and 10^-3
									// throw random on a subset of dimensions (e.g. 6 random 4 set)
									// the dimension where the side is set is picked at random (right?)
 default_random_engine setter;
 uniform_int_distribution<int> setterdist(0,Nvar-1);

 int indexstore[setdim];
 int checkset;

 for (int i=0; i<setdim; i++) indexstore[setdim] = -1;

 checkset = setterdist(setter);

 int Ns=0, Nsb=0, N2sb=0;
 int checkNs, checkNsb, checkN2sb;
 double R1, R2;

 // events generation, we use uniform distribution.

 for (int i=0; i<NEvents; i++) {
	 for (int k=0; k<Nvar; k++) { 
		 event = distribution(generator);
		 EV[i][k] = event;
	 }
 }


 // make root histogram
 TCanvas *c1 = new TCanvas("c1", "canv1", 600, 400);
 TCanvas *c2 = new TCanvas("c2", "canv2", 600, 400);
 TH1* h1 = new TH1D("h1", "R1", 100, 0.0, 10);
 TH1* h2 = new TH1D("h2", "R2", 100, 0.0, 10);
 
 // now I have NEvents in Nvar dimensions in memory

 for (int i=0; i<ratiocycles; i++) { 
	// reset values for each cycle
	VolumeOrig      = 1.;
	SidebandsVolume = 1.;
 	SecondsidebandsVolume = 1.;
	Ns = 0;
	Nsb = 0;
	N2sb = 0;

	for (int k = 0; k<setdim; k++) {								// pick the dimensions to initialize differently
		while (searchset(indexstore, setdim, checkset)) checkset = setterdist(setter); 	// avoids a dimension which is already set to be repeated
													// fin quando trova l'elemento ripete per averne uno nuovo
		indexstore[k] = checkset;
	}
	sort(indexstore, indexstore+setdim);
	for (int k=0; k<maxNvar; k++) {
	 // code to randomize blockmin/max values, uniformly between 0 and 1, but setting the width of 4 dimension to 1e-3
 	
		 lbound = distribution(generator);
 		 if (k!=indexstore[k]) rbound = distribution(generator);
		 else {
			if ((lbound + 1e-3)>1) {
			  adjust = 1. - lbound;		// this maintains the boundary effect of reducing volume 
			  rbound = lbound + adjust;
			}
		        rbound = lbound + set;
		 }
		 if (lbound>rbound) {adjust = lbound; lbound = rbound; rbound = adjust;} 
		 else if (lbound == rbound) rbound += 0.0000001;	
		 Blockmin[k] = lbound;
		 Blockmax[k] = rbound;
 	}
	// construct sidebands
	determineSB(Sidemin,Sidemax,Blockmin,Blockmax,Nvar);
 	determinesecondSB(secondsidemin, secondsidemax, Sidemin, Sidemax, Nvar);

        // now we look for events in each of our boxes
	for (int p = 0; p<NEvents; p++) {
	 checkNs = 0;
	 checkNsb = 0;
	 checkN2sb = 0;
		 for (int r = 0; r<Nvar; r++) {
			 if (EV[p][r] > Blockmin[r] && EV[p][r] < Blockmax[r]) checkNs++;
			 if ((EV[p][r] > Sidemin[r] && EV[p][r] < Sidemax[r])) checkNsb++;
			 if ((EV[p][r] > secondsidemin[r] && EV[p][r] < secondsidemax[r])) checkN2sb++;
		 } 
	 // to increment in a box we need the signal to be inside the boundaries for every dimension
	 if (checkNs == Nvar) Ns++;		
	 else if (checkNsb == Nvar) Nsb++;
	 else if (checkN2sb == Nvar) N2sb++;
	}
	// calculate volume
	for (int k=0; k<Nvar; k++) {
  	 	VolumeOrig      *= fabs(Blockmax[k]-Blockmin[k]);
  	 	SidebandsVolume *= fabs(Sidemax[k]-Sidemin[k]);
  	 	SecondsidebandsVolume *= fabs(secondsidemax[k]-secondsidemin[k]);
 	}
	// cout << "VOLUME: " << VolumeOrig << " and " << SidebandsVolume << endl;
	// calculate ratios
	sbbratio = SidebandsVolume/VolumeOrig;
	secsbbratio = SecondsidebandsVolume/VolumeOrig;
 	if (Nsb !=0) {
	 R1 = Ns/(Nsb*(SidebandsVolume-VolumeOrig)/VolumeOrig);
	 h1->Fill(R1);
	 cout << "iteration " << i << ", R1 = " << R1 << endl;
	}
 	if (N2sb !=0) { 
	 R2 = Ns/(N2sb*(SecondsidebandsVolume-2.*VolumeOrig)/VolumeOrig);
	 h2->Fill(R2);
	 cout << "iteration " << i << ", R2 = " << R2 << endl;
	}
        // if (sbbratio < 2 ) 
	// cout << "iteration " << i << ", ratio[2] = " << sbbratio << endl;
	// else if (secsbbratio < 3) 
	// cout << "iteration " << i << ", ratio[3] = " << secsbbratio << endl;
 }


 // Draw the histograms
 h1->Draw();
 c1->cd();
 h2->Draw();
 c2->cd();

 return;
}
