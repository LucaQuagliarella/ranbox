#include <iostream>
#include <cmath>
#include <random>
#include <vector>

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

void secsbrandeventsroot() {
 int Nvar = 10;
 double Sidemin[maxNvar] = {0.};
 double Sidemax[maxNvar] = {0.};
 double secondsidemin[maxNvar] = {0.};
 double secondsidemax[maxNvar] = {0.};
 double Blockmin[maxNvar];
 double Blockmax[maxNvar];
 double VolumeOrig      = 1.;
 double SidebandsVolume = 1.;
 double SecondsidebandsVolume = 1.;


 // make histograms for the distribution of the volumes
 long int ratiocycles = 10000;
 int NEvents = 1000000;
 /*
 double boxboxratio[ratiocycles];
 double sbboxratio[ratiocycles];
 double secsbboxratio[ratiocycles];
 */
 double EV [NEvents] [Nvar];
 double event;
 double bbratio, sbbratio, secsbbratio;
 double lbound,rbound,adjust;

 default_random_engine generator;
 uniform_real_distribution<double> distribution(0.0,1.0);		// we want volume of 10^-2 and 10^-3
									// throw random on a subset of dimensions (e.g. 6 random 4 set)

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
 TCanvas *c1 = new TCanvas("c1", "volume ratio plots", 600, 400);
 TCanvas *c2 = new TCanvas("c2", "volume ratio plots (second sideband)", 600, 400);
 TCanvas *c3 = new TCanvas("c3", "R1 dist", 600, 400);
 TCanvas *c4 = new TCanvas("c4", "R2 dist", 600, 400);
 TH1* h1 = new TH1D("h1", "volume ratio distribution, first sideband", 100, 0.0, 2.5);
 TH1* h2 = new TH1D("h2", "volume ratio distribution, second sideband", 100, 0.0, 3.5);
 TH1* h3 = new TH1D("h3", "R1", 100, 0.0, 5);
 TH1* h4 = new TH1D("h4", "R2", 100, 0.0, 5);
 
 // now I have NEvents in Nvar dimensions in memory

 for (int i=0; i<ratiocycles; i++) { 
	// reset values for each cycle
	VolumeOrig      = 1.;
	SidebandsVolume = 1.;
 	SecondsidebandsVolume = 1.;
	Ns = 0;
	Nsb = 0;
	N2sb = 0;
	for (int k=0; k<maxNvar; k++) {
		// code to randomize blockmin/max values, uniformly between 0 and 1
 	 lbound = distribution(generator);
 	 rbound = distribution(generator);
	 if (lbound>rbound) {adjust = lbound; lbound = rbound; rbound = adjust;} // con questo pezzo di codice ci si assicura che lbound sia 											    sempre minore di rbound
	 else if (lbound == rbound) rbound += 0.0000001;	
	 Blockmin[k] = lbound;
	 Blockmax[k] = rbound;
 	}
	// construct sidebands
	determineSB(Sidemin,Sidemax,Blockmin,Blockmax,Nvar);
 	determinesecondSB(secondsidemin, secondsidemax, Sidemin, Sidemax, Nvar);

        // now we look for events in each our boxes
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

	/*for (int k=0; k<Nvar; k++) {
	 // code to randomize blockmin/max values, uniformly between 0 and 1
	 cout << Blockmin[k] << endl;
	 cout << Blockmax[k] << endl;
	 cout << Sidemin[k] << endl;
	 cout << Sidemax[k] << endl;
	 cout << secondsidemin[k] << endl;
	 cout << secondsidemax[k] << endl;
 	}*/

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
	 h3->Fill(R1);
	 cout << "iteration " << i << ", R1 = " << R1 << endl;
	}
 	if (N2sb !=0) { 
	 R2 = Ns/(N2sb*(SecondsidebandsVolume-2.*VolumeOrig)/VolumeOrig);
	 h4->Fill(R2);
	 cout << "iteration " << i << ", R2 = " << R2 << endl;
	}
	// store ratios
	// sbboxratio[i] = sbbratio;
	// secsbboxratio[i] = secsbbratio;
        if (sbbratio < 2 ) 
	cout << "iteration " << i << ", ratio[2] = " << sbbratio << endl;
	else if (secsbbratio < 3) 
	cout << "iteration " << i << ", ratio[3] = " << secsbbratio << endl;

	// fill histograms
	h1->Fill(sbbratio);
        h2->Fill(secsbbratio);
	
 }


 // Draw the histograms
 h1->Draw();
 c1->cd();
 h2->Draw();
 c2->cd();
 h3->Draw();
 c3->cd();
 h4->Draw();
 c4->cd();

 return;
}
