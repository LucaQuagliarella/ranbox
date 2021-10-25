#include <iostream>
#include <cmath>
#include <random>

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


void secsbrandroot(int Nvar = 12) {
// int main() {
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
 long int ratiocycles = 100000;
 double bbratio, sbbratio, secsbbratio;
 double lbound,rbound,adjust;

 default_random_engine generator;
 uniform_real_distribution<double> distribution(0.0,1.0);

 // make root histogram
 TCanvas *c1 = new TCanvas("c1", "volume ratio plots", 600, 400);
 TCanvas *c2 = new TCanvas("c2", "volume ratio plots (second sideband)", 600, 400);
 // TCanvas *c3 = new TCanvas("c3", "lb rb scatterplot", 600, 400);
 TH1* h1 = new TH1D("h1", "volume ratio distribution, first sideband", 100, 0.0, 2.5);
 TH1* h2 = new TH1D("h2", "volume ratio distribution, second sideband", 100, 0.0, 3.5);

 for (int i=0; i<ratiocycles; i++) { 
	// reset values for each cycle
	VolumeOrig      = 1.;
	SidebandsVolume = 1.;
 	SecondsidebandsVolume = 1.;
	for (int k=0; k<maxNvar; k++) {
		// code to randomize blockmin/max values, uniformly between 0 and 1
 	 lbound = distribution(generator);
 	 rbound = distribution(generator);
	 if (lbound>rbound) {adjust = lbound; lbound = rbound; rbound = adjust;} // con questo pezzo di codice ci si assicura che lbound sia 											    sempre minore di rbound
	 else if (lbound == rbound) rbound += 0.0000001;	
	 Blockmin[k] = lbound;
	 Blockmax[k] = rbound;
	 // if (i % 2)    Blockmax[k] = 1.;	// tentativo estremo di forzare la distribuzione
 	}
	// construct sidebands
	determineSB(Sidemin,Sidemax,Blockmin,Blockmax,Nvar);
 	determinesecondSB(secondsidemin, secondsidemax, Sidemin, Sidemax, Nvar);
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
        if (sbbratio < 2 || secsbbratio < 3) 
	cout << "iteration " << i << ", ratio[2] = " << sbbratio << ", ratio[3] = " << secsbbratio << endl;
	// fill histograms
	h1->Fill(sbbratio);
        h2->Fill(secsbbratio);
 }

 // Draw the histograms
 h1->Draw();
 c1->cd();
 h2->Draw();
 // c1->Update();
 c2->cd();
 // return 0;
 return;
}
