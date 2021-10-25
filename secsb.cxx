#include <iostream>
#include <cmath>
#include <random>

using namespace std;

// Constants
// ---------
static const int maxNtr        = 10000; 
static bool force_gaussians    = false;      // In gaussian toy gen, force to pick gaussian dims for subspace 

static const int ND            = 50;         // total number of considered space dimensions
static const int maxNvar       = ND;         // max dimensionality of search box
static const int maxEvents     = 10000;
static const double sqrt2      = sqrt(2.);
static const int maxClosest    = 20;
static const double RegFactor  = 1.;  // But it could be larger if using useSB true, because of flukes
static const int N2var         = 105; // maxNvar*(maxNvar-1)/2;
static const double alpha      = 0.05;

// Control variables
// -----------------
static int maxJumps            = 10;         // max number of jumps in each direction to search for new minima
static double maxBoxVolume     = 0.5;        // max volume of search box
static double syst             = 0.;         // systematics on tau!
static int maxGDLoops          = 100;        // max number of GD loops
static bool verbose            = false;      // control printouts on screen
static bool plots              = true;       // whether to plot distributions
static bool RemoveHighCorr     = false;      // high-correlation vars are removed, NSEL are kept among NAD
static bool Multivariate       = true;       // whether signal shows correlations among features
static bool debug              = false;      // couts toggle
static bool compactify         = false;      // If true, empty space in domain of features is removed
static bool fixed_gaussians    = true;       // In gaussian toy generations, controls shape of signal
static bool narrow_gaussians   = true;       // In gaussian toy generations, controls what kind of signal is injected
static int Gaussian_dims       = 15;
static double maxRho           = 0.2;
static double maxHalfMuRange   = 0.35;
static int NseedTrials         = 1;          // Number of repetition of same subspace search, for tests of clustering
static double shrinking_factor = 0.9;        // decreasing factor for steps in direction not previously taken
static double widening_factor  = 1.0;        // widening factor for steps in same direction
static double InitialLambda    = 0.5;        // step size in search
static double sidewidth        = 0.5;
static double D2kernel         = 0.04;
static double secondsidewidth = 2./3.;

// Other control variables
// -----------------------
static double bignumber      = pow(10,20.);
static double smallnumber    = pow(10,-20.);
static double epsilon        = 0.01;        // coarseness of scan of multi-D unit cube
static double InvEpsilon     = 100.;        // inverse of epsilon


void determineSB (double Smi[maxNvar], double Sma[maxNvar], double Bmi[maxNvar], double Bma[maxNvar], int Nvar) {
  double minratio = bignumber;
  double AvailableRatio[maxNvar];
  // Find minimum ratio of available space to box space, among all directions
  for (int i=0; i<Nvar; i++) {
    Smi[i] = Bmi[i]*(1+sidewidth)-Bma[i]*sidewidth;
    if (Smi[i]<0.) Smi[i] = 0.;
    Sma[i] = Bma[i]*(1+sidewidth)-Bmi[i]*sidewidth;
    if (Sma[i]>1.) Sma[i] = 1.;
    AvailableRatio[i] = (1.-(Bma[i]-Bmi[i]))/(Bma[i]-Bmi[i]);
    if (AvailableRatio[i]<minratio) minratio = AvailableRatio[i];
  }
  // Order by available ratio
  int ind[maxNvar];
  for (int i=0; i<Nvar; i++) { ind[i]=i; };
  for (int times=0; times<Nvar; times++) {
    for (int i=Nvar-1; i>0; i--) {
      if (AvailableRatio[ind[i]]<AvailableRatio[ind[i-1]]) {
// Swap indices
int tmp  = ind[i];
ind[i]   = ind[i-1];
ind[i-1] = tmp;
      }
    }
  }
  // Now AvailableRatio[ind[Nvar-1]] is the largest, AvailableRatio[ind[0]] is the smallest
  double NeededRatioPerVar;
  double CurrentFactor = 1.;
  for (int i=0; i<Nvar; i++) {
    if (AvailableRatio[ind[i]]==0) continue; // can't use this dimension
    NeededRatioPerVar = pow(2./CurrentFactor,1./(Nvar-i))-1.;
    if (AvailableRatio[ind[i]]<NeededRatioPerVar) { // use all the space available for this var
      Smi[ind[i]] = 0.;
      Sma[ind[i]] = 1.;
      CurrentFactor = CurrentFactor*(1.+AvailableRatio[ind[i]]);
      if (i<Nvar-1) NeededRatioPerVar = pow(2./CurrentFactor,Nvar-i-1)-1.; // rescaled needed ratio for the others
    } else { // We can evenly share the volume in the remaining coordinates
      double distmin = Bmi[ind[i]];
      double deltax  = Bma[ind[i]]-Bmi[ind[i]];
      if (distmin>1.-Bma[ind[i]]) { // Upper boundary is closest
distmin = 1.-Bma[ind[i]];
if (2.*distmin/deltax>=NeededRatioPerVar) {
  Smi[ind[i]] = Bmi[ind[i]]-NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bmi[ind[i]]-NeededRatioPerVar*deltax/2.));
  Sma[ind[i]] = Bma[ind[i]]+NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bma[ind[i]]+NeededRatioPerVar*deltax/2.));
} else {
  Sma[ind[i]] = 1.;
  Smi[ind[i]] = 1.-deltax*(1.+NeededRatioPerVar); // epsilon*(int)(InvEpsilon*(1.-deltax*(1.+NeededRatioPerVar)));
}
CurrentFactor = CurrentFactor*(1.+NeededRatioPerVar);
      } else { // lower boundary is closest 
if (2.*distmin/deltax>=NeededRatioPerVar) {
  Smi[ind[i]] = Bmi[ind[i]]-NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bmi[ind[i]]-NeededRatioPerVar*deltax/2.));
  Sma[ind[i]] = Bma[ind[i]]+NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bma[ind[i]]+NeededRatioPerVar*deltax/2.));
} else {
  Smi[ind[i]] = 0.;
  Sma[ind[i]] = deltax*(1.+NeededRatioPerVar); // epsilon*(int)(InvEpsilon*(deltax*(1.+NeededRatioPerVar)));
}
CurrentFactor = CurrentFactor*(1.+NeededRatioPerVar);
      }
    }
  }
  return;
}


void determinesecondSB (double Smi[maxNvar], double Sma[maxNvar], double Bmi[maxNvar], double Bma[maxNvar], int Nvar) {
  double minratio = bignumber;
  double AvailableRatio[maxNvar];
  // Find minimum ratio of available space to box space, among all directions
  for (int i=0; i<Nvar; i++) {
    Smi[i] = Bmi[i]*(1+secondsidewidth)-Bma[i]*secondsidewidth;
    if (Smi[i]<0.) Smi[i] = 0.;
    Sma[i] = Bma[i]*(1+secondsidewidth)-Bmi[i]*secondsidewidth;
    if (Sma[i]>1.) Sma[i] = 1.;
    AvailableRatio[i] = (1.-(Bma[i]-Bmi[i]))/(Bma[i]-Bmi[i]);
    if (AvailableRatio[i]<minratio) minratio = AvailableRatio[i];
  }
  // Order by available ratio
  int ind[maxNvar];
  for (int i=0; i<Nvar; i++) { ind[i]=i; };
  for (int times=0; times<Nvar; times++) {
  	 for (int i=Nvar-1; i>0; i--) {
      	 	 if (AvailableRatio[ind[i]]<AvailableRatio[ind[i-1]]) {
// Swap indices
int tmp  = ind[i];
ind[i]   = ind[i-1];
ind[i-1] = tmp;
      }
    }
  }
  // Now AvailableRatio[ind[Nvar-1]] is the largest, AvailableRatio[ind[0]] is the smallest
  double NeededRatioPerVar;
  double CurrentFactor = 1.;
  for (int i=0; i<Nvar; i++) {
    if (AvailableRatio[ind[i]]==0) continue; // can't use this dimension
    NeededRatioPerVar = pow(3./2./CurrentFactor,1./(Nvar-i))-1.;
    if (AvailableRatio[ind[i]]<NeededRatioPerVar) { // use all the space available for this var
      Smi[ind[i]] = 0.;
      Sma[ind[i]] = 1.;
      CurrentFactor = CurrentFactor*(1.+AvailableRatio[ind[i]]);
      if (i<Nvar-1) NeededRatioPerVar = pow(3./2./CurrentFactor,Nvar-i-1)-1.; // rescaled needed ratio for the others
    } else { // We can evenly share the volume in the remaining coordinates
      double distmin = Bmi[ind[i]];
      double deltax  = Bma[ind[i]]-Bmi[ind[i]];
      if (distmin>1.-Bma[ind[i]]) { // Upper boundary is closest
distmin = 1.-Bma[ind[i]];
if (2.*distmin/deltax>=NeededRatioPerVar) {
  Smi[ind[i]] = Bmi[ind[i]]-NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bmi[ind[i]]-NeededRatioPerVar*deltax/2.));
  Sma[ind[i]] = Bma[ind[i]]+NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bma[ind[i]]+NeededRatioPerVar*deltax/2.));
} else {
  Sma[ind[i]] = 1.;
  Smi[ind[i]] = 1.-deltax*(1.+NeededRatioPerVar); // epsilon*(int)(InvEpsilon*(1.-deltax*(1.+NeededRatioPerVar)));
}
CurrentFactor = CurrentFactor*(1.+NeededRatioPerVar);
      } else { // lower boundary is closest 
if (2.*distmin/deltax>=NeededRatioPerVar) {
  Smi[ind[i]] = Bmi[ind[i]]-NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bmi[ind[i]]-NeededRatioPerVar*deltax/2.));
  Sma[ind[i]] = Bma[ind[i]]+NeededRatioPerVar*deltax/2.; // epsilon*(int)(InvEpsilon*(Bma[ind[i]]+NeededRatioPerVar*deltax/2.));
} else {
  Smi[ind[i]] = 0.;
  Sma[ind[i]] = deltax*(1.+NeededRatioPerVar); // epsilon*(int)(InvEpsilon*(deltax*(1.+NeededRatioPerVar)));
}
CurrentFactor = CurrentFactor*(1.+NeededRatioPerVar);
      }
    }
  }
  return;
}


int main () {
 double Nvar = 12.;
 double Sidemin[maxNvar] = {0.};
 double Sidemax[maxNvar] = {0.};
 double secondsidemin[maxNvar] = {0.};
 double secondsidemax[maxNvar] = {0.};
 double Blockmin[maxNvar];
 double Blockmax[maxNvar];
 double lbound, rbound, adjust;

 default_random_engine generator;
 uniform_real_distribution<double> distribution(0.0,1.0);

 for (int k=0; k<maxNvar; k++) {	// modify these parameters to modify the size of the signal box
					// to make a simple verification of the validity of the code i tried the simplest case, where the box
					// is exactly a cube in Nvar dimensions, so more attempts in less simple cases are needed to confirm.
 	 //Blockmin[k]=0.13;
	 //Blockmax[k]=0.45;
 	 lbound = distribution(generator);
 	 rbound = distribution(generator);
	 if (lbound>rbound) {adjust = lbound; lbound = rbound; rbound = adjust;} // con questo pezzo di codice ci si assicura che lbound sia 											    sempre minore di rbound
	 else if (lbound == rbound) rbound += 0.0000001;	
	 Blockmin[k] = lbound;
	 Blockmax[k] = rbound;
 }	
 determineSB(Sidemin,Sidemax,Blockmin,Blockmax,Nvar);
 double VolumeOrig      = 1.;
 double SidebandsVolume = 1.;
 double SecondsidebandsVolume = 1.;
 for (int k=0; k<Nvar; k++) {
  	 VolumeOrig      *= fabs(Blockmax[k]-Blockmin[k]);
  	 SidebandsVolume *= fabs(Sidemax[k]-Sidemin[k]);
 }
 determinesecondSB(secondsidemin, secondsidemax, Sidemin, Sidemax, Nvar);
 for (int k=0; k<Nvar; k++) {
  	 SecondsidebandsVolume *= fabs(secondsidemax[k]-secondsidemin[k]);
 } 

 cout << "box volume is: " << VolumeOrig << "\nsideband volume is: " << SidebandsVolume << "\n2nd sideband volume is: " << SecondsidebandsVolume << endl;

 for (int k=0; k<Nvar; k++) {
	 cout << "dimension " << k << " : boundaries" << endl;
         cout << "Bmin = " << Blockmin[k] << "; Bmax = " << Blockmax[k] << endl;
         cout << "Sidemin = " << Sidemin[k] << "; Sidemax = " << Sidemax[k] << endl;
         cout << "secsmin = " << secondsidemin[k] << "; secsbmax = " << secondsidemax[k] << endl;
	 cout << "******************" << endl;
 }
 return 0;
}
