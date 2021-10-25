#include <iostream>
#include <cmath>

// Signal vs square root of sum of background and corrective factor
double SvsBKG (int Non, int Noff, double SigboxVolume, double secondSBVolume) {
  if (Non==0 || SigboxVolume==0 || secondSBVolume==0) return 0;
  // if (Non==0 || Noff==0 || SigboxVolume==0 || secondSBVolume==0) return 0;
  // if (Non>0 && Noff==0) return 0.;
  double r = SigboxVolume/secondSBVolume;
  double SvB = Non/sqrt((Noff+1.)*r);
  return SvB;
}

// SvB with no sideband events, in case its needed 
double SvsBKGnosb(int Non, int Ntot, double SigboxVolume, double secondSBVolume) {
  if (Non==0|| SigboxVolume==0 || secondSBVolume==0) return 0.; 
  // if (Non>0 && Ntot-Non==0) return 0.;
  int Noff=Ntot-Non;
  double r = SigboxVolume/secondSBVolume;
  double SvB = Non/sqrt((Noff+1.)*r);
  return SvB;
}
