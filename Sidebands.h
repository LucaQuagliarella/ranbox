#include <iostream>
#include <cmath>

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

