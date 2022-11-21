#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

using namespace std;

const double mP=0.938272;
const double me=0.000511;
const double muP=2.79;
const double alpha = 0.0072973525;
const double hbarc=0.1973; // GeV*fm
const double nbGeVSq = hbarc*hbarc *1.E7;

double Gdip(double QSq)
{
  const double lamSq = 0.71;
  return TMath::Power(1. + QSq/lamSq,-2.);
}

double GE(double QSq)
{
  return Gdip(QSq);
}

double GM(double QSq)
{
  return muP*Gdip(QSq);
}

double dSigdOmegaRos(double E1, double cosTheta)
{
  double E3 = E1*mP/(mP+E1*(1.-cosTheta));
  double cosSqThetaOver2=0.5*(1.+cosTheta);
  double sinSqThetaOver2=0.5*(1.-cosTheta);

  double QSq = 2.*mP*(E1-E3);
  double tau = QSq/(4.*mP*mP);
  double epsilon = 1./(1. + 2.*(1.+tau)*sinSqThetaOver2/cosSqThetaOver2);
  
  double sigma_Mott = nbGeVSq*alpha*alpha*E3*cosSqThetaOver2/(4.*E1*E1*E1*sinSqThetaOver2*sinSqThetaOver2);
  
  return sigma_Mott * (1./(1.+tau)) * (GE(QSq)*GE(QSq) + tau/epsilon * GM(QSq)*GM(QSq));
}

const int requiredArgs=4;

int main(int argc, char ** argv)
{
  if (argc < requiredArgs)
    {
      cerr << "Wrong number of arguments. Usage is:\n"
           << "    generator /path/to/output/file [Nevents] [Beam energy (GeV)] [optional arguments...]\n"
           << "        Optional arguments:\n"
           << "            -t: minimum theta [deg.]\n"
           << "            -T: maximum theta [deg.]\n"
           << "\n\n";

      return -1;
    }

  // Set up the output file and the tree
  const int Nevents=atoi(argv[2]);
  TFile * outfile = new TFile(argv[1],"RECREATE");
  TTree * outtree = new TTree("T","ep generator tree");
  const double E1=atof(argv[3]);

  double cosTheta_e_max = TMath::Cos(1.*M_PI/180.);
  double cosTheta_e_min = TMath::Cos(179.*M_PI/180.);

  // optional arguments
  int c;
  while ((c = getopt (argc-requiredArgs+1, &argv[requiredArgs-1], "t:T:")) != -1)
    {
      switch (c)
        {
        case 't':
          cosTheta_e_max = cos(atof(optarg)*M_PI/180.);
          break;
        case 'T':
          cosTheta_e_min = cos(atof(optarg)*M_PI/180.);
          break;
	default:
          cerr << "Optional argument not found. Aborting...\n\n";
          return -1;
        }
    }

  // Set up a random number generator
  TRandom3 * myRand = new TRandom3(0); 
 
  // Tree info
  double mom_e,theta_e,phi_e,mom_p, theta_p, phi_p, weight;
  outtree->Branch("mom_e",&mom_e, "mom_e/D");
  outtree->Branch("theta_e",&theta_e, "theta_e/D");
  outtree->Branch("phi_e",&phi_e, "phi_e/D");
  outtree->Branch("mom_p",&mom_p, "mom_p/D");
  outtree->Branch("theta_p",&theta_p, "theta_p/D");
  outtree->Branch("phi_p",&phi_p, "phi_p/D");
  outtree->Branch("weight",&weight, "weight/D");
  
  // Loop over events
  for (int event=0 ; event < Nevents ; event++)
    {
      if (event %100000==0)
	{
	  cerr << "Working on event " << event << " out of " << Nevents << "...\n";
	}
      
      // First step is to generate an electron scattering angle
      double cosTheta_e = cosTheta_e_min + (cosTheta_e_max - cosTheta_e_min)*myRand->Rndm();
      theta_e = TMath::ACos(cosTheta_e);
      
      // This can be used to establish kinematics
      mom_e=E1*mP / (mP+E1*(1.-cosTheta_e));
      phi_e = myRand->Rndm()*2.*M_PI - M_PI;

      TVector3 p1(0.,0.,E1);
      TVector3 p3(mom_e*TMath::Sin(theta_e)*TMath::Cos(phi_e), mom_e*TMath::Sin(theta_e)*TMath::Sin(phi_e), mom_e*cosTheta_e);
      TVector3 p4 = p1 - p3;
      theta_p = p4.Theta();
      phi_p = p4.Phi();
      mom_p = p4.Mag();

      //Weight is: d^N sig / P(cosTheta, phi)
      weight = dSigdOmegaRos(E1, cosTheta_e) * (2.*M_PI) *  (cosTheta_e_max - cosTheta_e_min);

      // Sanitize
      if (weight<0.)
	weight=0.;
      
      outtree->Fill();
    }

  outfile->Write();
  outfile->Close();

  return 0;
}
