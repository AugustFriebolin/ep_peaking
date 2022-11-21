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

double delta_hard(double E1, double cosTheta)
{
  // We first need to establish QSq
  double E3 = E1*mP/(mP+E1*(1.-cosTheta));
  double QSq = 2.*mP*(E1-E3);
  
  // First evaluate the vacuum polarization, (only including electron term)
  double delta_vp = 1./(3.*M_PI) * ( -5./3. + TMath::Log(QSq/(me*me)));

  return 2.*alpha*(-(0.75/M_PI)*TMath::Log(QSq/me*me) + 1./M_PI - delta_vp);
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
  
  // Optional arguments
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
  double mom_e,theta_e,phi_e,mom_p, theta_p, phi_p, E1_eff, weight;
  outtree->Branch("mom_e",&mom_e, "mom_e/D");
  outtree->Branch("theta_e",&theta_e, "theta_e/D");
  outtree->Branch("phi_e",&phi_e, "phi_e/D");
  outtree->Branch("mom_p",&mom_p, "mom_p/D");
  outtree->Branch("theta_p",&theta_p, "theta_p/D");
  outtree->Branch("phi_p",&phi_p, "phi_p/D");
  outtree->Branch("E1_eff",&E1_eff, "E1_eff/D");
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
      
      // This gives us enough info to generate the energy radiated by the incoming electron
      double mom_e_no_rad = E1*mP/(mP+E1*(1.-cosTheta_e));
      double lambda_tilde_e = (alpha/M_PI) * (TMath::Log(4.*E1*E1/(me*me))-1. + 2.*TMath::Log(E1/mom_e_no_rad) + TMath::Log(0.5*(1.-cosTheta_e)));

      // Randomly determine initial radiation, get updated beam energy.
      // Goes like E_ISR ^ (lambda - 1)
      /*
	int[0,E] N * E^(lam - 1) dE = 1

	CDF(E) = N * (1/lam) E^lam

	CDF(Emax) = 1 = N * (1/lam) Emax^lam

	            N = lam / Emax^lam

		    ---->  PDF = lam E^{lam-1} / Emax^lam   =   (lam/E) * (E/Emax)^lam

	CDF(E) = (E/Emax)^lam

	E = Emax*(r)^{1/lam}

       */

      double E_rad_e = E1 * TMath::Power(myRand->Rndm(),1./lambda_tilde_e);
      E1_eff = E1 - E_rad_e;

      // Establish un-radiated out-going particle vectors
      double E3_eff = E1_eff*mP / (mP+E1_eff*(1.-cosTheta_e));
      phi_e = myRand->Rndm()*2.*M_PI - M_PI;
      TVector3 p1_eff(0.,0.,E1_eff);
      TVector3 p3_eff(E3_eff*TMath::Sin(theta_e)*TMath::Cos(phi_e), E3_eff*TMath::Sin(theta_e)*TMath::Sin(phi_e), E3_eff*cosTheta_e);
      TVector3 p4_eff = p1_eff - p3_eff;
      theta_p = p4_eff.Theta();
      phi_p = p4_eff.Phi();
      double mom4_eff = p4_eff.Mag();
      double E4_eff = TMath::Sqrt(mom4_eff*mom4_eff + mP*mP);

      // Calculate lepton and proton radiation
      double lambda_tilde_ep = (alpha/M_PI) * (TMath::Log(4.*E3_eff*E3_eff/(me*me)) -1. + 2.*TMath::Log(E1/mom_e_no_rad) + TMath::Log(0.5*(1.-cosTheta_e)));

      // Sanitize
      double E_rad_ep = 0.;

      if (lambda_tilde_ep > 0.) 
	E_rad_ep = E3_eff * TMath::Power(myRand->Rndm(),1./lambda_tilde_ep); 

      mom_e = E3_eff - E_rad_ep;

      /*
      if ((mom_e < 0.)||(isnan(mom_e)))
	{
	  cerr << E1_eff << " " << theta_e*180./M_PI << " " << lambda_tilde_ep << " " << E3_eff << " " << E_rad_ep << " " << mom_e << "\n";
	  //return -2;
	}
      */
      
      double lambda_tilde_p = (alpha/M_PI) * ((E4_eff/mom4_eff)*TMath::Log((E4_eff +mom4_eff)/(E4_eff-mom4_eff)) - 2.);
      double E_rad_p = (E4_eff-mP) * TMath::Power(myRand->Rndm(),1./lambda_tilde_p);
      double E4 = E4_eff - E_rad_p;
      mom_p = TMath::Sqrt(E4*E4 - mP*mP);
      
      // We have now filled all the kinematics. We only have to do the weight
      /*
	Weight is: d^N sig / P(cosTheta, phi, E_rad_e, E_rad_ep, E_rad_p)

	w = {dSig/dOmega} * (1-delta_hard) / (sqrt terms) * (2.*pi * delta cosTheta)
	
       */
      weight = dSigdOmegaRos(E1_eff, cosTheta_e) * (1. - delta_hard(E1_eff, cosTheta_e)) * 2.*M_PI * (cosTheta_e_max - cosTheta_e_min)
	/ (TMath::Power(TMath::Sqrt(E1*mom_e),lambda_tilde_e)
	   *TMath::Power(TMath::Sqrt(E1*mom_e),lambda_tilde_ep)
	   *TMath::Power(TMath::Sqrt(mP*mom_p),lambda_tilde_p));
      
      // Sanitize
      if (weight<0.)
	weight=0.;
            
      outtree->Fill();
    }

  outfile->Write();
  outfile->Close();

  return 0;
}
