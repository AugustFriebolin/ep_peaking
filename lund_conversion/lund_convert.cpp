#include "TFile.h"
#include "TTree.h"
#include "TVectorT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "ep_constants.h"

using namespace std;

const int NperFile=10000;
const int pCode=2212;

string int2str(int x)
{
    ostringstream temp;
    temp << x;
    return temp.str();
}

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "    lund_convert /path/to/input/rootfile /path/to/lund/basename\n\n"
	   << "  Aborting...\n";
      return -1;
    }

  TFile * infile = new TFile(argv[1]);
  TTree * intree = (TTree*)infile->Get("T");
  std::string basename(argv[2]);


  double theta_e, phi_e, mom_e, theta_p, phi_p, mom_p, weight;
  intree->SetBranchAddress("theta_e",&theta_e);
  intree->SetBranchAddress("phi_e",&phi_e);
  intree->SetBranchAddress("mom_e",&mom_e);
  intree->SetBranchAddress("theta_p",&theta_p);
  intree->SetBranchAddress("phi_p",&phi_p);
  intree->SetBranchAddress("mom_p",&mom_p);
  intree->SetBranchAddress("weight",&weight);

  int Nevents = intree->GetEntries();
  
  TVectorT<double> *runInfo = (TVectorT<double>*)infile->Get("runInfo");
  // runInfo[0]: Nevents, runInfo[2]=cosTheta_e_max, runInfo[3]=cosTheta_e_min
  const double Ebeam=(*runInfo)[1];
    
  intree->GetEvent(0);
  
  int fileNo=0;
  std::string outfilename=basename + "_" + int2str(fileNo) + ".lund";
  ofstream outfile(outfilename.c_str());

  // For documentation, see: https://gemc.jlab.org/gemc/html/documentation/generator/lund.html
  
  // Write first event
  // Start with the header
  //    [Number of particles] [Mass number of the target] [Atomic number of the target] [Target polarization] [Beam polarization]
  //    [Beam type (e=11, photon=22)] [Beam energy (GeV)] [Interacted nucleon (2212 or 2112)] [Process ID] [Event weight]
  outfile << "2 1 1 0. 0. 11 " << Ebeam << " 2212 1 " << weight << "\n";
  // Next write out individual particles
  //    [index] [Lifetime (ns)] [Propagate in Geant? (1 yes, 0 no)] [Particle ID] [Index of the parent] [Index of the first daughter]
  //    [Mom x (GeV)] [Mom y (GeV)] [Mom z (GeV)] [Energy (Gev)] [Mass (GeV)] [Vertex x (cm)] [Vertex y (cm)] [Vertex z (cm)]
  // Electron
  outfile << "1 1000. 1 11 0 0 "
	  << mom_e*sin(theta_e)*cos(phi_e) << " " << mom_e*sin(theta_e)*sin(phi_e) << " " << mom_e*cos(theta_e) << " "
	  << mom_e << " " << me << " 0. 0. 0.\n";
  // Nucleon
  outfile << "2 1000. 1 2212 0 0 "
	  << mom_p*sin(theta_p)*cos(phi_p) << " " << mom_p*sin(theta_p)*sin(phi_p) << " " << mom_p*cos(theta_p) << " "
	  << sqrt(mom_p*mom_p + mP*mP) << " " << mP << " 0. 0. 0.\n";
  
  for (int i=1 ; i<Nevents ; i++)
    {
      intree->GetEvent(i);

      if (i% NperFile ==0)
	{
	  cerr << "Finished writing file " << fileNo << "\n";
	  fileNo+=1;
	  cerr << "  Now beginning file " << fileNo << "\n";

	  outfile.close();
	  outfilename=basename + "_" + int2str(fileNo) + ".lund";
	  outfile.open(outfilename.c_str());
	}

      // Write out the event
        outfile << "2 1 1 0. 0. 11 " << Ebeam << " 2212 1 " << weight << "\n";
	outfile << "1 1000. 1 11 0 0 "
		<< mom_e*sin(theta_e)*cos(phi_e) << " " << mom_e*sin(theta_e)*sin(phi_e) << " " << mom_e*cos(theta_e) << " "
		<< mom_e << " " << me << " 0. 0. 0.\n";
	outfile << "2 1000. 1 2212 0 0 "
		<< mom_p*sin(theta_p)*cos(phi_p) << " " << mom_p*sin(theta_p)*sin(phi_p) << " " << mom_p*cos(theta_p) << " "
		<< sqrt(mom_p*mom_p + mP*mP) << " " << mP << " 0. 0. 0.\n";	
    }

  outfile.close();
  infile->Close();
    
  return 0;

}
