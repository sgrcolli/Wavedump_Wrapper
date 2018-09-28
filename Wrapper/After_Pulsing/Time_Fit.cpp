// Project include files.
#include "ns.h"

// Standard library include files.
#include <cmath>
#include <cstdlib>

// ROOT include files.
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TObject.h" 
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSpectrum.h"

//=========================================================================================================================================

int main(int argc,char **argv){


	// ******************
	// * Initialization *
	// ******************
	
	// Set up random number generator.
	randomSeedTime();
	
			
	int channel[4]={0,0,0,0};
	char answer;
	char histname[200]= "";
	int test;
	int Gain[4]={0,0,0,0};
	
	
	//Read in the HV data ====================================================================================
	string hvfile = "../HVScan.txt";
	ifstream file(hvfile.c_str());
	string hvdat;
	
	vector<int> PMT_number(125,0), HV(125,0);
	vector<vector<int>> HVstep;
	vector<int> step(5,0);
	for (int i=0; i<125; i++)
		HVstep.push_back(step);
	
	for (int i=0; i<125; i++){
		for (int j=0; j<7; j++){
			file >> hvdat;
			int pmt_info =atof(hvdat.c_str());
			if (j==0){
				PMT_number[i]=pmt_info;
			}
			if (j!=0 && j!=6){
				HVstep[i][j-1]=pmt_info;
			}
			if (j==6){
				HV[i]=pmt_info;
			}
			//printf("j %d, val %d \n",j,pmt_info);
		}
	}
	//========================================================================================================
	
	while(answer!='Y'&& answer!='y'){
		
		//Determing the PMT number and the applied Voltage=====================================================
		cout << "Input the PMT number in Channel 0 \n" ;
		cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
		cin  >> channel[0]; 
		cout <<endl;
		
		cout << "Input the PMT VOLTAGE in Channel 0 \n" ;
		cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
		cin  >> Gain[0]; 
		cout <<endl;
		
		cout << "Input the PMT number in Channel 1 \n" ;
		cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
		cin  >> channel[1]; 
		cout <<endl;
		
		cout << "Input the PMT VOLTAGE in Channel 1 \n" ;
		cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
		cin  >> Gain[1]; 
		cout <<endl;
		
		
		cout << "Input the PMT number in Channel 2 \n" ;
		cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
		cin  >> channel[2]; 
		cout <<endl;
		
		cout << "Input the PMT VOLTAGE in Channel 2 \n" ;
		cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
		cin  >> Gain[2]; 
		cout <<endl;
	
		cout << "Input the PMT number in Channel 3 \n";
		cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
		cin  >> channel[3]; 
		cout <<endl;
		
		cout << "Input the PMT VOLTAGE in Channel 3 \n" ;
		cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
		cin  >> Gain[3]; 
		cout <<endl;
		
		
		cout <<"Please verifiy the following: "<<endl;
		for (int i=0; i<4; i++){
			if (channel[i]<10)
				sprintf(histname,"NB000%d is in Channel %d Biased at %d Volts \n",channel[i], i, Gain[i]);
			if (channel[i]>=10 && channel[i] <100)
				sprintf(histname,"NB00%d is in Channel %d Biased at %d Volts \n",channel[i],  i, Gain[i]);
			if (channel[i]>=100)
				sprintf(histname,"NB0%d is in Channel %d  Biased at %d Volts \n",channel[i],  i, Gain[i]);
			cout << histname ;
		}
		
		cout <<"Is this correct? (y/n)  "<<endl;
		cin>>answer;
		cout <<answer<<endl;
		
	}
	
	
	
	
	// ***************
	// * Set up ROOT *
	// ***************
	// Create default ROOT application.
	TApplication *ta=new TApplication("ta",&argc,argv);
	
	// Data histogram.
	//TH1D* sData=newTH1D("data","Single Photon Energy;Channel;Counts",2000,-1000,9000);
	int PMT[4] = {0,0,0,0};
	
	
	vector<double> centroid(4,0), centroid_error(4,0),area(4,0), area_error(4,0);
	char filename[30];
	// Fitting the SPE Spectrum =======================================================================
	for (int r=0;r<4;r++){
			
		
	
		// Create canvas, allowing for window close.
	
	
		// *************************
		// * Create output spectra *
		// *************************
		TCanvas *tc=new TCanvas("Canvas","ROOT Canvas",1);
		tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
		tc->SetGrid();

		int mod =r;

		if (channel[r]<10)		
			sprintf(filename, "SPE/PMT_NB000%d_HV%d_AfterPulse.root",channel[r],  Gain[r]);
		if (channel[r]>=10 && channel[r] <100)
			sprintf(filename, "SPE/PMT_NB00%d_HV%d_AfterPulse.root",channel[r], Gain[r]);
		if (channel[r]>=100)
			sprintf(filename, "SPE/PMT_NB0%d_HV%d_AfterPulse.root",channel[r],  Gain[r]);
			
		TFile s(filename);
		s.ls();
		
		char root_name[30];
		sprintf(root_name, "TIME%d;1.root",mod);
		TH1D *speData = (TH1D*)s.Get(root_name);
		
		speData->GetYaxis()->SetTitle("Counts ");
		speData->GetYaxis()->SetTitleOffset(1.5);
		speData->GetXaxis()->SetTitle("Time [usec]");	
		
		// ****************
		// * Analyze data *
		// ****************
	
		// Save spectra.
		//ns().SaveData("rawspec.root");
	
		// Draw spectrum.
		speData->Draw();
		
		//Finding the location of the peaks
		TSpectrum *ss = new TSpectrum(6,3);
  		Int_t nfound = ss->Search(speData,6,"goff",0.0002);
  		
  		Double_t *peaks;
  		Double_t *intensity;
		peaks = ss->GetPositionX();
		intensity = ss->GetPositionY();
		//double peaks[5]={0.,0.,0.,0.,0.);
		//double 
		double maxI=0;
		int trigger=0;
   		for (int i =0; i<6;i++){
   			if (intensity[i]> maxI){
   				maxI = intensity[i];
   				trigger = i;
   			}
   		}
   		
   		
   		if (channel[r]<10)
			sprintf(histname,"PMT NB000%d: \n",channel[r]);
		if (channel[r]>=10 && channel[r] <100)
			sprintf(histname,"PMT NB00%d: \n",channel[r]);
		if (channel[r]>=100)
			sprintf(histname,"PMT NB0%d: \n",channel[r]);
		cout << histname ;
		for (int i=0; i<6; i++){
			if (i != trigger && peaks[i]>0.0){
				double ratio= intensity[i]/maxI;
				//for (int i=0; i<4; i++){
					
				//}
				printf("After-pulsing time of %f with ratio of %f \n", peaks[i],ratio);
			}
		
		}
		// Enter run loop.
		ta->Run("false");
	}
	// Making the HV fit ========================================================================
	
	
	return 0;
}





