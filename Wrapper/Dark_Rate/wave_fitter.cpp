#include "ns.h"
//Standard library include files.
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TObject.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraph.h"

using namespace std;





int main(int argc, char **argv)
{
	randomSeedTime();
		
	//Create default ROOT application.
	TApplication *ta= new TApplication("ta",&argc,argv);
	
	//Stores one waveform for processing
	TH1D* Wave = new TH1D("Wave","Waveform; Time (ns); ADC Counts",1024,0,1024);
	
	//Single Photoelectron Spectra with averaged accumulators
	//TH1D* SPE1 = new TH1D("SPE1","Single Photo-Electron; ADC Counts; Counts",1500.0,-500,40000); //pre-delay
	
	
	
	ifstream fin("../../Data/wave_0.dat");
	//================= Reads in the headers and assigns values for things=============
	for (int i=0; i<6; i++){
		//Read in the header for the script
		int header=0.;
		fin.read((char*)&header,sizeof(int));

	}

	//================= Reads in waveforms of length 1024 ==================

	int counter = 0;
	while (fin.is_open()&&fin.good()&&!fin.eof()){
		counter ++;
		
		
		//Records and ind. waveform into
		for (int i=0; i<1030; i++){
			//Read in result.
			float result=0.;
			fin.read((char*)&result,sizeof(float));
			if(!fin.eof())
				//printf("P: %f\n",result);
			if (i<1024){
				//inact an arbitrary offset in the data to make the peak
				double aoff = 2700;
				double flip_signal = (result-aoff)*-1.0;
				Wave->SetBinContent(i+1,flip_signal);
			}
		}
	
		
		//Set up rolling windows of 60 usec =======================================================
		
		
		
		
		
		// create loops for this background
		for (int i=60; i<=1024; i++){
			
			double background = 0;
			double pulse = 0;
			
			double maxBack= 0;
			int maxBin=0;
			for (int j=1; j<=60; j++){
				
				background += Wave->GetBinContent(i+j-60);	
				pulse += Wave->GetBinContent(i+j);
				
				//Also finding the max value in the window.
				if (maxBack < Wave->GetBinContent(i+j-60)){
				
					maxBack = Wave->GetBinContent(i+j-60);
					maxBin = i+j-60;
				}
			}
			
			
			//Determing if the charge exceeds 100 mV-ns
			double ChargeDiff = (pulse-background)/4096.0*1.0e3;
			//printf ("Charge Difference: %f \n",ChargeDiff);
			
			if (ChargeDiff>50.00){
				TH1D* Search = new TH1D("Search", "Search Window; Time (ns); ADC Counts", 200, 0., 200.);
				for (int s=0;s<200; s++){
					Search->SetBinContent(s+1, Wave->GetBinContent(maxBin+s+1));
				}	
				
				int binmax = Search->GetMaximumBin()+maxBin; 
				
				//Create your accumulators, only for trouble shooting, comment out when doing analysis 
				TH1D* AC0 = new TH1D("Ac0","Accumulator 0; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC1 = new TH1D("Ac1","Accumulator 1; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC2 = new TH1D("Ac2","Accumulator 2; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC3 = new TH1D("Ac3","Accumulator 3; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC4 = new TH1D("Ac4","Accumulator 4; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC5 = new TH1D("Ac5","Accumulator 5; Time (ns); ADC Counts",1024,0.,1024.);
				TH1D* AC6 = new TH1D("Ac6","Accumulator 6; Time (ns); ADC Counts",1024,0.,1024.);
				
				int gates[8] ={binmax-60,binmax-40,binmax-20,binmax,binmax+20,binmax+40,binmax+60,binmax+80};
				
				if(binmax+80<1024){
					//Define the accumulators
					double A0=0;double A1=0;double A2=0;double A3=0;double A4=0;double A5=0;double A6=0;
		
					//Shift the search window past the gates
					i = binmax+80;
				
			
					for (int b=1; b<1025; b++){
			
						int time = b;
						if (time>=gates[0] && time<gates[1]){
							AC0->SetBinContent(b,Wave->GetBinContent(b)); 
							A0+=Wave->GetBinContent(b);
						}	
						if (time>=gates[1] && time<gates[2]){
							AC1->SetBinContent(b,Wave->GetBinContent(b));
							A1+=Wave->GetBinContent(b);
						}
						if (time>=gates[2] && time<gates[3]){
							AC2->SetBinContent(b,Wave->GetBinContent(b)); 
							A2+=Wave->GetBinContent(b);
						}
						if (time>=gates[3] && time<gates[4]){
							AC3->SetBinContent(b,Wave->GetBinContent(b));
							A3+=Wave->GetBinContent(b);
						}
						if (time>=gates[4] && time<gates[5]){
							AC4->SetBinContent(b,Wave->GetBinContent(b));
							A4+=Wave->GetBinContent(b);
						}
						if (time>=gates[5] && time<gates[6]){
							AC5->SetBinContent(b,Wave->GetBinContent(b)); 
							A5+=Wave->GetBinContent(b);
						}
						if (time>=gates[6] && time<gates[7]){
							AC6->SetBinContent(b,Wave->GetBinContent(b));
							A6+=Wave->GetBinContent(b); 
						}
					}

				
		
					//SPE1->Fill(A2+A3+A4-(A0+A1)*3.0/2.0);
					//SPE2->Fill(A2+A3+A4-(A0+A1+A5+A6)*3.0/4.0);
		
					//=========================================== Plot Everything =============================================================================================
					//Create canvas allowing for window close
					TCanvas *tc0= new TCanvas("Canvas0","ROOT Canvas",1);
					tc0->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
					tc0->SetGrid();
		
					Wave->Draw("Same");

					AC0->Draw("Same");
					AC0->SetFillColor(kViolet+2);
					AC0->SetFillStyle(3002);

					AC1->Draw("Same");
					AC1->SetFillColor(kMagenta+2);
					AC1->SetFillStyle(3002);

					AC2->Draw("Same");
					AC2->SetFillColor(kBlue+2);
					AC2->SetFillStyle(3002);

					AC3->Draw("Same");
					AC3->SetFillColor(kCyan+2);
					AC3->SetFillStyle(3002);

					AC4->Draw("Same");
					AC4->SetFillColor(kGreen+2);
					AC4->SetFillStyle(3002);

					AC5->Draw("Same");
					AC5->SetFillColor(kOrange+2);
					AC5->SetFillStyle(3002);

					AC6->Draw("Same");
					AC6->SetFillColor(kYellow+2);
					AC6->SetFillStyle(3002);
		
					ta->Run("False");
				}
			}

		
		}
		
		
		

		
		
		
		
		
	}

	//Create canvas allowing for window close
	
	TCanvas *tc1= new TCanvas("Canvas1","ROOT Canvas",1);
	tc1->Connect("TCanvas1","Closed()","TApplication",gApplication, "Terminate()");
	tc1->SetGrid();
	
	//SPE1->Draw("Same");
	

	
	//SPE1->SaveAs("spe1.root");

	ta->Run();
	//closes the wave-dump file
	fin.close();
	
	return 0;
}

