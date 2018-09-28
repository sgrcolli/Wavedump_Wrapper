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






int main(int argc, char **argv)
{
	randomSeedTime();
	// Create default ROOT application.
		TApplication *ta=new TApplication("ta",&argc,argv);
		
	TH1D* Wave = new TH1D("Wave","Waveform; Time (nsec); Charge (mV-ns)",25005,-20000,30000);	
	TH1D* Search = new TH1D("Search", "Search Window; Time (ns); Charge (mV-ns)", 200, 0., 200.);
	TH2D* SPE2D = new TH2D("SPE2D", "Search Window; Time (us); Charge (mV-ns)", 1000, -20, 30.,1000,0,10000);
	for(int w=0;w<10001;w++){
	//for(int w=1;w<100;w++){
		//Stores one waveform for processing
		printf("\n\n\n Waveform %d \n\n\n",w);
		char filename[200]= "";
		if (w<10)
			sprintf(filename, "../../Data2/Data2_0000%d.csv",w);
		if (w>=10 && w<100)
			sprintf(filename, "../../Data2/Data2_000%d.csv",w);
		if (w>=100 && w<1000)
			sprintf(filename, "../../Data2/Data2_00%d.csv",w);
		if (w>=1000 && w<10000)
			sprintf(filename, "../../Data2/Data2_0%d.csv",w);
		if (w==10000)
			sprintf(filename, "../../Data2/Data2_%d.csv",w);	
		ifstream file(filename);
		
		string dummy;
	
		getline(file,dummy,',');getline(file,dummy,',');getline(file,dummy);
		for(int i=0;i<25007;i++){
			getline(file,dummy,',');
			double time = atof(dummy.c_str());
			//printf("Time: %f   \n ", time);
		
			getline(file,dummy,',');
			double value = atof(dummy.c_str());
			//printf("Channel: %f    ", value);
		
			getline(file,dummy);
			double channel = -1*(atof(dummy.c_str()));
			if (channel>2000)
				channel=1000;
			//printf("Channel: %f    \n", channel);
		
		
		  // printf("Time: %f   \n ", time);
		   Wave->SetBinContent(i,channel);
		}
		//closes the file
		file.close();
		
		
		//Rolling Window
		for(int i=30;i<24070;i++){
			double background = 0;
			double pulse =0;
			
			
			int maxBin = 0;
			double maxBack =0;
			
			for(int j=0; j<30; j++){
				
				background += Wave->GetBinContent(i+j-30);	
				pulse += Wave->GetBinContent(i+j);
				background += Wave->GetBinContent(i+j+30);	
				
				//Also finding the max value in the window.
				if (maxBack < Wave->GetBinContent(i+j)){
				
					maxBack = Wave->GetBinContent(i+j);
					maxBin = i+j;
				}
			}
			double ChargeDiff = (pulse-background*.5);
		
		
			if (ChargeDiff>150.00){
				
				for (int s=-100;s<100; s++){
					Search->SetBinContent(s+1, Wave->GetBinContent(maxBin+s+1));
				}	
				
				int binmax = Search->GetMaximumBin()+maxBin; 
				printf("Max Bin %d \n",binmax);
			
				
				int gates[8] ={binmax-30,binmax-20,binmax-10,binmax,binmax+10,binmax+20,binmax+30,binmax+40};
				
				if(binmax+80<25005){
					//Define the accumulators
					double A0=0;double A1=0;double A2=0;double A3=0;double A4=0;double A5=0;double A6=0;
		
					//Shift the search window past the gates
					i = binmax+300;
				
			
					for (int b=1; b<25005; b++){
			
						int time = b;
						if (time>=gates[0] && time<gates[1]){
							
							A0+=Wave->GetBinContent(b);
						}	
						if (time>=gates[1] && time<gates[2]){
							
							A1+=Wave->GetBinContent(b);
						}
						if (time>=gates[2] && time<gates[3]){
							 
							A2+=Wave->GetBinContent(b);
						}
						if (time>=gates[3] && time<gates[4]){
							
							A3+=Wave->GetBinContent(b);
						}
						if (time>=gates[4] && time<gates[5]){
							
							A4+=Wave->GetBinContent(b);
						}
						if (time>=gates[5] && time<gates[6]){
							
							A5+=Wave->GetBinContent(b);
						}
						if (time>=gates[6] && time<gates[7]){
							
							A6+=Wave->GetBinContent(b); 
						}
					}
					double charge =A2+A3+A4-(A0+A1+A5+A6)*3.0/4.0;
					double ap_time = Wave->GetXaxis()->GetBinCenter(binmax);
					//printf("Charge %f \n",charge);
					SPE2D->Fill(ap_time/1000.,charge);
					//SPE1->Fill(A2+A3+A4-(A0+A1)*3.0/2.0);
					//SPE2->Fill(A2+A3+A4-(A0+A1+A5+A6)*3.0/4.0);
		
				}
			}

		
		}
		
		//================= Reads in the headers and assigns values for things=============
	

	

	


	
	}
		//Create canvas allowing for window close
		TCanvas *tc0= new TCanvas("Canvas0","ROOT Canvas",1);
			tc0->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
			tc0->SetGrid();
	
		SPE2D->Draw("Same");
		SPE2D->SaveAs("spe2d.root");
		ta->Run();
	return 0;
}

