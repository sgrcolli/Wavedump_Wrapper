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
	
	//Stores one waveform for processing
	TH1D* Wave = new TH1D("Wave","Waveform; Time (ns); ADC Counts",1024,0,204.8);
	
	
	//Single Photoelectron Spectra with averaged accumulators
	TH1D **SPE=new TH1D*[4];	
	for (int w=0;w<4;w++){
		sprintf(histname, "SPE%d",w);
		SPE[w] = new TH1D(histname,"Single Photo-Electron; Charge (mV-ns); Counts",1500.0,-500.0,2000.0);
	}

	int totalwaves[4]={0,0,0,0};

       
	//================= Reads in the headers and assigns values for things=============
	

	//================= Reads in waveforms of length 1024 ==================

	//Include a counter to know the code is still running
	//cout << "got here " << endl;
	int counter = 0;
	for (int w=0; w<4; w++){

		char filename[200]= "";
		sprintf(filename,"../../Data/wave_%d.dat",w);
		ifstream fin(filename);
		
		for (int i=0; i<6; i++){
			//Read in the header for the script
			int header=0.;
			fin.read((char*)&header,sizeof(int));

		}
		counter = 0;
		while (fin.is_open()&&fin.good()&&!fin.eof()){
			counter++;
			//Are we there yet??
			if (counter%10000==0)
				printf("Waveform Progress: %d \n", counter);
	 
			//Records and ind. waveform into
			for (int i=0; i<1030; i++){
				//Read in result.
				float result=0.;
				fin.read((char*)&result,sizeof(float));
			
				if (i<1024){
					//inact an arbitrary offset in the data to make the peak
					double aoff = 2700;
					double flip_signal = (result-aoff)*-1.0;
					Wave->SetBinContent(i+1,flip_signal);
				
				}
			}
	
		
			//Determine the location of the peak
			int binmax = Wave->GetMaximumBin(); 
		     double maxtime = Wave->GetXaxis()->GetBinCenter(binmax);
			//printf("maxtime: %f\n",maxtime);
		
			int gates[8] ={binmax-300,binmax-200,binmax-100,binmax,binmax+100,binmax+200,binmax+300,binmax+400};
			


			//Peak must appear in reasonable location relative to the trigger
			if (maxtime>60.0 && maxtime<124.8){
				//Define the accumulators
				double A0=0;double A1=0;double A2=0;double A3=0;double A4=0;double A5=0;double A6=0;
				for (int i=1; i<=1024; i++){
				
					int time = i;
					if (time>=gates[0] && time<=gates[1]){
					
						A0+=Wave->GetBinContent(i);
					}	
					if (time>=gates[1] && time<=gates[2]){

						A1+=Wave->GetBinContent(i);
					}
					if (time>=gates[2] && time<=gates[3]){
					
						A2+=Wave->GetBinContent(i);
					}
					if (time>=gates[3] && time<=gates[4]){
					
						A3+=Wave->GetBinContent(i);
					}
					if (time>=gates[4] && time<=gates[5]){
					
						A4+=Wave->GetBinContent(i);
					}
					if (time>=gates[5] && time<=gates[6]){
						A5+=Wave->GetBinContent(i);
					}
					if (time>=gates[6] && time<=gates[7]){
				
						A6+=Wave->GetBinContent(i); 
					}

				}

			

				//Filling all the SPE
				double ADC_Counts = A2+A3+A4-(A0+A1+A5+A6)*3.0/4.0;
				double WaveCharge =  ADC_Counts*.2/4096.0*1.0e3;
				//printf("WaveCharge %f \n",WaveCharge);
				SPE[w]->Fill(WaveCharge);
			}
			
		}
		
		totalwaves[w]=counter;
		SPE[w]->Scale(1./(counter));
		//closes the wave-dump file
		fin.close();	
	}
	//Print out total number of waves for the relative quantum efficiency
	for (int i=0; i<4; i++)
		printf("Total Triggers from Wave %d: %d \n", i, totalwaves[i]);
	
	//After pulsing rate analysis==============================================================
	
	// Characterize the background
	
	//Characterize the triggers
	
	//Characterize the after pulses
	for (int i=0; i<4; i++)
		printf("Total Triggers from PMT %d: %d \n", Gain[i], totalwaves[i]);
	// Total Rate ==============================================================================
	//Create canvas allowing for window close
	TApplication *ta= new TApplication("ta",&argc,argv);
	TCanvas *tc1= new TCanvas("Canvas1","ROOT Canvas",1);
	tc1->Connect("TCanvas1","Closed()","TApplication",gApplication, "Terminate()");
	tc1->SetGrid();
	
	SPE[0]->Draw("Same");
	SPE[1]->Draw("Same");
	SPE[2]->Draw("Same");
	SPE[3]->Draw("Same");
	
	
	
	for (int i=0;i<4;i++){

		if (channel[i]<10)		
			sprintf(histname, "SPE_Root/PMT_NB000%d_HV%d_Analysis.root",channel[i],  Gain[i]);
		if (channel[i]>=10 && channel[i] <100)
			sprintf(histname, "SPE_Root/PMT_NB00%d_HV%d_Analysis.root",channel[i], Gain[i]);
		if (channel[i]>=100)
			sprintf(histname, "SPE_Root/PMT_NB0%d_HV%d_Analysis.root",channel[i],  Gain[i]);
			
		SPE[i]->SaveAs(histname);
	}
		
	
		
	ta->Run();
	
	
	return 0;
}

