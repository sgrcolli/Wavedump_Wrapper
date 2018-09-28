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

inline double linear(double x, double m, double b)
{
	return m*x+b;
}
inline double polynomial(double x, double p1, double p2, double p3, double p4, double p5)
{
	return p1*x*x*x*x+p2*x*x*x+p3*x*x+p4*x+p5;
}
double pfitt(double x, double *par){
	return polynomial(x, par[0],par[1],par[2],par[3],par[4]);
}
double fitfunc(double *x, double *par){
	return pfitt(x[0],par);
}

void addGraph(double (*func)(double,double*),int color ,double *p,double xmin,double xmax)
{
	//Extra graphs.
	int graphSteps=500;

	// Add error bands.
	TGraph *tg1=new TGraph(graphSteps);
	for (int i=0;i<graphSteps;i++)  {
		double x=xmin+(double)i*(xmax-xmin)/((double)graphSteps+1.);
		double y=func(x,p);
		tg1->SetPoint(i,x,y);
		
	}
	tg1->SetLineColor(color);
	tg1->SetLineWidth(2);
	tg1->Draw("SAME");
}



int main(int argc, char **argv)
{
	randomSeedTime();

	//Stores one waveform for processing
	TH1D* Wave = new TH1D("Wave","Waveform; Time (ns); ADC Counts",1024,0,204.8);
	
	//Single Photoelectron Spectra with averaged accumulators
	TH1D* SPE1 = new TH1D("SPE1","Single Photo-Electron; ADC Counts; Counts",1500.0,-500,40000); //pre-delay
	TH1D* SPE2 = new TH1D("SPE2","Single Photo-Electron; ADC Counts; Counts",1500.0,-500,40000); //post-delay
	TH1D* SPE3 = new TH1D("SPE3","Single Photo-Electron; ADC Counts; Counts",1500.0,-500,40000); //fitter spe
	
	ifstream fin("../Data/wave_1.dat");
	//================= Reads in the headers and assigns values for things=============
	for (int i=0; i<6; i++){
		//Read in the header for the script
		int header=0.;
		fin.read((char*)&header,sizeof(int));

	}

	//================= Reads in waveforms of length 1024 ==================

	
	while (fin.is_open()&&fin.good()&&!fin.eof()){
		//Create default ROOT application.
		TApplication *ta= new TApplication("ta",&argc,argv);
		//Create canvas allowing for window close
		TCanvas *tc= new TCanvas("Canvas","ROOT Canvas",1);
		tc->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
		tc->SetGrid();
		
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
		//Draw the wavform
		Wave->Draw("Same");

		//Determine the location of the peak
		int binmax = Wave->GetMaximumBin(); 
                double maxtime = Wave->GetXaxis()->GetBinCenter(binmax);
		//printf("maxtime: %f\n",maxtime);
		double delaygate = 20.0; //parameter for defining the delay in PSD
		//double gates[8] = {maxtime-60.,maxtime-40.,maxtime-20.,maxtime,maxtime+delaygate,maxtime+40.,maxtime+60.,maxtime+80.};
		int gates[8] ={binmax-300,binmax-200,binmax-100,binmax,binmax+100,binmax+200,binmax+300,binmax+400};
		//int gates[8] ={binmax-300,binmax-150,binmax-50,binmax,binmax+75,binmax+150,binmax+300,binmax+400};
		
                // Defining and plotting accumulators=========================================================================================================================
		//Create your accumulators, only for trouble shooting, comment out when doing analysis 
		TH1D* AC0 = new TH1D("Ac0","Accumulator 0; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC1 = new TH1D("Ac1","Accumulator 1; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC2 = new TH1D("Ac2","Accumulator 2; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC3 = new TH1D("Ac3","Accumulator 3; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC4 = new TH1D("Ac4","Accumulator 4; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC5 = new TH1D("Ac5","Accumulator 5; Time (ns); ADC Counts",1024,0.,204.8);
		TH1D* AC6 = new TH1D("Ac6","Accumulator 6; Time (ns); ADC Counts",1024,0.,204.8);
		//Only for background accumulators
		TH1D* ACF = new TH1D("AcF","Accumulator Fit; Time (ns); ADC Counts",1024,0.,204.8);
		
		//Define the accumulators
		double A0=0;double A1=0;double A2=0;double A3=0;double A4=0;double A5=0;double A6=0;

		//Peak must appear in reasonable location relative to the trigger
		if (maxtime>60.0 && maxtime<124.8){
			
			for (int i=1; i<1025; i++){
				
				int time = i;
				if (time>=gates[0] && time<gates[1]){
					AC0->SetBinContent(i,Wave->GetBinContent(i)); 
					ACF->SetBinContent(i,Wave->GetBinContent(i));
					A0+=Wave->GetBinContent(i);
				}	
				if (time>=gates[1] && time<gates[2]){
					AC1->SetBinContent(i,Wave->GetBinContent(i));
					ACF->SetBinContent(i,Wave->GetBinContent(i));
					A1+=Wave->GetBinContent(i);
				}
				if (time>=gates[2] && time<gates[3]){
					AC2->SetBinContent(i,Wave->GetBinContent(i)); 
					A2+=Wave->GetBinContent(i);
				}
				if (time>=gates[3] && time<gates[4]){
					AC3->SetBinContent(i,Wave->GetBinContent(i));
					A3+=Wave->GetBinContent(i);
				}
				if (time>=gates[4] && time<gates[5]){
					AC4->SetBinContent(i,Wave->GetBinContent(i));
					A4+=Wave->GetBinContent(i);
				}
				if (time>=gates[5] && time<gates[6]){
					AC5->SetBinContent(i,Wave->GetBinContent(i)); 
					ACF->SetBinContent(i,Wave->GetBinContent(i));
					A5+=Wave->GetBinContent(i);
				}
				if (time>=gates[6] && time<gates[7]){
					AC6->SetBinContent(i,Wave->GetBinContent(i));
					ACF->SetBinContent(i,Wave->GetBinContent(i));
					A6+=Wave->GetBinContent(i); 
				}

			}
			
			SPE1->Fill(A2+A3+A4-(A0+A1)*3.0/2.0);
			SPE2->Fill(A2+A3+A4-(A0+A1+A5+A6)*3.0/4.0);
		}
		
		
		//================================================ Fitting part =========================================================================================
		if (maxtime>60.0 && maxtime<124.8){		
			//Create canvas allowing for window close
			TCanvas *tc3= new TCanvas("Canvas3","ROOT Canvas",1);
			tc3->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
			tc3->SetGrid();
		
			double fitXMin = maxtime-60.0; 
			double fitXMax = maxtime+80.0;
		
			//Define the fit and parameter
			TF1 *tf11=new TF1("data1Fit",fitfunc,fitXMin,fitXMax,5);
		
			tf11->SetParameter(0,0);
			tf11->SetParameter(1,0);
			tf11->SetParameter(2,0);
			tf11->SetParameter(3,0);
			tf11->SetParameter(4,0);
			
			ACF->Draw("Same");
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

			//Perform fit
			TFitResultPtr tfrp1=ACF->Fit("data1Fit","RSE");

			// Print results. 
			tfrp1->Print();
			
			//Summing area under the curve
			double AF= 0;
			//Determing the SPE using the fits
			for (int i=gates[2]; i< gates[5];i++){
				AF+=pfitt(ACF->GetXaxis()->GetBinCenter(i),tf11->GetParameters());
			}
			SPE3->Fill(A2+A3+A4-AF);
		}
		
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

	//Create canvas allowing for window close
	TApplication *ta= new TApplication("ta",&argc,argv);
	TCanvas *tc1= new TCanvas("Canvas1","ROOT Canvas",1);
	tc1->Connect("TCanvas1","Closed()","TApplication",gApplication, "Terminate()");
	tc1->SetGrid();
	
	SPE1->Draw("Same");
	

	
	TCanvas *tc2= new TCanvas("Canvas2","ROOT Canvas",2);
	tc2->Connect("TCanvas2","Closed()","TApplication",gApplication, "Terminate()");
	tc2->SetGrid();
	
	
	SPE2->Draw("Same");
	
	TCanvas *tc3= new TCanvas("Canvas3","ROOT Canvas",3);
	tc3->Connect("TCanvas3","Closed()","TApplication",gApplication, "Terminate()");
	tc3->SetGrid();
	
	
	SPE3->Draw("Same");
	
	SPE3->SaveAs("spe3.root");
	SPE2->SaveAs("spe2.root");
	SPE1->SaveAs("spe1.root");

	ta->Run();
	//closes the wave-dump file
	fin.close();
	
	return 0;
}

