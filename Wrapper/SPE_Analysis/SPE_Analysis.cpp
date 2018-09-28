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

//Gain Curve Fitters ==========================================================================================================================
inline double PowerFunc(double x, double k, double n)
{
	return pow((k*x),n)*10.00;
}
double fitPow(double *x, double *k)
{
	return PowerFunc(x[0], k[0], k[1]);
}
// A Bunch of fitting stuff for the single SPE =================================================================================================


inline double gauss(double x,double sigma)
{
	// Normalized gaussian.
	return exp(-x*x/(2*sigma*sigma))/(sqrt(2.*kMathPi)*sigma);
      	//return exp(-x*x/(2*sigma*sigma))/(sqrt(2.*3.1415926)*sigma);

}


inline double expBackground(double x,double *p)
{
	return sqrt(p[0]*p[0])+sqrt(p[1]*p[1])*exp(-x*p[2]);
}

inline double photon1Signal(double x,double *p)
{
	return abs(p[3]*gauss(x-p[4],p[5]));
}

inline double photon2Signal(double x,double *p)
{
	return abs(p[3]*p[6]*gauss(x-2*p[4],sqrt(2.)*p[5]));
}


double fitFunc(double *x,double *p)
{
	return expBackground(x[0],p)+photon1Signal(x[0],p)+photon2Signal(x[0],p);
}

void addGraph(double (*func)(double,double*),int color,double *p,double xmin,double xmax)
{
	// Extra graphs.
	int graphSteps=500;
	
	// Add error bands.
	TGraph *tg1=new TGraph(graphSteps);
	for (int i=0;i<graphSteps;i++)	{
		double x=xmin+(double)i*(xmax-xmin)/((double)graphSteps+1.);
		double y=func(x,p);
		tg1->SetPoint(i,x,y);
	}
	tg1->SetLineColor(color);
	tg1->SetLineWidth(2);
	tg1->Draw("SAME");
}
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
			sprintf(filename, "SPE_Root/PMT_NB000%d_HV%d_Analysis.root",channel[r],  Gain[r]);
		if (channel[r]>=10 && channel[r] <100)
			sprintf(filename, "SPE_Root/PMT_NB00%d_HV%d_Analysis.root",channel[r], Gain[r]);
		if (channel[r]>=100)
			sprintf(filename, "SPE_Root/PMT_NB0%d_HV%d_Analysis.root",channel[r],  Gain[r]);
			
		TFile s(filename);
		s.ls();
		
		char root_name[30];
		sprintf(root_name, "SPE%d;1.root",mod);
		TH1D *speData = (TH1D*)s.Get(root_name);
		
		speData->GetYaxis()->SetTitle("Counts ");
		speData->GetYaxis()->SetTitleOffset(1.5);
		speData->GetXaxis()->SetTitle("Charge (mv*ns)");	
		
		// ****************
		// * Analyze data *
		// ****************
	
		// Save spectra.
		//ns().SaveData("rawspec.root");
	
		// Draw spectrum.
		speData->Draw();
	
		
		// Fit the data.
		//double fitXMin=75.;
		//found 75 to be too high, 30 was more effective
		double fitXMin=30.0; //Ron Collins, 19/09/18
		
		double fitXMax=2000.;
	
		
		TF1 *tf11=new TF1("data1Fit",fitFunc,fitXMin,fitXMax,7);
		tf11->SetLineWidth(3);

		// // Exponential background.
		// tf11->SetParameter(0,4.77335);
		// tf11->SetParLimits(0,0,1000);
		// tf11->SetParameter(1,91.6516);		
		// tf11->SetParameter(2, 0.0000656504);
		// //tf11->SetParameter(2, 2.0);
		// //tf11->SetParLimits(2,2,1000);

		// // Single photon gaussian.
		// tf11->SetParameter(3,195784);
		// tf11->SetParameter(4,158.745);
		// tf11->SetParameter(4,350.745);
		// tf11->SetParameter(5,161.693);


		// // Double and triple photon gaussian.
		// tf11->SetParameter(6, 0.223033);

		// Exponential background.
		tf11->SetParameter(0,2.45E-6);
		//tf11->SetParLimits(0,0,1000);

		tf11->SetParameter(1,3.6E-3);		
		//tf11->SetParLimits(1, 0, 0.0005);

		tf11->SetParameter(2, 0.023);
		//tf11->SetParLimits(2, 0,0.025);

		// Single photon gaussian.
		tf11->SetParameter(3,-0.416);
		tf11->SetParameter(4,321);
		tf11->SetParameter(5,99.2);


		// Double and triple photon gaussian.
		tf11->SetParameter(6, -0.19);

	
		// Perform fit.
		TFitResultPtr tfrp1=speData->Fit("data1Fit","RSE");
	
		// Print results.
		tfrp1->Print();
		
		area[r] =(abs(tf11->GetParameter(3))+2.0*abs(tf11->GetParameter(3)*tf11->GetParameter(6)))/1.67;
		area_error[r]=sqrt((tf11->GetParError(3))*(tf11->GetParError(3)));

		centroid[r] =tf11->GetParameter(4);
		centroid_error[r]=sqrt((tf11->GetParError(4))*(tf11->GetParError(4))+centroid[r]*centroid[r]/area[r])+centroid[r]*.03;
					
		
		
		//Area of the curve
		printf("Area: %f \n",area[r]);
		//cout<< "area error" << AreaError[r] << endl;
		printf("Centroid: %f \n",centroid[r]);
	
		//addGraph(skewedBackground,kGreen+2,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(expBackground,kMagenta+2,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(photon1Signal,kGreen+2,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(photon2Signal,kCyan+2,tf11->GetParameters(),fitXMin,fitXMax);
		//addGraph(photon3Signal,kOrange+3,tf11->GetParameters(),fitXMin,fitXMax);
	
		// Update canvas.
		tc->Update();
		tc->Paint();
		tc->Draw();
		tc->Modified();
		
		//Local Minima and Maxima
		int binmax = speData->GetMaximumBin(); 
		double minVal = speData->GetBinContent(binmax);
		int minBin =0;
		int pbin = speData->FindBin(350);
		for (int i=binmax;i<=pbin;i++){	
				double scope = speData->GetBinContent(i);
				if(minVal>scope) {
						minVal = scope;
						minBin = i;
						//printf("min: %f \n",minVal);
				}
			}  
		double maxVal = speData->GetBinContent(minBin);
		for (int i=minBin;i<=1500;i++){	
				double scope = speData->GetBinContent(i);
				if(maxVal<scope) {
						maxVal = scope;
						//printf("max: %f \n",maxVal);
				}
			}  
			
   		double PV =maxVal/minVal;
   		printf("Peak to Valley: %f \n",PV);
		// Enter run loop.
		ta->Run("false");
	}
	// Making the HV fit ========================================================================
	
	
	return 0;
}





