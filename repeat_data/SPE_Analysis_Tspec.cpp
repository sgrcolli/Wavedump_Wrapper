// Project include files.
#include "ns.h"

// Standard library include files.
#include <cmath>
#include <cstdlib>
#include <iostream>

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

//Taken from the fitLED-old2.C file 
// Ron was here
//#include "RooRealVar.h"
//#include "RooAddPdf.h"
//#include "RooGaussian.h"
//#include "RooDataHist.h"
//#include "RooFitResult.h"
//#include "RooPlot.h"
//#include "RooDerivative.h"
//#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include <tuple> //unsure if this is right
//Ron stopped being here

//using namespace RooFit; //unsure if this is right, Ron was here

typedef std::tuple<double,double,double,double,double> InitParams;//Ron was here

InitParams initializeFit(TH1D* h);
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

inline double ChoozExpBackground(double x,double *p)
{
  return (p[8]/p[10]) *exp(- (x-p[9])/p[10] );
}

inline double pedestalFit(double x,double *p)
{
  return (p[0]/p[1])*sqrt(2.*kMathPi)*exp(-(pow(x-p[2],2))/(2*pow(p[1],2)));
}

inline double SPE_Fit(double x,double *p)
{
  return (p[3]/p[4])*sqrt(2.*kMathPi)*exp(-(pow(x-p[5],2))/(2*pow(p[4],2)));
}

inline double second_Photon_Fit(double x,double *p)
{
  //return  p[6]*(p[3]/(sqrt(.2)*p[4]))*sqrt(2.*kMathPi)*exp(-(pow(x-p[5],2))/(2*pow((sqrt(2.)*p[4]),2)));
  //return p[7] * (  ( (pow (p[4], 2) )/2 ) * exp(-p[6]) * (1/(sqrt(4. * kMathPi) * p[4] )) * exp( - (pow(x - 2., 2))/(4. * pow(p[4], 2) ) )   );
  return (p[3]/p[6])*sqrt(2.*kMathPi)*exp(-(pow(x-p[7],2))/(2*pow(p[6],2)));

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
  //cout << "x0 is " << x[0] << endl;
  return expBackground(x[0],p)+photon1Signal(x[0],p)+photon2Signal(x[0],p);
}

double fitFuncTSPEC(double *x,double *p)
{
  if (x[0] >  p[9]){
    return pedestalFit(x[0],p) + SPE_Fit(x[0],p) + second_Photon_Fit(x[0],p) + ChoozExpBackground(x[0],p);
  }
  else return pedestalFit(x[0],p) + SPE_Fit(x[0],p) + second_Photon_Fit(x[0],p);
    
  // double i = x[0];
  // double N_Ped = p[0];
  // double sigma_Ped = p[1];
  // double mean_Ped = p[2];
  // double N_SPE = p[3];
  // double simga_SPE = p[4];
  // double mean_SPE = p[5];
  // double dev2PE = p[6]
  // double mean2PE = p[7]
  // double Nexp = p[8];
  // double entryMin = p[9];
  // double tau = p [10];
  
  // double histogramMean = p[6];
  // double histogramEnteries = p[7] ;

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

  //std::cout << "I have " << argc << " arguments" << std::endl;
  std::cout << "first argument is " << argv[1] << std::endl;
  std::cout << "second argument is " << argv[2] << std::endl;


	// ******************
	// * Initialization *
	// ******************
	
	// Set up random number generator.
	randomSeedTime();
	
			
	int channel[5]={0,0,0,0,0};
	char answer;
	char histname[200]= "";
	int test;
	int Gain[5]={0,0,0,0,0};
	
	channel[0] = atoi(argv[1]);
	Gain[0] = atoi(argv[2]);

	channel[1] = atoi(argv[3]);
	Gain[1] = atoi(argv[4]);

	channel[2] = atoi(argv[5]);
	Gain[2] = atoi(argv[6]);

	channel[3] = atoi(argv[7]);
	Gain[3] = atoi(argv[8]);

	channel[4] = atoi(argv[9]);
	Gain[4] = atoi(argv[10]);


	//Read in the HV data ====================================================================================
	string hvfile = "../../HVScan.txt";
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
	//cout << "got here" << endl;
	//========================================================================================================
	
	// while(answer!='Y'&& answer!='y'){
		
	// 	//Determing the PMT number and the applied Voltage=====================================================
	// 	cout << "Input the PMT number in Channel 0 \n" ;
	// 	cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
	// 	cin  >> channel[0]; 
	// 	cout <<endl;
		
	// 	cout << "Input the PMT VOLTAGE in Channel 0 \n" ;
	// 	cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
	// 	cin  >> Gain[0]; 
	// 	cout <<endl;
		
	// 	cout << "Input the PMT number in Channel 1 \n" ;
	// 	cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
	// 	cin  >> channel[1]; 
	// 	cout <<endl;
		
	// 	cout << "Input the PMT VOLTAGE in Channel 1 \n" ;
	// 	cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
	// 	cin  >> Gain[1]; 
	// 	cout <<endl;
		
		
	// 	cout << "Input the PMT number in Channel 2 \n" ;
	// 	cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
	// 	cin  >> channel[2]; 
	// 	cout <<endl;
		
	// 	cout << "Input the PMT VOLTAGE in Channel 2 \n" ;
	// 	cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
	// 	cin  >> Gain[2]; 
	// 	cout <<endl;
	
	// 	cout << "Input the PMT number in Channel 3 \n";
	// 	cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
	// 	cin  >> channel[3]; 
	// 	cout <<endl;
		
	// 	cout << "Input the PMT VOLTAGE in Channel 3 \n" ;
	// 	cout << "ENTER THE VOLTAGE IN VOLTS \n" <<endl;
	// 	cin  >> Gain[3]; 
	// 	cout <<endl;
		
		
		// cout <<"Please verifiy the following: "<<endl;
		// for (int i=0; i<5; i++){
		// 	if (channel[i]<10)
		// 		sprintf(histname,"NB000%d is in Channel %d Biased at %d Volts \n",channel[i], i, Gain[i]);
		// 	if (channel[i]>=10 && channel[i] <100)
		// 		sprintf(histname,"NB00%d is in Channel %d Biased at %d Volts \n",channel[i],  i, Gain[i]);
		// 	if (channel[i]>=100)
		// 		sprintf(histname,"NB0%d is in Channel %d  Biased at %d Volts \n",channel[i],  i, Gain[i]);
		// 	cout << histname ;
		// }
		
		// cout <<"Is this correct? (y/n)  "<<endl;
		// cin>>answer;
		// cout <<answer<<endl;
		
		// }
	
	
	
	
	// ***************
	// * Set up ROOT *
	// ***************
	// Create default ROOT application.
	TApplication *ta=new TApplication("ta",&argc,argv);
	//cout << "got here 2" << endl;
	
	// Data histogram.
	//TH1D* sData=newTH1D("data","Single Photon Energy;Channel;Counts",2000,-1000,9000);
	int PMT[5] = {0,0,0,0,0};
	
	
	vector<double> centroid(5,0), centroid_error(5,0),area(5,0), area_error(5,0); //this might need to be 4,0 not sure
	char filename[30];
	// Fitting the SPE Spectrum =======================================================================
	for (int r=0;r<5;r++){
			
		
	
		// Create canvas, allowing for window close.
	
	
		// *************************
		// * Create output spectra *
		// *************************
		TCanvas *tc=new TCanvas("Canvas","ROOT Canvas",1);
		tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
		tc->SetGrid();

		//int mod =r;

		if (channel[r]<10)		
			sprintf(filename, "PMT_NB000%d_HV%d_Analysis.root",channel[r],  Gain[r]);
		if (channel[r]>=10 && channel[r] <100)
			sprintf(filename, "PMT_NB00%d_HV%d_Analysis.root",channel[r], Gain[r]);
		if (channel[r]>=100)
			sprintf(filename, "PMT_NB0%d_HV%d_Analysis.root",channel[r],  Gain[r]);
		//Ron was here	
		TFile s(filename);
		s.ls();
		TH1D *speData = 0;

		for (int mod = 0; mod < 4 ; mod++) {
		  char root_name[30];
		  sprintf(root_name, "SPE%d;1.root",mod);
		  //TH1D *speData = (TH1D*)s.Get(root_name);
		  speData = (TH1D*)s.Get(root_name);
		  if (speData == 0 ){
		    //cout << mod << endl;
		    //cout << "no SPE data!" << endl;
		  }
		  else break;
		}
		
		InitParams fitResults = initializeFit(speData);

		double PedestalPos =  std::get<0>(fitResults);
		double SPE_Ped     =  std::get<1>(fitResults); // SPE relative to the pedestal 
		double PedWidth    =  std::get<2>(fitResults);
		double SPE_Width   =  std::get<3>(fitResults);
		double SN_Ratio    =  std::get<4>(fitResults); 

		cout << endl;
		cout << "Pedestal position at "            <<  PedestalPos << endl;
		cout << "SPE relative to pedestal "        <<  SPE_Ped     << endl;
		cout << "width of pedestal is "            <<  PedWidth    << endl;
		cout << "width of SPE is "                 <<  SPE_Width   << endl;
		cout << "The ratio of signal to noise is " <<  SN_Ratio    << endl;
		//cout << endl;

		//Ron stopped being herer
		speData->GetYaxis()->SetTitle("Counts ") ;
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
	
		//fitXMin = PedestalPos + SPE_Ped + PedWidth; //Remember that the pedestal width is netative Ron09/10/2018
		fitXMin = -20;
		//Ron was here
		double meanSPE = PedestalPos + SPE_Ped;
		double meanPed = PedestalPos;
		double devSPE  = SPE_Width;
		double devPed  = -PedWidth;
		double N_SPE   = speData->GetBinContent(speData->FindBin(PedestalPos + SPE_Ped));
		double N_Ped   = speData->GetBinContent(speData->FindBin(PedestalPos));
		cout << "N SPE is " << N_SPE << endl;
		cout << "N pedestal is " <<  N_Ped << endl;
		double histogramMean = speData->GetMean();
		double histogramEntries = speData->GetEntries();		
		//double histogramEnteries = speData->Integral();
		double MinValDat = speData -> GetMinimum();

		double sec_PE_Mean  = 2.0*meanSPE + meanPed;
		double sec_PE_stdev = sqrt(2.0*devSPE*devSPE + devPed*devPed); 
		//double sec_PE_Mean_min = sqrt(2)*meanSPE;

		double FindMinPedPos  = speData->FindBin(PedestalPos) ;
		double FindMindSpePos = speData->FindBin(PedestalPos + SPE_Ped);
		double valMin = speData -> GetBinContent(FindMinPedPos);
		double entryMin = 0;
		double entryDat = MinValDat;
		
		cout << "FindMinPedPos is "  << FindMinPedPos << endl;
		cout << "FindMindSpePos is " << FindMindSpePos<< endl;

		for(uint entry = 0; entry < histogramEntries; entry++){
		  double CurrentBinWidth = speData -> GetBinWidth(entry);
		  entryDat = entryDat + CurrentBinWidth;
		  //cout << "histogram x value is " << entryDat << endl;
		  if(entry > FindMinPedPos && entry < FindMindSpePos){
		    double currentVal =  speData -> GetBinContent(entry); 
		    if (currentVal < valMin){
		      valMin = currentVal;
		      //entryDat = entry; 
		      entryMin = speData->GetBinCenter(entry);
		    }
		  }
		}

		//double binConv = FindMinPedPos/PedestalPos;
		double Nexp = valMin;
		double iMin = entryMin;
		double tau  = 1;
		cout << "Nexp is " << Nexp << endl;
		cout << "iMin is " << iMin << endl;
		cout << "tau is "  << tau  << endl;

		cout << endl;
		//Ron Stopped being here

		//TF1 *tf11=new TF1("data1Fit",fitFunc,fitXMin,fitXMax,7);
		TF1 *tf11=new TF1("data1Fit",fitFuncTSPEC, fitXMin,fitXMax,11); //Ron was here
		tf11->SetParameter(0, N_Ped);
		tf11->SetParameter(1, devPed);
		tf11->SetParameter(2, meanPed);
		tf11->SetParLimits(2, meanPed-devPed, meanPed+devPed);
		tf11->SetParameter(3, N_SPE);
		tf11->SetParameter(4, devSPE);
		tf11->SetParameter(5, meanSPE);
		tf11->SetParLimits(5, meanSPE-devSPE, meanSPE+devSPE); //if you want to force limits 
		tf11->SetParameter(6, sec_PE_stdev);
		tf11->SetParameter(7, sec_PE_Mean);
		tf11->SetParLimits(7, sec_PE_Mean-devSPE, sec_PE_Mean+devSPE ); //if you want to force limits 
		tf11->SetParameter(8,  Nexp);
		tf11->SetParameter(9,  iMin);
		tf11->SetParameter(10, tau);
		tf11->SetParLimits(10, 0, 1);

		//tf11->SetParameter(6, histogramMean);
		//tf11->SetParameter(7, histogramEnteries);
		//tf11->SetParameter(6, -0.416);//values taken from old code,-0.416 
		//tf11->SetParameter(7, -0.19); //values taken from old code,-0.19
		//tf11->SetParameter(8, 321);   //values taken from old code,321
		//tf11->SetParameter(9, 99.2);  //values taken from old code,99.2

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
		//tf11->SetParameter(0,2.45E-6); // worked on the old one
		//tf11->SetParLimits(0,0,1000);

		//tf11->SetParameter(1,3.6E-3); // worked on the old one		
		//tf11->SetParLimits(1, 0, 0.0005);

		//tf11->SetParameter(2, 0.023); // worked on the old one
		//tf11->SetParLimits(2, 0,0.025);

		// Single photon gaussian.
		//tf11->SetParameter(3,-0.416); // worked on the old one
		//tf11->SetParameter(4,321); // worked on the old one
		//tf11->SetParameter(5,99.2); // worked on the old one


		// Double and triple photon gaussian.
		//tf11->SetParameter(6, -0.19); // worked on the old one

	
		// Perform fit.
		TFitResultPtr tfrp1=speData->Fit("data1Fit","RSE");
	
		cout << "got here 3" << endl;

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

		//addGraph(expBackground,kMagenta+2,tf11->GetParameters(),fitXMin,fitXMax);
		//addGraph(photon1Signal,kGreen+2,tf11->GetParameters(),fitXMin,fitXMax);
		//addGraph(photon2Signal,kCyan+2,tf11->GetParameters(),fitXMin,fitXMax);

		addGraph(SPE_Fit,kCyan+2,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(pedestalFit,kGreen+2,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(second_Photon_Fit,kOrange+1,tf11->GetParameters(),fitXMin,fitXMax);
		addGraph(ChoozExpBackground,kMagenta+2,tf11->GetParameters(),iMin,fitXMax);
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


//Everything below this point was taken directly from the fiteLED-old2.C from Mathew needham
//Ron was here 
//#include <tuple>
//typedef std::tuple<double,double,double,double,double> InitParams;

//using namespace RooFit;

InitParams initializeFit(TH1D* h){

  TSpectrum *s = new TSpectrum(2,3);
  Int_t nfound = s->Search(h,2,"goff",0.0002);
  std::cout << "found peaks " << nfound << std::endl;
  
  /*** returns positions of Pedestal and Signal approx.****/
  Double_t *peaks;
  peaks = s->GetPositionX();
   std::cout << peaks[0] << " " << peaks[1]  << std::endl;
  
  /*** Find the valley ***/
  Int_t SigSig = h->FindBin(peaks[1]); //Signal bin approx **used to find valley position
  Int_t PedPed  = h->GetMaximumBin(); //pedestal bin
   int PedMax = h->GetMaximum(); // pedestal events

  int Compare[2] = {PedMax,PedMax}; // array to store 2 numbers for comparison to determine which is larger
  Int_t ValleyBin = 0; //to hold valley bin number


  /****Finding the Valley*****/
  for(int i = PedPed; i< SigSig ;i++){
    Compare[0] = h->GetBinContent(i);

    if(Compare[0] < Compare[1]){
      Compare[1] = Compare[0]; //finds minimum number of events
      ValleyBin = i; //sets bin number for minimum
    }//end if
  }//end for

 
  double ssignal = (peaks[1] -  h->GetBinCenter(ValleyBin))/2.5;
  double sped = (h->GetBinCenter(ValleyBin) - peaks[0])/2.5;
  
  // hack
  //  if (sped > 30) sped = 30;
  
   /*****Finding Pedestal events****/
  int Noise = 0; //to hold pedestal events
  for(int i = 0;i< ValleyBin;i++){
    Noise += h->GetBinContent(i); // total noise events
  }//end for

  /** Approximate fraction of photons per event**/
  int Events = h->GetEntries();  // # of events
  double Ratio = (double) Noise / (double) Events; // Ratio of Random events to photon events

  std::cout << "Peaks " << peaks[0] << " " << peaks[1] << " " << sped <<   " " << ssignal << " "   <<  Ratio << std::endl;

  //return InitParams(5,45, 2,5,0.4);
  return InitParams(peaks[0],peaks[1]- peaks[0],sped,ssignal,Ratio);
}

//Ron stopped being here
