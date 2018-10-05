// Project include files.
//#include "ns.h"

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

double fitFunc(double *x,double *p);
double kMathPi = 3.1415926; //3.1415926
//void addGraph(double (*func)(double,double*),int color,double *p,double xmin,double xmax);


inline double gauss(double x,double sigma)
{
	// Normalized gaussian.
	return exp(-x*x/(2*sigma*sigma))/(sqrt(2.*kMathPi)*sigma);
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


void SPE_fit_tester_macro(){
  string path_to_read = "SPE_Root/";
  //string file_to_read = "PMT_NB0167_HV1530_Analysis.root";
  //string file_to_read = "PMT_NB0169_HV1450_Analysis.root";
  //string file_to_read = "PMT_NB0170_HV1550_Analysis.root";
  //string file_to_read = "PMT_NB0171_HV1610_Analysis.root";
  string file_to_read = "PMT_NB0166_HV1370_Analysis.root"; 

  file_to_read = path_to_read + file_to_read;
  TFile *FileLoad  = new TFile(file_to_read.c_str(),"read");

  TH1D*TestHist=0;
  FileLoad->GetObject("SPE2",TestHist); 
  if (TestHist == 0) cout << "TEST HIST IS ZERO" << endl;

  TestHist -> GetYaxis()->SetTitle("Counts ");
  TestHist -> GetYaxis()->SetTitleOffset(1.5);
  TestHist -> GetXaxis()->SetTitle("Charge (mv*ns)");

  // Draw spectrum.
  TestHist -> Draw();

  // Fit the data.
  double fitXMin = 40.;
		
  double fitXMax = 2000.;

  TF1 *tf11=new TF1("data1Fit",fitFunc,fitXMin,fitXMax,7);
  tf11->SetLineWidth(3);

  // Exponential background.
  tf11->SetParameter(0,4.77335);
  tf11->SetParameter(1,91.6516);		
  tf11->SetParameter(2, 0.0000656504);
  //tf11->SetParameter(2, 2.0);


  // Single photon gaussian.
  tf11->SetParameter(3,195784);
  tf11->SetParameter(4,158.745);
  tf11->SetParameter(4,350.745);
  tf11->SetParameter(5,161.693);

  // Double and triple photon gaussian.
  tf11->SetParameter(6, 0.223033);

  // Perform fit.
  TFitResultPtr tfrp1 = TestHist -> Fit("data1Fit","RSE");

  addGraph(expBackground,kMagenta+2,tf11->GetParameters(),fitXMin,fitXMax);
  addGraph(photon1Signal,kGreen+2,tf11->GetParameters(),fitXMin,fitXMax);
  addGraph(photon2Signal,kCyan+2,tf11->GetParameters(),fitXMin,fitXMax);
  // Print results.
  tfrp1->Print();

}




