#include<iostream>
#include<fstream>
#include<vector>
#include"TMath.h"

#include"TFile.h"
#include"TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDerivative.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
// Result

typedef struct {

double value;
double error;


} ValueWithError;

typedef struct {

  ValueWithError ped;
  ValueWithError fped;
  ValueWithError pemean;
  ValueWithError pewidth;
  ValueWithError npe;
  ValueWithError peak;
  ValueWithError valley;
  ValueWithError peakToValley;
  ValueWithError fvalley;
  double rate;
  
} Result;

void fillValueWithError(ValueWithError* val,TH1F* h){
  val->value = h->GetMean();
  val->error = h->GetRMS();
}

void fillValueWithError(ValueWithError* val,RooRealVar* var){
  val->value = var->getVal();
  val->error = var->getError();
}


#include <tuple>
typedef std::tuple<double,double,double,double,double> InitParams;

using namespace RooFit;

InitParams initializeFit(TH1F* h){

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

//  std::get<0>(t) 

const char* varname(std::string header, std::string var){
  return (header + var).c_str();
}

RooAddPdf* makePMTPDF(RooRealVar* counts,double pmval, double psval, double mval, double sval, double f1peval, double dmval = 15, double dsval = 20, double fdval = 0.9){

  std::cout << pmval << " " << psval << " "  << mval << " "  << sval  << " "  <<f1peval << std::endl;


 RooRealVar* dm = new RooRealVar("dmean","dmean",dmval, 0, 100 );   // pedestal position
 RooRealVar* ds = new RooRealVar("dsigma","dsigma",dsval,0,200 );  // pedestal sigma --not sure how to set generically

 RooGaussian* dgauss = new  RooGaussian("dgauss","dgauss", *counts, *dm, *ds);

  
 RooRealVar* pedm = new RooRealVar("pedmean","pedmean",pmval, 0, pmval + 3*psval );   // pedestal position
 RooRealVar* peds = new RooRealVar("pedsigma","pedsigma",psval,0,5*psval );  // pedestal sigma --not sure how to set generically

 RooGaussian* pedgauss = new  RooGaussian("pedgauss","pedgauss", *counts, *pedm, *peds);

 // 1 pe
 RooRealVar* m0 = new RooRealVar("mean0","mean0",mval,  mval - 10*sval,  mval  + 10*sval);   // pedestal position
 RooRealVar* s0 = new RooRealVar("sigma0","sigma0", sval,0, 10*sval);  // pedestal sigma --not sure how to set generically



 RooRealVar* npe = new RooRealVar ("npe","npe", f1peval,0.2,2); // number of photo electrons

 std::string form1 = "mean0 + pedmean"; 
 RooFormulaVar* m =  new RooFormulaVar("mean","mean",form1.c_str(),RooArgSet(*m0,*pedm));
 std::string form2 =  "sqrt(sigma0*sigma0 + pedsigma*pedsigma)" ; 
 RooFormulaVar* s =  new RooFormulaVar("sigma","sigma",form2.c_str(),RooArgSet(*s0,*peds)); 
 RooGaussian* gauss = new  RooGaussian("gauss1","gauss1", *counts, *m, *s);

 std::string form3 = "2*mean0 + pedmean"; 
 RooFormulaVar* m2 =  new RooFormulaVar("mean2","mean2",form3.c_str(),RooArgSet(*m0,*pedm));
 std::string form4 =  "sqrt(2*sigma0*sigma0 + pedsigma*pedsigma)" ; 
 RooFormulaVar* s2 =  new RooFormulaVar("sigma2","sigma2",form4.c_str(),RooArgSet(*s0,*peds)); 
 RooGaussian* gauss2 = new  RooGaussian("gauss2","gauss2", *counts, *m2, *s2);

 std::string form5 = "3*mean0 +pedmean"; 
 RooFormulaVar* m3 =  new RooFormulaVar("mean3","mean3",form5.c_str(),RooArgSet(*m0,*pedm));
 std::string form6 =  "sqrt(3*sigma0*sigma0 + pedsigma*pedsigma)" ; 
 RooFormulaVar* s3 =  new RooFormulaVar("sigma3","sigma3",form6.c_str(),RooArgSet(*s0,*peds)); 
 RooGaussian* gauss3 = new  RooGaussian("gauss3","gauss3", *counts, *m3, *s3);

 std::string form7 = "4*mean0 +pedmean"; 
 RooFormulaVar* m4 =  new RooFormulaVar("mean4","mean4",form7.c_str(),RooArgSet(*m0,*pedm));
 std::string form8 =  "sqrt(4*sigma0*sigma0 + pedsigma*pedsigma)" ; 
 RooFormulaVar* s4 =  new RooFormulaVar("sigma4","sigma4",form8.c_str(),RooArgSet(*s0,*peds)); 
 RooGaussian* gauss4 = new  RooGaussian("gauss4","gauss4", *counts, *m4, *s4);

 
 //RooRealVar* frac_ped1 = new RooRealVar("frac_ped1", "frac_ped1", 0.5,0,1);
 RooRealVar* fd = new RooRealVar("fd","fd",fdval,0,1.01);

 RooFormulaVar*  frac_ped1 = new RooFormulaVar("frac_ped1", "frac_ped1",  "TMath::Exp(-npe)", *npe); //Poisson for 0 Photo electrons
 RooFormulaVar* frac_pe1 = new RooFormulaVar("frac_pe1", "frac_pe1",  "TMath::Poisson(1,npe) ", *npe); // ""   for 1   ""
 RooFormulaVar* frac_pe2 = new RooFormulaVar("frac_pe2", "frac_pe2",  "TMath::Poisson(2,npe) ",  *npe);
 RooFormulaVar* frac_pe3 = new RooFormulaVar("frac_pe3", "frac_pe3",  "TMath::Poisson(3,npe) ", *npe);



 /* RooRealVar* frac_ped1 = new RooRealVar("frac_ped1","frac_ped1", TMath::Exp(-0.66),  0,  1  ); 
 RooRealVar* frac_pe1 = new RooRealVar("frac_pe1","frac_pe1", TMath::Poisson(1,0.66),  0,  1  ); 
 RooRealVar* frac_pe2 = new RooRealVar("frac_pe2","frac_pe2", TMath::Poisson(2,0.66),  0,  1  ); 
 RooRealVar* frac_pe3 = new RooRealVar("frac_pe3","frac_pe3", TMath::Poisson(3,0.66),  0,  1  ); 
 */
 
 RooAddPdf* smodel = new RooAddPdf("smodel","smodel",RooArgList(*pedgauss,*gauss,*gauss2,*gauss3,*gauss4), RooArgList(*frac_ped1,*frac_pe1,*frac_pe2,*frac_pe3));

 
 RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*smodel,*dgauss),*fd);
 
 //RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*pedgauss,*gauss), RooArgList(*frac_ped1));
 return model;
}


RooAddPdf* makePMTPDF(RooRealVar* counts,const RooArgList& fitpars){

 RooRealVar* fpedmean = (RooRealVar*)fitpars.find("pedmean");
 RooRealVar* pedsigma = (RooRealVar*)fitpars.find("pedsigma");
 //RooRealVar* pedsigma2 = (RooRealVar*)fitpars.find("pedsigma2");
 RooRealVar* fmean = (RooRealVar*)fitpars.find("mean0");
 RooRealVar* fsigma = (RooRealVar*)fitpars.find("sigma0");
 RooRealVar* f1pe = (RooRealVar*)fitpars.find("npe");
 RooRealVar* dm = (RooRealVar*)fitpars.find("dmean");
 RooRealVar* ds = (RooRealVar*)fitpars.find("dsigma");
 RooRealVar* fd = (RooRealVar*)fitpars.find("fd");

 
 return makePMTPDF( counts,fpedmean->getVal(),pedsigma->getVal(), fmean->getVal(), fsigma->getVal(), f1pe->getVal(), dm->getVal(), ds->getVal(),fd->getVal());
}

RooAddPdf* makePMTPDF(RooRealVar* counts,InitParams& params){
  return makePMTPDF(counts,std::get<0>(params),std::get<2>(params),std::get<1>(params),std::get<3>(params),1-std::get<4>(params));
}



Result* propagateAndFill(RooRealVar* counts,RooAddPdf* model ,RooFitResult* fres){

 Result* res = new Result();
  
 const RooArgList& fitpars = fres->floatParsFinal(); 
 RooRealVar* fpedmean = (RooRealVar*)fitpars.find("pedmean");
 RooRealVar* fmean = (RooRealVar*)fitpars.find("mean0");
 RooRealVar* fsigma = (RooRealVar*)fitpars.find("sigma0");
 RooRealVar* fnpe = (RooRealVar*)fitpars.find("npe");
 // RooRealVar* fped = (RooRealVar*)fitpars.find("f");

 if (fpedmean) fillValueWithError(&res->ped,fpedmean);
 if (fmean) fillValueWithError(&res->pemean,fmean);
 if (fsigma) fillValueWithError(&res->pewidth,fsigma);
 if (fnpe) fillValueWithError(&res->npe,fnpe);
 //fillValueWithError(&res->fped,fped);

 
 //now the complicated ones that require sampling the fitted pdf/covariance
 TH1F* histo = new TH1F("valley", "valley", 200, res->ped.value,  res->pemean.value ); histo->Sumw2();
 TH1F* histo2 = new TH1F("peak", "peak", 200,  res->pemean.value - res->pewidth.value ,  res->pemean.value +3* res->pewidth.value ); histo2->Sumw2();
 TH1F* histo3 = new TH1F("peakToValley", "peakToValley", 100,  0,10 ); histo3->Sumw2();
 TH1F* histo4 = new TH1F("f", "f", 100,  0.,1 ); histo4->Sumw2();

 RooArgSet nset(*counts) ;
 
 
 for (int i = 0; i < 100; ++i){
   const RooArgList& sample = fres->randomizePars();
   // std::string name = "pdf_" + std::to_string(i);
   RooAddPdf* pdf = makePMTPDF(counts,sample);
   TF1* fmodel = pdf->asTF( *counts,fitpars,*counts);
   double vpos = fmodel->GetMinimumX(10,20);
   double ppos = fmodel->GetMaximumX(20,45);
   histo->Fill(vpos);
   histo2->Fill(ppos);
   std::cout << "****" <<  fmodel->GetMaximum(20,45) <<  " " << fmodel->GetMinimum(10,20) << std::endl;
   histo3->Fill(fmodel->Eval(ppos)/fmodel->Eval(vpos));

   //RooArgSet* theset = pdf->getComponents();
   // RooAbsPdf* signalmodel =  (RooAbsPdf*)theset->find("gmodel");;
   counts->setRange("signal",vpos, 520) ;
   
   RooAbsReal* igx_s2 =  pdf->createIntegral(*counts,NormSet(*counts),Range("signal")) ;
   histo4->Fill(igx_s2->getVal());
  
 }

 fillValueWithError(&res->valley,histo);
 fillValueWithError(&res->peak,histo2);
 fillValueWithError(&res->peakToValley,histo3);
 fillValueWithError(&res->fvalley,histo4);
 
 
 // histo4->Draw();
 return res;
}

Result* fitModel(std::stting file, double maxval = 5 ,double max = 1000){

  //  Result* res = new Result();
  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  TH1F *histo =new TH1F("h","h",100, -1., maxval); histo->Sumw2();
  
  TFile *input=new TFile(file.c_str());
  TTree *tree = (TTree*)input->Get("tree");
  double val; tree->SetBranchAddress("charge", &val);
  for (int i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    if (val < maxval) histo->Fill(val);
  }
  
  
  histo->Draw("HISTO");
  //return 0;
  InitParams params = initializeFit(histo);

  //   return 0;
  
  RooRealVar* counts = new RooRealVar("charge", "charge", -5, maxval); 
  RooDataHist data("data", "dataset", *counts , histo);
   
  RooAddPdf* model = makePMTPDF(counts,params);
  
  RooFitResult* fres = model->fitTo(data,Save());

 
 RooPlot* frame = counts->frame();
 data.plotOn(frame);
 //model->plotOn(frame, Components("smodel"),LineColor(kBlue), LineStyle(2));
 // model->plotOn(frame, Components("gauss2"),LineColor(kGreen), LineStyle(2));
 // model->plotOn(frame, Components("dgauss"),LineColor(kBlack), LineStyle(2));
 model->plotOn(frame, Components("gauss1"),LineColor(kBlue), LineStyle(2));
 model->plotOn(frame, Components("gauss2"),LineColor(kGreen), LineStyle(2));
 // model->plotOn(frame, Components("dgauss"),LineColor(kBlack), LineStyle(2));
 model->plotOn(frame, LineColor(kRed));

 model->plotOn(frame, LineColor(kRed));

  std::cout << " plotting " << std::endl;
 TAxis* xachse = frame->GetXaxis();  TAxis* yachse = frame->GetYaxis();
 xachse->SetTitleFont (132);
 yachse->SetTitleFont (132);
 xachse->SetLabelFont (132);
 yachse->SetLabelFont (132);
 yachse->SetTitleOffset(1.1); //was 0.13
 xachse->SetTitleSize(0.06); //was 0.13
 yachse->SetTitleSize(0.06); //was 0.13
 xachse->SetLabelSize(0.06); //was 0.13
 yachse->SetLabelSize(0.06); //was 0.13
 yachse->SetTitle("Entries/ Count");
 xachse->SetTitle("Counts");
 //frame->SetMaximum(max);
 frame->Draw();
 
 
 // Result* res = propagateAndFill(counts,model,fres);

 //Result* res;
 return 0;
}

