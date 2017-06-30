#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"
#include "boost/python/list.hpp"
#include "boost/python/str.hpp"
//#include "boost/filesystem.hpp"
#include <iostream>
#include <stdint.h>
#include "TString.h"
#include <string>
#include <boost/python/exception_translator.hpp>
#include <exception>
#include "../interface/pythonToSTL.h"
#include "friendTreeInjector.h"
#include "TROOT.h"
#include "colorToTColor.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include <algorithm>

using namespace boost::python; //for some reason....

TH1F* getProfile(TH2F* scatterHist, TString profileType)
{
  TH1F *hProfile = new TH1F("hProfile","",
                         scatterHist->GetNbinsX(),
                         scatterHist->GetXaxis()->GetXmin(),
                         scatterHist->GetXaxis()->GetXmax());
  hProfile->Sumw2();

  TH1F *htemp = new TH1F("htemp","",
                         scatterHist->GetNbinsY(),
                         scatterHist->GetYaxis()->GetXmin(),
                         scatterHist->GetYaxis()->GetXmax());
  htemp->Sumw2();

  const int nq = 2;
  double yq[2],xq[2];
  xq[0] = 0.5;
  xq[1] = 0.90;
  double yq_IQM[2],xq_IQM[2];
  xq_IQM[0] = 0.25;
  xq_IQM[1] = 0.75;
  // Loop over x bins
  for(int xBin = 1; xBin <= scatterHist->GetNbinsX(); xBin++) {
    htemp->Reset();
    // Copy y bin content in this x bins to htemp
    for(int yBin = 1; yBin <= scatterHist->GetNbinsY(); yBin++) {
      htemp->SetBinContent(yBin,scatterHist->GetBinContent(xBin,yBin));
      htemp->SetBinError(yBin,scatterHist->GetBinError(xBin,yBin));
    }
    htemp->SetEntries(htemp->GetEffectiveEntries());
    htemp->GetXaxis()->SetRange(0,-1);

    if(htemp->GetEntries()>100){
      double mean = htemp->GetMean(); 
      double meanerror = htemp->GetMeanError();
      double width = htemp->GetRMS();
      if(width < 0.1) width = 0.1;
      if(htemp->GetSumOfWeights() <= 0) {
        continue; 
      } else {
        htemp->Fit("gaus","QNI","", mean - 3 * width,mean + 3 * width);
        TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
        mean = f->GetParameter(1);
        meanerror = f->GetParError(1);
        width = f->GetParameter(2);
        if(width < 0.05) width = 0.05;
        if( (htemp->Fit(f,"QI","goff",mean - 1.5 * width, mean + 1.5 * width) == 0) ) { //removed N option to store GaussFit with histogram
          mean = f->GetParameter(1);
          meanerror = f->GetParError(1);
          width = f->GetParameter(2);
              
          if(profileType=="GaussFitMean"){
          hProfile->SetBinContent(xBin,mean);
          hProfile->SetBinError(xBin,meanerror);
          }
          else if (profileType=="GaussFitWidth"){
            hProfile->SetBinContent(xBin,width);
            hProfile->SetBinError(xBin,f->GetParError(2));
          }
          else if  (profileType=="GaussFitWidthNormByMean"){ 
            hProfile->SetBinContent(xBin,width/mean);
            hProfile->SetBinError(xBin,f->GetParError(2)/mean);
          }
          else if  (profileType=="GaussFitChi2"){ 
            hProfile->SetBinContent(xBin,f->GetChisquare() / f->GetNumberFreeParameters());
            hProfile->SetBinError(xBin,0.01);
          }
          else if  (profileType=="GaussFitProb"){ 
            hProfile->SetBinContent(xBin,f->GetProb());
            hProfile->SetBinError(xBin,0.01);
          }
          mean = htemp->GetMean();
          meanerror = htemp->GetMeanError();
          width = htemp->GetRMS();
          htemp->ComputeIntegral(); //needed for TH1::GetQuantiles() according to documentation
          htemp->GetQuantiles(nq,yq,xq);

          if(profileType=="Mean"){
          hProfile->SetBinContent(xBin,mean);
          hProfile->SetBinError(xBin,meanerror);
          }
          else if (profileType=="StandardDeviationNormByMean"){
            hProfile->SetBinContent(xBin,width/mean);
            hProfile->SetBinError(xBin,htemp->GetRMSError()/mean);
          }
          else if (profileType=="StandardDeviation"){
            hProfile->SetBinContent(xBin,width);
            hProfile->SetBinError(xBin,htemp->GetRMSError());
          }
          else if (profileType=="Median"){
            hProfile->SetBinContent(xBin,yq[0]);
            hProfile->SetBinError(xBin,0.0001);//not meaningful
          }
          htemp->GetQuantiles(nq,yq_IQM,xq_IQM);
          Int_t IQ_low_bin_i = htemp->FindBin(yq_IQM[0]);
          Int_t IQ_hig_bin_i = htemp->FindBin(yq_IQM[1]);
          htemp->GetXaxis()->SetRange(IQ_low_bin_i,IQ_hig_bin_i);
          mean = htemp->GetMean();
          meanerror = htemp->GetMeanError();
          if(profileType=="IQMean"){
          hProfile->SetBinContent(xBin,mean);
          hProfile->SetBinError(xBin,meanerror);
          }
        }



      }//endweightelse
    }
  }
  return hProfile;


}


void makePlots(
        const boost::python::list intextfiles,
        const boost::python::list names,
        const boost::python::list variables,
        const boost::python::list cuts,
        const boost::python::list colors,
        std::string outfile,
        std::string xaxis,
        std::string yaxis,
        bool normalized,
        std::string profileType="None",
        float OverrideMin=-1e100,
        float OverrideMax=1e100) {


    std::vector<TString>  s_intextfiles=toSTLVector<TString>(intextfiles);
    std::vector<TString>  s_vars = toSTLVector<TString>(variables);
    std::vector<TString>  s_names = toSTLVector<TString>(names);
    std::vector<TString>  s_colors = toSTLVector<TString>(colors);
    std::vector<TString>  s_cuts = toSTLVector<TString>(cuts);

    //reverse to make the first be on top
    std::reverse(s_intextfiles.begin(),s_intextfiles.end());
    std::reverse(s_vars.begin(),s_vars.end());
    std::reverse(s_names.begin(),s_names.end());
    std::reverse(s_colors.begin(),s_colors.end());
    std::reverse(s_cuts.begin(),s_cuts.end());

    TString tProfileType = profileType;
    if(!(
          (tProfileType=="None") ||
          (tProfileType=="GaussFitMean") ||
          (tProfileType=="GaussFitWidth") ||
          (tProfileType=="GaussFitWidthNormByMean") || 
          (tProfileType=="GaussFitChi2") || 
          (tProfileType=="GaussFitProb") || 
          (tProfileType=="Mean") ||
          (tProfileType=="StandardDeviationNormByMean") ||
          (tProfileType=="StandardDeviation") ||
          (tProfileType=="Median") ||
          (tProfileType=="IQMean") ))
      throw std::runtime_error("makePlots: profileType not known; None|| GaussFitMean|| GaussFitWidth|| GaussFitWidthNormByMean) || GaussFitChi2 || GaussFitProb || Mean|| StandardDeviationNormByMean|| StandardDeviation|| Median|| IQMean");
    bool makeProfile = false;
    if(tProfileType!="None") makeProfile = true;



    TString toutfile=outfile;
    if(!toutfile.EndsWith(".pdf"))
        throw std::runtime_error("makePlots: output files need to be pdf format");


    if(!s_names.size())
        throw std::runtime_error("makePlots: needs at least one legend entry");
    /*
     * Size checks!!!
     */
    if(s_intextfiles.size() !=s_names.size()||
            s_names.size() != s_vars.size() ||
            s_names.size() != s_colors.size()||
            s_names.size() != s_cuts.size())
        throw std::runtime_error("makePlots: input lists must have same size");

    //make unique list of infiles
    std::vector<TString> u_infiles;
    std::vector<TString> aliases;
    TString oneinfile="";
    bool onlyonefile=true;
    for(const auto& f:s_intextfiles){
        if(oneinfile.Length()<1)
            oneinfile=f;
        else
            if(f!=oneinfile)
                onlyonefile=false;
    }
    for(const auto& f:s_intextfiles){
        //if(std::find(u_infiles.begin(),u_infiles.end(),f) == u_infiles.end()){
        u_infiles.push_back(f);
        TString s="";
        s+=aliases.size();
        aliases.push_back(s);
        //	std::cout << s <<std::endl;
        //}
    }



    friendTreeInjector injector;
    for(size_t i=0;i<u_infiles.size();i++){
        if(!aliases.size())
            injector.addFromFile((TString)u_infiles.at(i));
        else
            injector.addFromFile((TString)u_infiles.at(i),aliases.at(i));
    }
    injector.createChain();

    TChain* c=injector.getChain();
    std::vector<TH1F*> allhistos;
    TLegend * leg=new TLegend(0.2,0.75,0.8,0.88);
    leg->SetBorderSize(0);

    leg->SetNColumns(3);
    leg->SetFillStyle(0);

    TString addstr="";
    if(normalized)
        addstr="normalized";
    //    if(makeProfile)
    //        addstr+="prof";
    float max=-1e100;
    float min=1e100;

    TString tfileout=toutfile;
    tfileout=tfileout(0,tfileout.Length()-4);
    tfileout+=".root";

    TFile * f = new TFile(tfileout,"RECREATE");
    gStyle->SetOptStat(0);

    for(size_t i=0;i<s_names.size();i++){
        TString tmpname="hist_";
        tmpname+=i;
        c->Draw(s_vars.at(i)+">>"+tmpname,s_cuts.at(i),addstr);

        TH1F *histo;
        if(makeProfile){
          std::cout << "pulling in profile of type" << tProfileType << std::endl;
          TH2F *scatterHist = (TH2F*) gROOT->FindObject(tmpname); 
          histo = getProfile(scatterHist,tProfileType );
        }
        else histo = (TH1F*) gROOT->FindObject(tmpname);
        histo->SetLineColor(colorToTColor(s_colors.at(i)));
        histo->SetLineStyle(lineToTLineStyle(s_colors.at(i)));
        histo->SetTitle(s_names.at(i));
        histo->SetName(s_names.at(i));

        histo->SetFillStyle(0);
        histo->SetLineWidth(2);

        float tmax=histo->GetMaximum();
        float tmin=histo->GetMinimum();
        if(tmax>max)max=tmax;
        if(tmin<min)min=tmin;
        if(makeProfile &&OverrideMin!=-1e100){
          //std::cout << "overriding min/max"<< std::endl;
          max = OverrideMax;
          min = OverrideMin;
        }
        //std::cout << "min" << min << " max" << max << std::endl;

        allhistos.push_back(histo);

        histo->Write();


    }
    for(size_t i=allhistos.size();i;i--){
        leg->AddEntry(allhistos.at(i-1),s_names.at(i-1),"l");
    }

    TCanvas cv("plots");

    allhistos.at(0)->Draw("AXIS");
    allhistos.at(0)->GetYaxis()->SetRangeUser(min,1.3*max); //space for legend on top

    allhistos.at(0)->GetXaxis()->SetTitle(xaxis.data());
    allhistos.at(0)->GetYaxis()->SetTitle(yaxis.data());

    allhistos.at(0)->Draw("AXIS");
    for(size_t i=0;i<s_names.size();i++){
        allhistos.at(i)->Draw("same,hist");
    }
    leg->Draw("same");

    cv.Write();
    cv.Print(toutfile);

    f->Close();

}


void makeProfiles(
        const boost::python::list intextfiles,
        const boost::python::list names,
        const boost::python::list variables,
        const boost::python::list cuts,
        const boost::python::list colors,
        std::string outfile,
        std::string xaxis,
        std::string yaxis,
        std::string profileType,
        bool normalized,float minimum, float maximum) {

  makePlots(
        intextfiles,
        names,
        variables,
        cuts,
        colors,
        outfile,
        xaxis,
        yaxis,
        normalized,
        profileType,
        minimum,
        maximum);
}  



// Expose classes and methods to Python
BOOST_PYTHON_MODULE(c_makePlots) {
    //__hidden::indata();//for some reason exposing the class prevents segfaults. garbage collector?
    //anyway, it doesn't hurt, just leave this here
    def("makePlots", &makePlots);
    def("makeProfiles", &makeProfiles);

}

