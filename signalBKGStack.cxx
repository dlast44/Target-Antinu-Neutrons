//File: signalBKGStack.cxx
//Info: This is a script to run a loop over all MC breakdown plots in a single histos file and save nice plots from them. Primarily used for validation against older plots.
//
//Usage: signalBKGStack <histos_file> <output_directory> <plot_label>
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//C++ includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/stat.h>

//ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

TCanvas* DrawIntType(string name_QE, TFile* inFile, TString sample){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  MnvH1D* h_QE_Sig = (MnvH1D*)inFile->Get((TString)name_QE);
  h_QE_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  string name_sig = (string)h_QE_Sig->GetName();
  name_sig.erase(name_sig.length()-3,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-12,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_QE_Sig->GetTitle();
  TString Xtitle = h_QE_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_QE_Sig->GetYaxis()->GetTitle();
  /*
  cout << title << endl;
  title.erase(0, 8);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_RES_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_RES");
  h_RES_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_DIS_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_DIS");
  h_DIS_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_2p2h_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_2p2h");
  h_2p2h_Sig->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Sig->SetFillColor(TColor::GetColor("#44AA99"));

  MnvH1D* h_Other_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_Other");
  h_Other_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_QE_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_IntType_QE");
  h_QE_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillStyle(3444);

  MnvH1D* h_RES_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_IntType_RES");
  h_RES_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillStyle(3444);

  MnvH1D* h_DIS_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_IntType_DIS");
  h_DIS_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillStyle(3444);

  MnvH1D* h_2p2h_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_IntType_2p2h");
  h_2p2h_Bkg->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_IntType_Other");
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillStyle(3444);

  THStack* h = new THStack();
  h->Add(h_Other_Bkg);
  h->Add(h_2p2h_Bkg);
  h->Add(h_DIS_Bkg);
  h->Add(h_RES_Bkg);
  h->Add(h_QE_Bkg);

  h->Add(h_Other_Sig);
  h->Add(h_2p2h_Sig);
  h->Add(h_DIS_Sig);
  h->Add(h_RES_Sig);
  h->Add(h_QE_Sig);

  h->Draw("hist");
  c1->Update();

  h->SetTitle(sample);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitle(Ytitle);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.075);
  
  size_t pos = 0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.045);
  }

  h->Draw("hist");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(h_QE_Sig,"Sig. + QE");
  leg->AddEntry(h_QE_Bkg,"Bkg. + QE");

  leg->AddEntry(h_RES_Sig,"Sig. + RES");
  leg->AddEntry(h_RES_Bkg,"Bkg. + RES");

  leg->AddEntry(h_DIS_Sig,"Sig. + DIS");
  leg->AddEntry(h_DIS_Bkg,"Bkg. + DIS");

  leg->AddEntry(h_2p2h_Sig,"Sig. + 2p2h");
  leg->AddEntry(h_2p2h_Bkg,"Bkg. + 2p2h");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->Draw();
  c1->Update();
  return c1;
}

TCanvas* DrawTargetType(string name_C, TFile* inFile, TString sample){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  MnvH1D* h_C_Sig = (MnvH1D*)inFile->Get((TString)name_C);
  h_C_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  string name_sig = (string)h_C_Sig->GetName();
  name_sig.erase(name_sig.length()-2,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-15,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_C_Sig->GetTitle();
  TString Xtitle = h_C_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_C_Sig->GetYaxis()->GetTitle();
  /*
  cout << title << endl;
  title.erase(0, 7);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_Fe_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_Fe");
  h_Fe_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Pb_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_Pb");
  h_Pb_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_O_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_O");
  h_O_Sig->SetLineColor(TColor::GetColor("#332288"));
  h_O_Sig->SetFillColor(TColor::GetColor("#332288"));

  MnvH1D* h_H_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_H");
  h_H_Sig->SetLineColor(TColor::GetColor("#DDCC77"));
  h_H_Sig->SetFillColor(TColor::GetColor("#DDCC77"));

  //h_Prot_Sig->SetFillColor(TColor::GetColor("#999933"));

  MnvH1D* h_Other_Sig = (MnvH1D*)inFile->Get((TString)name_sig+"_Other");
  h_Other_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  //h_None_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_C_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_C");
  h_C_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillStyle(3444);

  MnvH1D* h_Fe_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_Fe");
  h_Fe_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillStyle(3444);

  MnvH1D* h_Pb_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_Pb");
  h_Pb_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillStyle(3444);

  MnvH1D* h_O_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_O");
  h_O_Bkg->SetLineColor(TColor::GetColor("#332288"));
  h_O_Bkg->SetFillColor(TColor::GetColor("#332288"));
  h_O_Bkg->SetFillStyle(3444);

  MnvH1D* h_H_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_H");
  h_H_Bkg->SetLineColor(TColor::GetColor("#DDCC77"));
  h_H_Bkg->SetFillColor(TColor::GetColor("#DDCC77"));
  h_H_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)inFile->Get((TString)name_bkg+"_bkg_TargetType_Other");
  h_Other_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillStyle(3444);

  THStack* h = new THStack();
  h->Add(h_Other_Bkg);
  h->Add(h_H_Bkg);
  h->Add(h_O_Bkg);
  h->Add(h_Pb_Bkg);
  h->Add(h_Fe_Bkg);
  h->Add(h_C_Bkg);

  h->Add(h_Other_Sig);
  h->Add(h_H_Sig);
  h->Add(h_O_Sig);
  h->Add(h_Pb_Sig);
  h->Add(h_Fe_Sig);
  h->Add(h_C_Sig);

  h->Draw("hist");
  c1->Update();

  h->SetTitle(sample+" Target Type Breakdown");//+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitle(Ytitle);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.075);
  
  size_t pos=0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"None");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.045);
  }

  h->Draw("hist");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(h_C_Sig,"Sig. + C");
  leg->AddEntry(h_C_Bkg,"Bkg. + C");

  leg->AddEntry(h_Fe_Sig,"Sig. + Fe");
  leg->AddEntry(h_Fe_Bkg,"Bkg. + Fe");

  leg->AddEntry(h_Pb_Sig,"Sig. + Pb");  
  leg->AddEntry(h_Pb_Bkg,"Bkg. + Pb");
  
  leg->AddEntry(h_O_Sig,"Sig. + O");
  leg->AddEntry(h_O_Bkg,"Bkg. + O");

  leg->AddEntry(h_H_Sig,"Sig. + H");
  leg->AddEntry(h_H_Bkg,"Bkg. + H");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->Draw();
  c1->Update();
  return c1;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 4) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string fileName=string(argv[1]);
  string outDir=string(argv[2]);
  TString label=argv[3];

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = fileName;
  size_t pos=0;

  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "Input need be .root file." << endl;
    return 4;
  }

  cout << "Input Signal file name parsed to: " << fileNameStub << endl;

  TFile* inFile = new TFile(fileName.c_str(),"READ");

  TList* keyList = inFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    //cout << key->GetName() << endl;
    pos=0;
    string name=(string)key->GetName();
    if((pos = name.find("_sig_IntType_QE")) != string::npos){
      TCanvas* c1 = DrawIntType(name,inFile,label);
      name.erase(name.length()-15,name.length());
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked_test.pdf");
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked_test.png");
      cout << "" << endl;
      delete c1;
    }
    else if ((pos = name.find("_sig_TargetType_C")) != string::npos){
      TCanvas* c1 = DrawTargetType(name,inFile,label);
      name.erase(name.length()-17,name.length());
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked_test.pdf");
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked_test.png");
      cout << "" << endl;
      delete c1;
    }
    /*
    else if ((pos = name.find("_sig_Type_C")) != string::npos){
      cout << "FOUND: " << name << endl;
      name.erase(name.length()-18,name.length());
      cout << "Parsed to: " << name  << endl;
    }
    */
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
