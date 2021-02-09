// C++ includes
#include <iostream>
#include <vector>

//ROOT Includes
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void MakeHist(const char *rockFile, const char *airFile)
//              const char *energyFile)
{
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1);//(kSolar);
  can->SetRightMargin(1);

  TFile *rock = new TFile(rockFile, "READ");
  TFile *air = new TFile(airFile, "READ");
  //TFile *energy = new TFile(energyFile, "READ");

  TTree *rockTree = (TTree*)rock->Get("IDTree");
  TTree *airTree = (TTree*)air->Get("IDTree");
  //TTree *energyTree = (TTree*)energy->Get("IDTree");

  //rockTree->Draw("digiHitQ:digiHitT>>rockQvsT(50,0,4000,50,0,100", "", "goff");
  //TH2D *rockQvsT = (TH2D*)gDirectory->Get("rockQvsT");
  //if (&normalise) rockQvsT->Scale(1/rockQvsT->Integral());
  //
  //airTree->Draw("digiHitQ:digiHitT>>airQvsT(50,0,4000,50,0,100", "", "goff");
  //TH2D *airQvsT = (TH2D*)gDirectory->Get("airQvsT"); 
  //if (&normalise) airQvsT->Scale(1/airQvsT->Integral());

  //// plot the difference histogram
  //TH2D *HitDiff = new TH2D(*airQvsT);
  //HitDiff->SetNameTitle("HitDiff", "Difference in Hits Between Air and Rock;\
  //Hit Time (ns);Charge");
  //if (!(HitDiff->GetSumw2N() > 0)) // ensure proper error propagation
  //  HitDiff->Sumw2(kTRUE);
  //Double_t scale_factor = -1.;
  //HitDiff->Add(rockQvsT, scale_factor);

  //double minval = abs(HitDiff->GetMinimum());
  //double maxval = abs(HitDiff->GetMaximum());
  //double benchMark = ((minval > maxval) ? minval : maxval) * 0.1;
  //
  //cout << "Bench mark: " << benchMark << endl;

  //// Disregard the histogram bins that are smaller than the bench mark, which
  //// is 10% of the abs(max) or abs(min) values depends on which one is bigger
  //for (int x = 0; x < 50; ++x)
  //{
  //  for (int y = 0; y < 50; ++y)
  //  {
  //    if (abs(HitDiff->GetBinContent(x, y)) <= benchMark)
  //      HitDiff->SetBinContent(x, y, 0);
  //  }
  //}

  //HitDiff->GetZaxis()->SetLabelSize(0.02);
  //HitDiff->Draw("COLZ1");

  //can->SaveAs(outputFile);
  
  rockTree->Draw("digiHitQHistoID>>>ht_rock(100, 0, 100)");
  TH1D *ht_rock = (TH1D*)gDirectory->Get("ht_rock");
  //ht_rock->SetNameTitle("ht_rock", "Charge;Charge (p.e.);");
  ht_rock->SetLineColor(2);

  airTree->Draw("digiHitQHistoID>>>ht_air(100, 0, 100)");
  TH1D *ht_air = (TH1D*)gDirectory->Get("ht_air");
  ht_air->SetLineColor(4);
 
  //energyTree->Draw("digiHitT>>ht_energy(100, 0, 4000)");
  //TH1D *ht_energy = (TH1D*)gDirectory->Get("ht_energy");
  //ht_energy->SetLineColor(6);

  ht_rock->Draw("PLC");
  ht_air->Draw("SAME PLC");
  //ht_energy->Draw("SAME PLC");

  TLegend *legend = new TLegend(0.55,0.75,0.9,0.9);
  legend->AddEntry(ht_rock, "Rock 100 GeV", "l");
  legend->SetTextSize(0.03);
  legend->AddEntry(ht_air, "Air 100 GeV", "l");
  legend->SetTextSize(0.03);
  //legend->AddEntry(ht_energy, "Rock 1000 GeV", "l");
  //legend->SetTextSize(0.03);
  legend->Draw();
}
