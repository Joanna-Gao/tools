// To use this file, do root -l 'HistOverlay.C("rockFileName.root",
// "airFileName.root", "outputFileName.pdf")'

// C++ includes
#include  <iostream>

//ROOT Includes
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void HistOverlay(const char *rockFile, const char *airFile,
                 const char *outputFile)
{
  TCanvas *can = new TCanvas("can", "can", 1500, 600);
  can->Divide(2,1);
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetPalette(kSolar);
  can->SetLogy();
  can->SetRightMargin(1);

  TFile *rock = new TFile(rockFile, "READ");
  TFile *air = new TFile(airFile, "READ");

  //*********************
  // For comparison plots
  //*********************

  ////TH1D *rockHit = (TH1D*)rock->Get("sumQHisto");
  //TH2D *rockHit = (TH2D*)rock->Get("digiQvsTHistoID");
  //rockHit->Rebin2D();
  //rockHit->SetXTitle("Hit Time (ns)");
  //rockHit->SetYTitle("Charge");
  //rockHit->GetYaxis()->SetTitleOffset(1);
  //rockHit->SetMarkerStyle(kOpenCircle);
  //rockHit->SetMarkerColor(29);
  //rockHit->SetLineColor(29);
  ////rockHit->SetFillColor(29);  // set colour for box plots
  ////can->cd(1);
  //rockHit->Draw("BOX PMC PLC");

  ////TH1D *airHit = (TH1D*)air->Get("sumQHisto");
  //TH2D *airHit = (TH2D*)air->Get("digiQvsTHistoID");
  //airHit->Rebin2D();
  //airHit->SetMarkerStyle(kOpenSquare);
  //airHit->SetMarkerColor(46);
  //airHit->SetLineColor(46);
  ////airHit->SetFillColor(46);  // set colour for box plots 
  ////can->cd(3);
  //airHit->Draw("BOX SAME PMC PLC"); // PFC for 2D PLC PMC

  //TLegend *legend = new TLegend(0.55,0.75,0.9,0.9);
  //legend->AddEntry(rockHit, "Rock Surrounding the Tank", "l");
  //legend->SetTextSize(0.03);
  //legend->AddEntry(airHit, "Air Surrounding the Tank", "l");
  //legend->SetTextSize(0.03);
  //legend->Draw();

  //*********************


  // renormalise the rock and air plot and plot the difference in each bin

  TH2D *rockHitNorm = (TH2D*)rock->Get("digiQvsTHistoID");
  rockHitNorm->Rebin2D();
  //rockHitNorm->Scale(1/rockHitNorm->Integral());
  //can->cd(1);
  //rockHitNorm->Draw("BOX");
  
  TH2D *airHitNorm = (TH2D*)air->Get("digiQvsTHistoID");
  airHitNorm->Rebin2D();
  //airHitNorm->Scale(1/airHitNorm->Integral());

  // plot the difference histogram
  TH2D *HitDiff = new TH2D(*airHitNorm);
  HitDiff->SetNameTitle("HitDiff", "Difference in Hits Between Air and Rock;\
  Hit Time (ns);Charge");
  if (!(HitDiff->GetSumw2N() > 0)) // ensure proper error propagation
    HitDiff->Sumw2(kTRUE);
  Double_t scale_factor = -1.;
  HitDiff->Add(rockHitNorm, scale_factor);
  double minval = abs(HitDiff->GetMinimum());
  double maxval = abs(HitDiff->GetMaximum());
  double benchMark = ((minval > maxval) ? minval : maxval) * 0.1;
  
  cout << "Bench mark: " << benchMark << endl;

  // Disregard the histogram bins that are smaller than the bench mark, which
  // is 10% of the abs(max) or abs(min) values depends on which one is bigger
  for (int x = 0; x < 50; ++x)
  {
    for (int y = 0; y < 50; ++y)
    {
      if (abs(HitDiff->GetBinContent(x, y)) <= benchMark)
        HitDiff->SetBinContent(x, y, 0);
    }
  }

  HitDiff->GetZaxis()->SetLabelSize(0.02);
  can->cd(1);
  HitDiff->Draw("BOX");
  can->cd(2);
  HitDiff->Draw("COLZ1");

  //TH2D *HitDiffCheck = (TH2D*)rockHitNorm->Clone();
  //HitDiffCheck->Reset();
  //double rockContent, airContent;
  //// Subtract the bin content manually for cross check
  //for (int x = 0; x < 50; ++x)
  //{
  //  for (int y = 0; y < 50; ++y)
  //  {
  //    rockContent = rockHitNorm->GetBinContent(x, y);
  //    airContent = airHitNorm->GetBinContent(x, y);
  //    HitDiffCheck->SetBinContent(x, y, (airContent - rockContent));
  //  }
  //}
  //can->cd(2);
  //HitDiffCheck->Draw("BOX");

  //can->SaveAs(outputFile);
}
