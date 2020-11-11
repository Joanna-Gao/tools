// To use this file, do root -l 'HistOverlay.C("rockFileName.root","airFileName.root")'

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
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(kSolar);
  can->SetLogy();
  can->cd();

  TFile *rock = new TFile(rockFile, "READ");
  TFile *air = new TFile(airFile, "READ");

  //TH1D *rockHit = (TH1D*)rock->Get("sumQHisto");
  TH2D *rockHit = (TH2D*)rock->Get("digiQvsTHistoID");
  rockHit->SetXTitle("Hit Time (ns)");
  rockHit->SetYTitle("Charge");
  rockHit->GetYaxis()->SetTitleOffset(1);
  //rockHit->SetMarkerStyle(kOpenCircle);
  //rockHit->SetLineColor(2);
  rockHit->Draw("PFC PLC PMC");

  //TH1D *airHit = (TH1D*)air->Get("sumQHisto");
  TH2D *airHit = (TH2D*)air->Get("digiQvsTHistoID");  
  airHit->SetLineColor(2);
  //airHit->SetMarkerStyle(kOpenSquare);
  airHit->Draw("SAME PFC PLC PMC"); // PFC for 2D PLC PMC

  TLegend *legend = new TLegend(0.55,0.75,0.9,0.9);
  legend->AddEntry(rockHit, "Rock Surrounding the Tank", "l");
  legend->SetTextSize(0.03);
  legend->AddEntry(airHit, "Air Surrounding the Tank", "l");
  legend->SetTextSize(0.03);
  legend->Draw();

  //can->SaveAs(outputFile);
}
