#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void simple_plot_overlap_air_rock(){
  
 //prelim1
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitle(0);
  gStyle->SetPalette(1);
  gStyle->SetFillColor(1);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  //gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  //gStyle->SetTitleSize(0.06,"Z");
  gStyle->SetPadBorderMode(0);
  //gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09); 
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);

  int nb=1;
  
  //open Rock data
  TFile *f = new TFile("HistoRockRandomMu10GeV1000Evnts.root");
  TTree *rockTree = (TTree*)f->Get("IDTree"); 
  rockTree->Draw("digiHitTHistoID>>hr3t(100,0,4000)");
  rockTree->Draw("digiHitQHistoID>>hr3q(100,0,100)");
  TH1F *hr3t=(TH1F*)gDirectory->Get("hr3t");
  TH1F *hr3q=(TH1F*)gDirectory->Get("hr3q");
cout<<"alright";
  //TH1F *hr3t = new TH1F("hr3t","hr3t",100,1000,2500);
  //TH1F *hr3q = new TH1F("hr3q","hr3q",100,0,50);
  //IDTree->Project("hr3t","digiHitT");cout<<"done 1"<<endl;
  //IDTree->Project("hr3q","digiHitQ");cout<<"done 2"<<endl;
  hr3t->GetXaxis()->SetTitle("time (ns)");
  hr3t->GetYaxis()->SetTitle("Events");
  hr3t->GetXaxis()->SetTitleSize(0.06);
  hr3t->GetYaxis()->SetTitleSize(0.06);
  hr3t->GetYaxis()->SetTitleOffset(0.2);
  hr3q->GetXaxis()->SetTitle("charge (PE)");
  hr3q->GetYaxis()->SetTitle("Events");
  hr3q->GetXaxis()->SetTitleSize(0.06);
  hr3q->GetYaxis()->SetTitleSize(0.06);
  hr3q->GetYaxis()->SetTitleOffset(0.2);
  //hr3t->SetMaximum(fbmx);
  //hr3t->GetXaxis()->SetRange(1,95);
  hr3t->SetLineWidth(2);hr3q->SetLineWidth(2);
  hr3t->SetLineColor(2);hr3q->SetLineColor(2);
  hr3t->SetLineStyle(1);hr3q->SetLineStyle(1);
  //hr3t->SetFillColor(1);hr3q->SetFillColor(2);
  //hr3t->SetFillStyle(1);hr3q->SetFillStyle(3004);
  
  //open Air data
  TFile *ff = new TFile("HistoAirRandomMu10GeV1000Evnts.root");
  TTree *airTree = (TTree*)ff->Get("IDTree");
  airTree->Draw("digiHitTHistoID>>ha3t(100,0,4000)");
  airTree->Draw("digiHitQHistoID>>ha3q(100,0,100)");
  TH1F *ha3t=(TH1F*)gDirectory->Get("hr3t");
  TH1F *ha3q=(TH1F*)gDirectory->Get("hr3q");

  //TH1F *ha3t=(TH1F*) ff->Get("digiHitTHistoID");ha3t->SetName("ha3t");
  //TH1F *ha3q=(TH1F*) ff->Get("digiHitQHistoID");ha3q->SetName("ha3q");
  //TH1F *ha3t = new TH1F("ha3t","ha3t",100,1000,2500);
  //TH1F *ha3q = new TH1F("ha3q","ha3q",100,0,50);
  //IDTree->Project("ha3t","digiHitT");cout<<"done 3"<<endl;
  //IDTree->Project("ha3t","digiHitQ");cout<<"done 4"<<endl;
  ha3t->GetXaxis()->SetTitle("time (ns)");
  ha3t->GetYaxis()->SetTitle("Events");
  ha3t->GetXaxis()->SetTitleSize(0.06);
  ha3t->GetYaxis()->SetTitleSize(0.06);
  ha3q->GetXaxis()->SetTitle("charge (PE)");
  ha3q->GetYaxis()->SetTitle("Events");
  ha3q->GetXaxis()->SetTitleSize(0.06);
  ha3q->GetYaxis()->SetTitleSize(0.06);
  //ha3t->SetMaximum(fbmx);
  //ha3t->GetXaxis()->SetRange(1,95);
  ha3t->SetLineWidth(2);ha3q->SetLineWidth(2);
  ha3t->SetLineColor(4);ha3q->SetLineColor(4);
  ha3t->SetLineStyle(1);ha3q->SetLineStyle(1);
  //ha3t->SetFillColor(1);ha3q->SetFillColor(2);
  //ha3t->SetFillStyle(1);ha3q->SetFillStyle(3004);
  
 //plots
 TCanvas *c101 = new TCanvas("c101","c101",10,10,800,400);
 c101->Divide(1,2);
 c101->cd(1);c101->SetLogy();
 hr3t->GetXaxis()->SetRange(25,100);
 hr3t->Draw();ha3t->Draw("same");
 c101->cd(2);c101->SetLogy();
 hr3q->Draw();ha3q->Draw("same");
 

 TLegend *leg1 = new TLegend(0.6,0.5,0.8,0.8);//left,bottom,right,up
 leg1->AddEntry(hr3t,Form("Rock, %d GeV",nb),"l");
 leg1->AddEntry(ha3t,Form("Air, %d GeV",nb),"l");
 leg1->SetFillColor(0);
 leg1->SetLineColor(0);
 //leg1->SetTextFont(22);
 leg1->SetTextSize(0.05);
 c101->cd(1);leg1->Draw("");
//
// //save histos
// TFile *fsave = new TFile("nuwro.root","RECREATE");
// hwtt->Write();hwqe->Write();hwrs->Write();hwds->Write();


}
