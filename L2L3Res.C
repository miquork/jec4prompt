// Purpose: Read input JES data as TProfile2D and produce JEC L2L3Res fit
//          Focus on JEC stability and robustness
#include "TFile.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"

void L2L3Res() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/L2L3Res");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/L2L3Res");

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  ////////////////////////////
  // Load input data and MC //
  ////////////////////////////
  
  // Load input file
  int run = 392175;
  TFile *f = new TFile("rootfiles/J4PHists_runs392175to392175_photonjet.root","READ");
  assert(f && !f->IsZombie());
  TFile *fm = new TFile("rootfiles/reweighted_J4PHists_photonjet_GJ-4Jets.root","READ");
  assert(fm && !fm->IsZombie());
  TFile *fmm = new TFile("../jecsys3/rootfiles/Prompt2024/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root","READ");
  assert(fmm && !fmm->IsZombie());
  
  
  curdir->cd();

  // Load input data (MPF, DB) from file
  //TProfile2D *p2m0 = (TProfile2D*)f->Get("MPF_2D"); assert(p2m0);
  TProfile2D *p2m0 = (TProfile2D*)f->Get("DB_2D"); assert(p2m0);
  TProfile2D *p2m0m = (TProfile2D*)fm->Get("DB_2D"); assert(p2m0m);
  TProfile2D *p2m0mm = (TProfile2D*)fmm->Get("Gamjet2/p2m2"); assert(p2m0mm);
  TProfile2D *p2corrmm = (TProfile2D*)fmm->Get("Gamjet2/p2corr"); assert(p2corrmm);
  
  
  ////////////////////////////////////////////////////////////
  // Initial plot of data and MC with rebinning and cutting //
  ////////////////////////////////////////////////////////////

  // Define pT rebinning to ensure sufficient statistics
  // Original binning
  // {5, 7, 9, 11, 13, 15, 17, 20, 24, 28, 32, 36, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000}
  const Double_t v13[] =
    {20, 28, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153,
     174, 220, 300, 430, 638, 1032, 2000, 7000};
  const int n13 = sizeof(v13)/sizeof(v13[0])-1;
  
  // Find bin edges in eta
  int i1 = p2m0->GetXaxis()->FindBin(-1.305);
  int i2 = p2m0->GetXaxis()->FindBin(+1.305)-1;
  double eta1 = p2m0->GetXaxis()->GetBinLowEdge(i1);
  double eta2 = p2m0->GetXaxis()->GetBinLowEdge(i2)+1;

  cout << Form("Fitting %1.3f #LT eta < %1.3f reference region\n",eta1,eta2);
  
  // Offline bin edges are different than JEC4Prompt, so redo
  int i1mm = p2m0mm->GetXaxis()->FindBin(-1.305);
  int i2mm = p2m0mm->GetXaxis()->FindBin(+1.305)-1;
  double eta1mm = p2m0mm->GetXaxis()->GetBinLowEdge(i1mm);
  double eta2mm = p2m0mm->GetXaxis()->GetBinLowEdge(i2mm)+1;

  // Slice out |eta|<1.305 and rebin
  TProfile *p1m0 = p2m0->ProfileY("p1m0",i1,i2);
  TProfile *p1m0_rebin = (TProfile*)p1m0->Rebin(n13,"p1m0_rebinned",v13);
  TH1D *h1m0 = p1m0_rebin->ProjectionX("h1m0");
  TH1D *h1m0_cut = (TH1D*)h1m0->Clone("h1m0_cut");
  h1m0_cut->GetXaxis()->SetRangeUser(40,300);

  TProfile *p1m0m = p2m0m->ProfileY("p1m0",i1,i2);
  TProfile *p1m0m_rebin = (TProfile*)p1m0m->Rebin(n13,"p1m0m_rebinned",v13);
  TH1D *h1m0m = p1m0m_rebin->ProjectionX("h1m0m");
  TH1D *h1m0m_cut = (TH1D*)h1m0m->Clone("h1m0m_cut");
  h1m0m_cut->GetXaxis()->SetRangeUser(40,300);

  // Don't rebin offline MC because bin edges don't match with JEC4Prompt
  TProfile *p1m0mm = p2m0mm->ProfileY("p1m0mm",i1mm,i2mm);
  TH1D *h1m0mm = p1m0mm->ProjectionX("h1m0mm");

  // Scale out previous L2Relative in offline MC to compare to JEC4Prompt
  TProfile *p1corrmm = p2corrmm->ProfileY("pcorrmm",i1mm,i2mm);
  TH1D *h1corrmm = p1corrmm->ProjectionX("h1corrmm");
  TH1D *h1m0mm_corr = (TH1D*)h1m0mm->Clone("h1m0mm_corr");
  h1m0mm_corr->Divide(h1corrmm);

  // Create canvas for nice background
  double xmin(15), xmax(4500);
  //TH1D *h = tdrHist("h","JES with MPF",0.82,1.18,"Tag photon p_{T,#gamma} (GeV)",15,4500);
  TH1D *h = tdrHist("h","JES with DB",0.62,1.38,
		    "Tag photon p_{T,#gamma} (GeV)",15,4500);
  lumi_136TeV = Form("Run %d, 1 fb^{-1}",run);
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();
  drawCustomLogXLabels(h);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(xmin,1,xmax,1);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.18, 0.75,"|#eta| < 1.305");

  tdrDraw(p1m0mm,"HIST",kNone,kBlue-9,kSolid,-1,kNone,0);
  tdrDraw(h1m0mm_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p1m0m_rebin,"HIST",kNone,kGreen+2,kSolid,-1,kNone,0);
  tdrDraw(p1m0,"Pz",kOpenSquare,kRed-9,kSolid,-1,kNone,0);
  tdrDraw(p1m0_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(h1m0_cut,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

  TLegend *leg = tdrLeg(0.40,0.90-0.05*6,0.65,0.90);
  leg->AddEntry(h1m0_cut,"Cut and rebinned data","PLE");
  leg->AddEntry(p1m0_rebin,"Rebinned data","PLE");
  leg->AddEntry(p1m0,"Original data","PLE");
  leg->AddEntry(p1m0m_rebin,"Summer24 MC (online)","PLE");
  leg->AddEntry(h1m0mm_corr,"Summer24 MC (offline)","PLE");
  leg->AddEntry(p1m0mm,"Summer24 MC with JES","PLE");

  gPad->RedrawAxis();
    
  c1->SaveAs("pdf/L2L3Res_c1_Eta13.pdf");

  
  ////////////////////////////////////////////
  // Do data/MC ratio for fitting residuals //
  ////////////////////////////////////////////

  TH1D *h1m0r = (TH1D*)h1m0_cut->Clone("h1m0r");
  //h1m0r->Divide(h1m0m);
  // Do ratio manually due to different binning
  for (int i = 1; i != h1m0r->GetNbinsX()+1; ++i) {
    double pt = h1m0r->GetXaxis()->GetBinCenter(i);
    int j = h1m0m->GetXaxis()->FindBin(pt);
    double ymc = h1m0m->GetBinContent(j);
    //double ymc = h1m0mm_corr->Interpolate(pt);
    h1m0r->SetBinContent(i, ymc>0 ? h1m0->GetBinContent(i)/ymc : 0);
    h1m0r->SetBinError(i, ymc>0 ? h1m0->GetBinError(i)/ymc : 0);
  } // for i
  
  //TH1D *h2 = tdrHist("h2","JES with MPF over MC",0.92,1.08,
  TH1D *h2 = tdrHist("h2","Residual JES with DB over MC",0.92,1.08,
		     "Tag photon p_{T,#gamma} (GeV)",15,4500);
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);
  gPad->SetLogx();
  drawCustomLogXLabels(h2);

  l->DrawLine(xmin,1,xmax,1);
  tex->DrawLatex(0.18, 0.75,"|#eta| < 1.305");

  tdrDraw(h1m0r,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

  c2->SaveAs("pdf/L2L3Res_c2_Eta13_Ratio.pdf");
  
  
  // Loop over |eta| bins, rebinning data to keep uncertainties controlled
  
  
  // Monitoring plots for data-fit


  // Monitoring plots for fit chi2 vs |eta|


  // Production of L2L3Res text file
  // (Later, add also mapping from pTtag to <pTprobe,raw> before this)
  
} // L2L3Res
