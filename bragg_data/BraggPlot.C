#include <iostream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include <sstream>
#include <string>

TGraph *gr;
TGraph *gr_csda;

double func(Double_t *x, Double_t *) 
{
  return gr->Eval(x[0]);
}

double func_csda(Double_t *x, Double_t *) 
{
  return gr_csda->Eval(x[0]);
}


void BraggPlot()
{
  ifstream infile;
  infile.open("apdata_lookup.txt");

  double E_proj, E_csda, range;
  double LAr_density = 1.784E-3; //density in g/cm^3

  int npoints = 132;
  double energy[132], energy_csda[132], resrange[132], dEdx[132], dEdx_csda[132];

  int counter = 0;
  while (infile >> E_proj >> E_csda >> range)
  {
    energy[counter] = E_proj;
    energy_csda[counter] = E_csda;
    resrange[counter] = range/LAr_density;

    if (counter > 0)
    {
      dEdx[counter] = (energy[counter] - energy[counter-1])/(resrange[counter] - resrange[counter-1]);
      dEdx_csda[counter] = (energy_csda[counter] - energy_csda[counter-1])/(resrange[counter] - resrange[counter-1]);
    }
    else
    {
      dEdx[counter] = energy[counter]/resrange[counter];
      dEdx_csda[counter] = energy_csda[counter]/resrange[counter];
    }
 
    
    counter++;
  }

  //TGraph *gr = new TGraph(70, resrange, dEdx);
  
  gr = new TGraph(70, resrange, dEdx);
  gr->SetLineColor(kBlue);
  gr->SetMarkerColor(kBlue);
  gr_csda = new TGraph(70, resrange, dEdx_csda);
  gr_csda->SetLineColor(kRed);
  gr_csda->SetMarkerColor(kRed);

  double xmin = 0;
  double xmax = 45;
  TF1 *fn = new TF1("fn", func, xmin, xmax, 0);
  TF1 *fn_csda = new TF1("fn_csda", func_csda, xmin, xmax, 0); 
 
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_csda);

  TCanvas *test = new TCanvas();
  test->SetLogy();
  mg->Draw("AC");
  
  test->Update();

  TFile *fout = new TFile("bragg_curves.root", "RECREATE");
  gr->Write();
  gr_csda->Write();
  fn->Write();
  fn_csda->Write();
  fout->Close();
 
}

