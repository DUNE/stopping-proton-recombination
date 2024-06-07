#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>

TFile *braggfile = new TFile("bragg_data/bragg_curves.root");
TF1   *braggpeak_rev=(TF1 *)braggfile->Get("fn");
TF1   *braggpeak_rev_RS=(TF1 *)braggfile->Get("fn_RS");

struct FitResult
{
  double trackshift;
  double scaling;
  double chi2;
  int nparams;
};

//Raw PSTAR Bragg curve for Ar
Double_t floating_bragg (Double_t *x, Double_t *par)
{
  return 100*par[0]*braggpeak_rev->Eval(*x + par[1]);
}

//Same curve resampled to 6mm intervals
Double_t floating_bragg_RS (Double_t *x, Double_t *par)
{
  return 100.*par[0]*braggpeak_rev_RS->Eval(*x + par[1]);
}

FitResult ExtractTrackShift(TGraph* &dEdx_gr, double min_rr)
{
  //Extract lower limit to fit
  /*
  double max_dedx = 0;
  double max_dedx_index = 99;
  for (int i = 0; i < dEdx_gr->GetN(); i++)
  {
    double dedx = dEdx_gr->GetPointY(i);
    if (dedx > max_dedx)
    {
      max_dedx = dedx;
      max_dedx_index = i;
    }
  }
  double min_rr = dEdx_gr->GetPointX(max_dedx_index+1);
  */

  //Define function to fit to
  TF1 *fl_br;
  double max_rr = 45; //The residual range the APS data extends to
  fl_br = new TF1("fl_br", floating_bragg, 0, 42, 2);
  //fl_br = new TF1("fl_br", floating_bragg_RS, 0, 42, 2);

  //Configure fitter
  fl_br->SetParameter(0, 0.6);
  fl_br->SetParError(0, 0.01);

  fl_br->SetParameter(1, 0.);
  fl_br->SetParError(1, 0.01);

  //Perform fit 
  dEdx_gr->Fit("fl_br", "WQN", "SAME", min_rr, (max_rr - 3.));

  FitResult result;
  if (fl_br != NULL)
  {
    result.trackshift = fl_br->GetParameter(1);
    result.scaling    = fl_br->GetParameter(0);
    result.chi2       = fl_br->GetChisquare();
    result.nparams    = fl_br->GetNpar();
  }
  else
  {
    result.trackshift = -9999;
    result.scaling    = -9999;
    result.chi2       = -9999;
    result.nparams    = -9999;
  }
  return result;
}

