
#include "TROOT.h"
#include "TImage.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TF2.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TColor.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <vector>
#include <string>
#include <sstream>

using namespace std;


//questionable runs: 4308, 4302 4032 4013, 4016, 4017 4068.82, 4073.62

void doFit(float*** vals,float** wyFactors, float& amp1, float& amp1Err,float& amp2, float& amp2Err,int numAngBins,float r, TH1D* hChi2)
{

  cout <<"r: "<< r <<endl;
  TGraph2DErrors g;
  for(int iAngBin=0;iAngBin<numAngBins;iAngBin++)
    {
      for(int iAngBin2=0;iAngBin2<numAngBins;iAngBin2++)
	{
	  int pointIndex=iAngBin*numAngBins+iAngBin2;
	  double N1 = vals[iAngBin][iAngBin2][0];
	  double N2 = vals[iAngBin][iAngBin2][1];
	  float y=0.0;
	  float ey=0.0;
	  float x1=0.0;
	  float x2=0.0;
	  float ex=0.0;
	  if ((N1 + r * N2) > 0) {
	    y = (N1 - r * N2) / (N1 + r * N2);
	    ey = N1 * 4 * r * r * N2 * N2 / ((N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2));
	    ey += N2 * 4 * r * r * N1 * N1
	      / ((N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2));
	    ey = sqrt(ey);


	    ///phi dependent wy factor:
	    cout <<"correcting with wy factors: "<< wyFactors[iAngBin][iAngBin2] <<" for ang bin " << iAngBin<<", " << iAngBin2 <<endl;
	    y/=wyFactors[iAngBin][iAngBin2];
	    ey/=wyFactors[iAngBin][iAngBin2];
	    //

	  } else {
	    //	    cout <<"no counts for phi bin " << iAngBin << ", " << iAngBin2 << endl;
	    y = 0;
	    //no counts... uncertainty high
	    ey = 10000;
	  }
	  x1 = (iAngBin + 0.5) * 2 * M_PI / numAngBins;
	  x2 = (iAngBin2 + 0.5) * 2 * M_PI / numAngBins;

	  g.SetPoint(pointIndex,x1,x2,y);
	  g.SetPointError(pointIndex,ex,ex,ey);
	}
    }
  TF2 f2("f2","[0]*sin(x)+[1]*sin(y-x)",0,2*M_PI,0,2*M_PI);
  f2.SetParameters(0.0,0.0);
  g.Fit(&f2);
  hChi2->Fill(f2.GetChisquare()/f2.GetNDF());

  cout <<"chi2/ndf: "<< f2.GetChisquare()/f2.GetNDF() <<" chi2: "<< f2.GetChisquare()<<" ndf: "<< f2.GetNDF()<<endl;
  amp1=f2.GetParameter(0);
  amp2=f2.GetParameter(1);

  amp1Err=f2.GetParError(0);
  amp2Err=f2.GetParError(1);

}


int main(int argc, char** argv)
{

  float yMin[]={-0.035,-0.015,-0.01};
  //  float yMax[]={0.05,0.055,0.05};
  //with m cut
  float yMax[]={0.05,0.1,0.18};



  float beamPolarization=0.86;
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit();
  gStyle->SetOptStat(2211);

  gStyle->SetOptTitle(0);

  //lets assume that 16 is maxAngBins,
  int maxAngBins=16;
  int numOutbendingRuns=68;
  int maxKinBins=50+numOutbendingRuns;
  float* meanKin=new float[maxKinBins];
  float* wyFactor=new float[maxKinBins];
   float** wyFactor2D=new float*[maxAngBins];
  float*** vals=new float**[maxAngBins];
  for(int i=0;i<maxAngBins;i++)
    {
      vals[i]=new float*[maxAngBins];
      wyFactor2D[i]=new float[maxAngBins];
      for(int j=0;j<maxAngBins;j++)
	{
	  vals[i][j]=new float[2];
	  for(int k=0;k<2;k++)
	    {
	      vals[i][j][k]=0.0;
	    }
	}
    }

  if(argc<2)
    {
      cout <<"need input file! "<<endl;
      exit(0);
    }
  cout <<"opening " << argv[1]<<endl;
  ifstream file;
  file.open(argv[1]);
  string line;
  int binIndex=0;  
  ///this tokenizes lines
  TH1D* hChi2=new TH1D("hChi2","hChi2",20,0,5);
  hChi2->GetXaxis()->SetTitle("#chi^{2}/NDF");
  //  hChi2->GetXaxis()->SetTitleSize(20);
  while(getline(file,line))
    {
      //and this space
      stringstream ssBinning(line);
      string binningName,sNumKinBins,sNumAngBins;

      string sKin;
      string sWy;
      float kin=0;
      float wy=0;
      getline(ssBinning,binningName,' ');
      getline(ssBinning,sNumKinBins,' ');
      getline(ssBinning,sNumAngBins,' ');
      int numKinBins=stoi(sNumKinBins);
      int numAngBins=stoi(sNumAngBins);

      for(int iKinBin=0;iKinBin<numKinBins;iKinBin++)
	{
	  ssBinning>>sKin;
	  kin=stof(sKin);
	  if(isnan(kin))
	    kin=-1;
	  meanKin[iKinBin]=kin;
	  ssBinning>>sWy;
	  wy=stof(sWy);
	  if(isnan(wy))
	    wy=-1000;
	  wyFactor[iKinBin]=wy;
	  cout <<" mean " << binningName <<" bin " << iKinBin <<": " << kin << endl;
	  cout <<" wyfactor " << binningName <<" bin " << iKinBin <<": " << wy << endl;
	}

      cout <<"binningName: "<< binningName<<", " << numKinBins<<", " << numAngBins<<endl;
      binIndex++;
      ///next line
      getline(file,line);
      stringstream ssVals(line);
      string sVal;
      string srVal;
      float val;
      float rVal;

      float y1[1000];
      float y2[1000];
      float x[1000];

      float ey1[1000];
      float ey2[1000];
      float ex[1000];
      //in case some kin bins are empty
      int graphIndex=0;
      for(int iKinBin=0;iKinBin<numKinBins;iKinBin++)
	{
	  ssVals>>srVal;
	  rVal=stof(srVal);
	  cout <<" r val is : "<<rVal<< endl;
	  for(int iAngBin=0;iAngBin<numAngBins;iAngBin++)
	    {
	      for(int iAngBin2=0;iAngBin2<numAngBins;iAngBin2++)
		{
		  //first wy 2D
		  ssVals>>srVal;
		  float wy2D=stof(srVal);
		  wyFactor2D[iAngBin][iAngBin2]=wy2D;
		  //cout <<"wy2D: "<< wy2D <<endl;
		  for(int ipol=0;ipol<2;ipol++)
		    {
		      ssVals >> sVal;
		      val=stof(sVal);
		      vals[iAngBin][iAngBin2][ipol]=val;
		    }
		}
	    }
	  ////do fit
	  float amp1=0.0;
	  float amp1Err=0.0;

	  float amp2=0.0;
	  float amp2Err=0.0;

	  if(meanKin[iKinBin]!=-1)
	    {
	      cout <<"do fit for binIndex: "<< binIndex<<" mean kin: "<< meanKin[iKinBin]<<endl;
	      doFit(vals,wyFactor2D,amp1,amp1Err,amp2,amp2Err,numAngBins,rVal,hChi2);
	      //already pointby poin
	      //	      y1[graphIndex ]=amp1/wyFactor[iKinBin];
	      //	      ey1[graphIndex]=amp1Err/wyFactor[iKinBin];
	      y1[graphIndex]=amp1;
	      ey1[graphIndex]=amp1Err;
	      y2[graphIndex]=amp2;
	      ey2[graphIndex]=amp2Err;

	      //and correct for beam polarization
	      y1[graphIndex]/=beamPolarization;
	      ey1[graphIndex]/=beamPolarization;
	      y2[graphIndex]/=beamPolarization;
	      ey2[graphIndex]/=beamPolarization;

	      x[graphIndex]=meanKin[iKinBin];
	      ex[graphIndex]=0.0;
	      graphIndex++;

	      cout <<"index : "<< graphIndex-1<<" y1: "<< y1[graphIndex-1];
	      cout <<" y2: "<< y2[graphIndex-1] <<", ey1: "<< ey1[graphIndex-1];
	      cout <<", ey2: "<< ey2[graphIndex-1] <<" x: " << x[graphIndex-1];
	      cout<<" ex: "<< ex[graphIndex-1] <<endl;
	    }
	  else
	    {
	      cout <<"skipping kinbin: "<< iKinBin<<endl;
	    }

	  ///////

	}
      TCanvas c("can","can",1000,800);
      c.cd();
      cout <<"using " << graphIndex<< " points for graph " <<endl;
      for(int i=0;i<graphIndex;i++)
	{
	  cout <<"x: "<< x[i] <<" y1: "<< y1[i] <<" y2: "<< y2[i] <<endl;
	}
      TGraphErrors g1(graphIndex,x,y1,ex,ey1);
      TGraphErrors g2(graphIndex,x,y2,ex,ey2);
      string xaxisLabel;
	  g1.GetYaxis()->SetRangeUser(yMin[binIndex-1],yMax[binIndex-1]);
      if(binIndex==1)
	{
	  xaxisLabel="M_{Inv} [GeV/c^{2}]";

	}
      if(binIndex==2)
	{
	  xaxisLabel="z";

	}
      if(binIndex==3)
	xaxisLabel="x";
      if(binIndex==4)
	{
	  xaxisLabel="run number";
		g1.GetYaxis()->SetRangeUser(-0.2,0.2);
		g2.GetYaxis()->SetRangeUser(-0.2,0.2);
	}


      ///from example

      //	g1->GetXaxis()->SetTitle("#phi (^{o})");
      //	g1->GetXaxis()->SetNdivisions(206,kFALSE);
      //	g1->GetXaxis()->SetLabelSize(0.05);
      //	g1->GetYaxis()->SetLabelSize(0.05);
      //	g1->GetXaxis()->SetTitleSize(0.05);
      //	g1->GetXaxis()->SetTitleOffset(0.9);

      ////



	g1.SetMarkerStyle(24);
	g2.SetMarkerStyle(24);
	g1.GetYaxis()->SetLabelSize(0.04);
	g1.GetYaxis()->SetTitleSize(0.06);
	g1.GetYaxis()->SetTitle("A_{LU}^{sin(#phi_{R})}");
	g2.GetYaxis()->SetTitleSize(0.06);
	g2.GetYaxis()->SetLabelSize(0.04);
	g2.GetYaxis()->SetTitle("A_{LU}^{sin(#phi_{H}-#phi_{R})}");



	g1.GetXaxis()->SetLabelSize(0.04);
	g1.GetXaxis()->SetTitle(xaxisLabel.c_str());
	g2.GetXaxis()->SetLabelSize(0.04);
	g2.GetXaxis()->SetTitle(xaxisLabel.c_str());

	float meanAmp=0.0;
	if(binIndex==4)
	  {
	    //since we do it only for runs, but large number here...
	    TF1 f1("f1","[0]",x[0]-100,x[graphIndex]+100);
	    f1.SetParameter(0,0.0);
	    g1.Fit(&f1);
	    //this is the first fit to the sin(phiR)!!
	    meanAmp=f1.GetParameter(0);
	    f1.SetParameter(0,0.0);
	    g2.Fit(&f1);



	  }
	//	      c.Divide(2,1);
	//	      TVirtualPad* pad1=c.cd(1);
	//	      pad1->SetLeftMargin(0.2);
	///////////
	//	TH2F *H = new TH2F("H","Example BSA",100,0,g1.,100,-0.32,0.32);
	g1.GetYaxis()->SetTitle("A_{LU}^{sin(#phi_{R})}");
	g1.GetXaxis()->SetTitle(xaxisLabel.c_str());
	//	H->GetXaxis()->SetTitle("#phi (^{o})");
	//	g1.GetXaxis()->SetNdivisions(206,kFALSE);
	g1.GetXaxis()->SetLabelSize(0.05);
	g1.GetYaxis()->SetLabelSize(0.05);
	g1.GetXaxis()->SetTitleSize(0.05);
	g1.GetXaxis()->SetTitleOffset(0.9);

	TImage *i1 = TImage::Open("clasPic.png");
	float d = 0.23;
	float e = 0.89;
	float ratio = 1.0*i1->GetWidth()/(1.0*i1->GetHeight());
	cout << ratio << endl;


	g1.Draw("AP");
	//	TVirtualPad* cPad=gROOT->GetSelectedPad();
	TPad *p1 = new TPad("i1", "i1",e - d * ratio , e - d, e, e);
	p1->Draw();
	p1->cd();
	i1->Draw("SAME");
	g1.Draw("PSAME");

	//	H->Draw();
	/////////
	g1.SetMarkerStyle(20);
	g1.SetLineWidth(2);
	c.cd();



	//      g1.Draw("Psame");
      gPad->Update();
      TLine l(gPad->GetUxmin(),0.0,gPad->GetUxmax(),0.0);
      l.Draw();

      TText *t;
      TPaveText* prelim;
      TPaveText* pavetitle;
      TPaveText* fitres;


      TLatex reaction;


      //already incremented
      //m,z,x
      if(binIndex==1)
	{
	  reaction.DrawLatex(0.33,0.04,"e p#rightarrow e' #pi^{+} #pi^{-} +X");
	  t=new TText(0.45,0.02,"CLAS Preliminary");
	  prelim = new TPaveText(0.5,0.02,1.0,0.1);
	  pavetitle = new TPaveText(30,0.33,330,0.39);
	  fitres = new TPaveText(20,-0.3,180,-0.05);
	}
      if(binIndex==2)
	{
	  reaction.DrawLatex(0.33,0.04,"e p#rightarrow e' #pi^{+} #pi^{-} +X");
	  t=new TText(0.4,0.01,"CLAS Preliminary");
	  prelim = new TPaveText(0.4,0.01,0.8,0.1);
	  pavetitle = new TPaveText(30,0.33,330,0.39);
	  fitres = new TPaveText(20,-0.3,180,-0.05);
	}
      if(binIndex==3)
	{
	  reaction.DrawLatex(0.15,0.04,"e p#rightarrow e' #pi^{+} #pi^{-} +X");
	  t=new TText(0.2,0.03,"CLAS Preliminary");
	  prelim = new TPaveText(0.2,0.012,0.5,0.03);
	  pavetitle = new TPaveText(30,0.33,330,0.39);
	  fitres = new TPaveText(20,-0.3,180,-0.05);
	}


	prelim->SetFillStyle(0);
	prelim->SetLineWidth(0);
	prelim->SetBorderSize(0);
	prelim->SetTextColorAlpha(1,0.2);
	prelim->AddText("PRELIMINARY");
	prelim->Draw();

	pavetitle->SetLineWidth(0);
	pavetitle->SetBorderSize(0);
	pavetitle->SetFillStyle(0);
	pavetitle->AddText("Example BSA");
	pavetitle->Draw();

	fitres->SetBorderSize(0);
	fitres->SetFillStyle(0);
	fitres->AddText("x_{B} = 0.5");
	fitres->AddText("Q^{2} = 8.5 GeV^{2}");
	fitres->AddText("-t = 0.3 GeV^{2}");





	//	fitres->AddText(Form("A^{sin#phi} = %1.2f #pm %1.3f",fit->GetParameter(0),fit->GetParError(0)));
	//	fitres->AddText(Form("A^{sin2#phi} = %1.2f #pm %1.3f",fit->GetParameter(1),fit->GetParError(1)));
	//	fitres->AddText(Form("#chi^{2} / NDF = %1.2f / %d",fit->GetChisquare(),fit->GetNDF()));
	fitres->Draw();

	//      t->SetTextFont(43);
	//      t->SetTextSize(30);
      ///      t->Draw();
      gPad->Update();
      //            pad1=c.cd(2);
      //            pad1->SetLeftMargin(0.3);
      //          g2.Draw("AP");
      char buffer[300];

      sprintf(buffer,"asym2DFit_out_%s.png",binningName.c_str());
      c.SaveAs(buffer);
      sprintf(buffer,"asym2DFit_out_%s.pdf",binningName.c_str());
      c.SaveAs(buffer);
      sprintf(buffer,"asym2DFit_out_%s.root",binningName.c_str());
      c.SaveAs(buffer);
      //do some more studies of the run dependendence
	if(binIndex==4)
	  {
	    cout << "mean Amp: "<< meanAmp<<endl;
	    TH1D pulls1D("pulls1D","pulls1D",20,-3.5,3.5);
	    for(int i=0;i<graphIndex;i++)
	      {
		if(ey1[i]!=0)
		  {
		    y1[i]=(y1[i]-meanAmp)/ey1[i];
		    pulls1D.Fill(y1[i]);
		  }
		cout <<"run: " << x[i] <<" pull: "<< y1[i] <<endl;
	      }
	    //	    pulls1D.Fit("gaus");
	    //	    c.Divide(1,1);
	    c.cd(0);
	    TGraph g3(graphIndex,x,y1);
	    g3.GetXaxis()->SetTitle("run number");
	    g3.GetYaxis()->SetTitle("pull");
	    g3.SetMarkerStyle(20);
	    g3.SetMarkerSize(1);
	    g3.Draw("AP");

	    c.SaveAs("pullsRuns.png");
	    c.SaveAs("pullsRuns.pdf");
	    pulls1D.GetXaxis()->SetTitle("pull");
	    pulls1D.Draw();
	    c.SaveAs("pullsRuns1D.png");
	    c.SaveAs("pullsRuns1D.pdf");



	  }

    }


  TCanvas c2;
  hChi2->GetXaxis()->SetTitle("fit #chi^{2}");
  hChi2->Draw();
  c2.SaveAs("fitChi2.png");
  c2.SaveAs("fitChi2.pdf");


  file.close();
  cout <<"Hello World!" <<endl;
  return 0;
}
