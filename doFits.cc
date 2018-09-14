#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TF2.h"
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


void doFit(float*** vals,float& amp1, float& amp1Err,float& amp2, float& amp2Err,int numAngBins,float r)
{


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
  amp1=f2.GetParameter(0);
  amp2=f2.GetParameter(1);

  amp1Err=f2.GetParError(0);
  amp2Err=f2.GetParError(1);

}


int main(int argc, char** argv)
{
  //lets assume that 16 is maxAngBins,
  int maxAngBins=16;
  int maxKinBins=50;
  float* meanKin=new float[maxKinBins];
  float*** vals=new float**[maxAngBins];
  for(int i=0;i<maxAngBins;i++)
    {
      vals[i]=new float*[maxAngBins];
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
  
  ///this tokenizes lines
  while(getline(file,line))
    {
      //and this space
      stringstream ssBinning(line);
      string binningName,sNumKinBins,sNumAngBins;
      string sKin;
      float kin=0;
      getline(ssBinning,binningName,' ');
      getline(ssBinning,sNumKinBins,' ');
      getline(ssBinning,sNumAngBins,' ');
      int numKinBins=stoi(sNumKinBins);
      int numAngBins=stoi(sNumAngBins);

      for(int iKinBin=0;iKinBin<numKinBins;iKinBin++)
	{
	  ssBinning>>sKin;
	  kin=stof(sKin);
	  meanKin[iKinBin]=kin;
	  cout <<" mean " << binningName <<" bin " << iKinBin <<": " << kin << endl;
	}

      cout <<"binningName: "<< binningName<<", " << numKinBins<<", " << numAngBins<<endl;
      ///next line
      getline(file,line);
      stringstream ssVals(line);
      string sVal;
      string srVal;
      float val;
      float rVal;

      float y1[10];
      float y2[10];
      float x[10];

      float ey1[10];
      float ey2[10];
      float ex[10];
      for(int iKinBin=0;iKinBin<numKinBins;iKinBin++)
	{
	  ssVals>>srVal;
	  rVal=stof(srVal);
	  cout <<" r val is : "<<rVal<< endl;
	  for(int iAngBin=0;iAngBin<numAngBins;iAngBin++)
	    {
	      for(int iAngBin2=0;iAngBin2<numAngBins;iAngBin2++)
		{
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
	  doFit(vals,amp1,amp1Err,amp2,amp2Err,numAngBins,rVal);
	  y1[iKinBin]=amp1;
	  ey1[iKinBin]=amp1Err;
	  y2[iKinBin]=amp2;
	  ey2[iKinBin]=amp2Err;
	  x[iKinBin]=meanKin[iKinBin];
	  ex[iKinBin]=0.0;
	  ///////

	}
      TCanvas c;
      TGraphErrors g1(numKinBins,x,y1,ex,ey1);
      TGraphErrors g2(numKinBins,x,y2,ex,ey2);
      c.Divide(2,1);
      c.cd(1);
      g1.Draw("ALP");
      c.cd(2);
      g2.Draw("ALP");
      char buffer[300];
      sprintf(buffer,"asym2DFit_out_%s.png",binningName.c_str());
      c.SaveAs(buffer);

    }


  file.close();
  cout <<"Hello World!" <<endl;
  return 0;
}
