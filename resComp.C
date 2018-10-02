//to compare the results with different cuts
///first set (yX1) is the EB PID, second set with Stefan's PID cuts



void resComp()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit();
  gStyle->SetOptStat(2211);
  float xM[3];
  float xZ[3];
  float xX[3];

  float yM1[3];
  float yM2[3];

  float yZ1[3];
  float yZ2[3];

  float yX1[3];
  float yX2[3];


  float eyM1[3];
  float eyM2[3];

  float eyZ1[3];
  float eyZ2[3];

  float eyX1[3];
  float eyX2[3];


  xM[0]= 0.346893;
  xM[1]=0.610165;
  xM[2]=1.02478;

  xZ[0]=0.28375;
  xZ[1]=0.459325;
  xZ[2]= 0.7111;


  xX[0]=0.152663;
  xX[1]=0.277083;
  xX[2]=0.480122;


  yM1[0]=-0.01472;
  eyM1[0]=0.0148798;

  yM1[1]=-0.013896;
  eyM1[1]=0.00697;

  yM1[2]=0.022354;
  eyM1[2]=0.00617;

  yM2[0]=-0.01078;
  eyM2[0]=0.0175;

  yM2[1]=-0.01904;
  eyM2[1]=0.00813;

  yM2[2]=0.02082;
  eyM2[2]=0.006838;


  yZ1[0]=-0.02416;
  eyZ1[0]=0.0308;

  yZ1[1]=0.001307;
  eyZ1[1]=0.00534;

  yZ1[2]=0.01819;
  eyZ1[2]=0.00771;

  yZ2[0]=-0.0311179;
  eyZ2[0]=0.0400882;

  yZ2[1]=-0.00149162;
  eyZ2[1]=0.00614001;

  yZ2[2]=0.0189443;
  eyZ2[2]=0.00851503;


  yX1[0]=0.0103609;
  eyX1[0]=0.00694;

  yX1[1]=0.002518;
  eyX1[1]=0.006215;

  yX1[2]=0.01036;
  eyX1[2]=0.01266;


  yX2[0]=0.00619881;
  eyX2[0]=0.00789337;

  yX2[1]=0.00194366;
  eyX2[1]=0.00707007;

  yX2[2]=0.0170197;
  eyX2[2]=0.0142986;


  float pullsM_Y[3];
  float pullsM_EY[3];

  float pullsX_Y[3];
  float pullsX_EY[3];

  float pullsZ_Y[3];
  float pullsZ_EY[3];

  for(int i=0;i<3;i++)
    {
      pullsM_Y[i]=2*(yM1[i]-yM2[i])/(eyM1[i]+eyM2[i]);
      pullsZ_Y[i]=2*(yZ1[i]-yZ2[i])/(eyZ1[i]+eyZ2[i]);
      cout <<"pullsZ, i: " << i << " y1: "<< yZ1[i] << " yZ2: " << yZ2[i] <<" diff: "<< yZ1[i]-yZ2[i] <<" uncert: "<< eyZ1[i] <<" and " << eyZ2[i] <<" pull: "<< 2*(yZ1[i]-yZ2[i])/(eyZ1[i]+eyZ2[i])<<endl;
      cout <<"diff factor : "<< (eyZ1[i]+eyZ2[i])<<endl;
      pullsX_Y[i]=2*(yX1[i]-yX2[i])/(eyX1[i]+eyX2[i]);
      cout <<"pullsX, i: " << i << " y1: "<< yX1[i] << " yX2: " << yX2[i] <<" diff: "<< yX1[i]-yX2[i] <<" uncert: "<< eyX1[i] <<" and " << eyX2[i] <<" pull: "<< 2*(yX1[i]-yX2[i])/(eyX1[i]+eyX2[i])<<endl;
      cout <<"diff factor : "<< (eyX1[i]+eyX2[i])<<endl;

      pullsM_EY[i]=2*(eyM1[i]-eyM2[i])/(eyM1[i]+eyM2[i]);
      pullsZ_EY[i]=2*(eyZ1[i]-eyZ2[i])/(eyZ1[i]+eyZ2[i]);
      pullsX_EY[i]=2*(eyX1[i]-eyX2[i])/(eyX1[i]+eyX2[i]);


    }


  TGraph pullsM(3,xM,pullsM_Y);
  TGraph pullsX(3,xX,pullsX_Y);
  TGraph pullsZ(3,xZ,pullsZ_Y);

  TGraph pullsEM(3,xM,pullsM_EY);
  TGraph pullsEX(3,xX,pullsX_EY);
  TGraph pullsEZ(3,xZ,pullsZ_EY);

  TCanvas c;
  pullsM.GetXaxis()->SetTitle("M [GeV]");
  pullsM.GetYaxis()->SetTitle("pulls EB-StefanPID");

  pullsZ.GetXaxis()->SetTitle("z");
  pullsZ.GetYaxis()->SetTitle("pulls between diff PID/Fid cuts");

  pullsX.GetXaxis()->SetTitle("x");
  pullsX.GetYaxis()->SetTitle("pulls EB-StefanPID");


  pullsEM.GetXaxis()->SetTitle("M [GeV]");
  pullsEM.GetYaxis()->SetTitle("uncertainty pulls, diff PID/Fid cuts");

  pullsEZ.GetXaxis()->SetTitle("z");
  pullsEZ.GetYaxis()->SetTitle("uncertainty pulls, diff PID/Fid cuts");

  pullsEX.GetXaxis()->SetTitle("x");
  pullsEX.GetYaxis()->SetTitle("pulls EB-StefanPID");


  pullsM.SetMarkerStyle(24);
  pullsM.GetYaxis()->SetLabelSize(0.04);
  pullsM.GetYaxis()->SetTitleSize(0.04);

  pullsZ.SetMarkerStyle(24);
  pullsZ.GetYaxis()->SetLabelSize(0.04);
  pullsZ.GetYaxis()->SetTitleSize(0.04);

pullsX.SetMarkerStyle(24);
  pullsX.GetYaxis()->SetLabelSize(0.04);
  pullsX.GetYaxis()->SetTitleSize(0.04);

  pullsEM.SetMarkerStyle(24);
  pullsEM.GetYaxis()->SetLabelSize(0.04);
  pullsEM.GetYaxis()->SetTitleSize(0.04);

  pullsEZ.SetMarkerStyle(24);
  pullsEZ.GetYaxis()->SetLabelSize(0.04);
  pullsEZ.GetYaxis()->SetTitleSize(0.04);

  pullsEX.SetMarkerStyle(24);
  pullsEX.GetYaxis()->SetLabelSize(0.04);
  pullsEX.GetYaxis()->SetTitleSize(0.04);

  pullsM.Draw("AP");
  c.SaveAs("pullsM.png");
  c.SaveAs("pullsM.pdf");
  pullsZ.Draw("AP");
  c.SaveAs("pullsZ.png");
  c.SaveAs("pullsZ.pdf");
  pullsX.Draw("AP");
  c.SaveAs("pullsX.png");
  c.SaveAs("pullsX.pdf");

  pullsEM.Draw("AP");
  c.SaveAs("pullsEM.png");
  c.SaveAs("pullsEM.pdf");
  pullsZ.Draw("AP");
  c.SaveAs("pullsEZ.png");
  c.SaveAs("pullsEZ.pdf");
  pullsEX.Draw("AP");
  c.SaveAs("pullsEX.png");
  c.SaveAs("pullsEX.pdf");


}
