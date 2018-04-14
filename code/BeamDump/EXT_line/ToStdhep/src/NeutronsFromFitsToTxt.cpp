#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <cstdlib>

#include <chrono>
#include <random>


using namespace std;

  double const X_low  = -2000;
  double const X_high = 6000;
  double const Y_low  = -2000;
  double const Y_high = 2000;
  double const CX_low1  =-0.01; 
  double const CX_high1 =0.022;
  double const CX_low2  = -1; 
  double const CX_high2 =  1;
  double const CY_low1  =-0.006;
  double const CY_high1 = 0.006;
  double const CY_low2  = -1;
  double const CY_high2 =  1;
  double const CZ_low1  =-0.99999;
  double const CZ_high1 =-0.99998;
  double const CZ_low2  = -1;
  double const CZ_high2 =  -0.9995;
  double const CZ_low3  = -0.9995;
  double const CZ_high3 =  0;
  double const E_low  = 0;
  //double const E_high = 0.006;
  double const E_high = 0.06;
  double const T_low  = 0;
  double const T_high = 0.000012;
  //double const T_high = 0.0012;
  TF1 fit_X ("X", "pol5"  ,X_low,  X_high);
  TF1 fit_Y ("Y", "pol4"  ,Y_low,  Y_high);
  TF1 fit_CX1("CX1","pol6"  ,CX_low1, CX_high1);
  TF1 fit_CX2("CX2","gaus(0) + gaus(3)"  ,CX_low2, CX_high2);
  TF1 fit_CY1("CY1","pol0"  ,CY_low1, CY_high1);
  TF1 fit_CY2("CY2","gaus(0) + gaus(3)",CY_low2, CY_high2);
  TF1 fit_CZ1("CZ1","pol1",CZ_low1, CZ_high1);
  TF1 fit_CZ2("CZ2","expo",CZ_low2, CZ_high2);
  TF1 fit_CZ3("CZ3","expo",CZ_low3, CZ_high3);
  TF1 fit_E ("E", "pol1"  ,E_low,  E_high);
  TF1 fit_T ("T", "gaus(0) + gaus(3)"  ,T_low,  T_high);

  double const X_ymax = 1.5;
  double const Y_ymax = 1.8;
  double const CX_ymax1 = 4;
  double const CX_ymax2 = 0.4;
  double const CY_ymax1 = 6;
  double const CY_ymax2 = 0.4;
  double const CZ_ymax1 = 3;
  double const CZ_ymax2 = 1;
  double const CZ_ymax3 = 0.6;
  double const E_ymax = 0.6;
  double const T_ymax = 35;

  TH1D histo_X("histo_X","X",100,X_low,X_high);
  TH1D histo_Y("histo_Y","Y",100,Y_low,Y_high);
  TH1D histo_CX("histo_CX","CX",100,-1,1);
  TH1D histo_CY("histo_CY","CY",100,-1,1);
  TH1D histo_CZ("histo_CZ","CZ",100,-1,0);
  TH1D histo_E("histo_E","E",100,E_low,E_high);
  TH1D histo_T("histo_T","T",100,T_low,T_high);
  FILE * pfile = NULL;

int main()
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::minstd_rand0 generator (seed);

  fit_X.SetParameter(0,    0.223719);
  fit_X.SetParameter(1,-3.03467e-05);
  fit_X.SetParameter(2, 1.24832e-07);
  fit_X.SetParameter(3, 2.75243e-11);
  fit_X.SetParameter(4,-1.65776e-14);
  fit_X.SetParameter(5, 1.57655e-18);

  fit_Y.SetParameter(0,    0.681069);
  fit_Y.SetParameter(1,  7.6696e-06);
  fit_Y.SetParameter(2,-3.23118e-08);
  fit_Y.SetParameter(3, 1.80031e-12);
  fit_Y.SetParameter(4, 3.41852e-14);

  fit_CX1.SetParameter(0, 0.693777    );
  fit_CX1.SetParameter(1, 36.3197     );
  fit_CX1.SetParameter(2, 11046.8     );
  fit_CX1.SetParameter(3, 1.53333e+06 );
  fit_CX1.SetParameter(4, -6.12379e+07);
  fit_CX1.SetParameter(5, -3.16164e+09);
  fit_CX1.SetParameter(6, 6.69334e+10 );

  fit_CX2.SetParameter(0, 1.45563e+00);
  fit_CX2.SetParameter(1,0);
  //fit_CX2.SetParameter(1,-1.02674e-02);
  fit_CX2.SetParameter(2, 4.71313e-02);
  fit_CX2.SetParameter(3, 2.53096e-01);
  fit_CX2.SetParameter(4,0);
  //fit_CX2.SetParameter(4,-2.53408e-02);
  fit_CX2.SetParameter(5, 1.89370e-01);

  fit_CY1.SetParameter(0,4.24616);

  fit_CY2.SetParameter(0, 4.40000e-01);
  fit_CY2.SetParameter(1, 0.00000e+00);
  fit_CY2.SetParameter(2, 1.56600e-01);
  fit_CY2.SetParameter(3, 8.95000e+00);
  fit_CY2.SetParameter(4,-5.80000e-04);
  fit_CY2.SetParameter(5, 1.53000e-02);

  fit_CZ1.SetParameter(0,-6039.17);
  fit_CZ1.SetParameter(1,-6041.23);

  fit_CZ2.SetParameter(0,-9.83103e+03);
  fit_CZ2.SetParameter(1,-9.83092e+03);

  fit_CZ3.SetParameter(0,-1.34096e+01);
  fit_CZ3.SetParameter(1,-1.36575e+01);
  
  fit_E.SetParameter(0,0.0659346);
  fit_E.SetParameter(1,0.0835921);

  fit_T.SetParameter(0,3.40000e+01);
  fit_T.SetParameter(1,2.18000e-06);
  fit_T.SetParameter(2,2.30000e-07);
  fit_T.SetParameter(3,3.90000e+00);
  fit_T.SetParameter(4,2.18000e-06);
  fit_T.SetParameter(5,7.70000e-07);


  int events = 100000;//Update me!
  double x_X(0);
  double y_X(0);
  double x_Y(0);
  double y_Y(0);
  double x_CX(0);
  double y_CX(0);
  double x_CY(0);
  double y_CY(0);
  double x_CZ(0);
  double y_CZ(0);
  double x_E(0);
  double y_E(0);
  double x_T(0);
  double y_T(0);
  bool set_X(false);
  bool set_Y(false);
  bool set_CX(false);
  bool set_CY(false);
  bool set_CZ(false);
  bool set_E(false);
  bool set_T(false);
  bool fit_failed(false);
  pfile = fopen("TheOutput.txt","w");
  for(int i = 0; i < events; ++i)
  {
    //cout << i << endl;
    fit_failed = false;
    if(!set_X)
    {
      x_X = (X_high-X_low)*double(generator())/double(generator.max()) + X_low;
      y_X = X_ymax*double(generator())/double(generator.max());
    }
    if(!set_Y)
    {
      x_Y = (Y_high-Y_low)*double(generator())/double(generator.max()) + Y_low;
      y_Y = Y_ymax*double(generator())/double(generator.max());
    }
    if(!set_CX)
    {
      x_CX = (CX_high2-CX_low2)*double(generator())/double(generator.max()) + CX_low2;
      if( x_CX > CX_low1 && x_CX < CX_high1 )
      {
        y_CX = CX_ymax1*double(generator())/double(generator.max());
      }
      else
      {
        y_CX = CX_ymax2*double(generator())/double(generator.max());
      }
    }
    if(!set_CY)
    {
      x_CY = (CY_high2-CY_low2)*double(generator())/double(generator.max()) + CY_low2;
      if( x_CY > CY_low1 && x_CY < CY_high1 )
      {
        y_CY = CY_ymax1*double(generator())/double(generator.max());
      }
      else
      {
        y_CY = CY_ymax2*double(generator())/double(generator.max());
      }
    }
    if(!set_CZ)
    {
      x_CZ = (CZ_high2-CZ_low2)*double(generator())/double(generator.max()) + CZ_low2;
      if( x_CZ > CZ_low1 && x_CZ < CZ_high1 )
      {
        y_CZ = CZ_ymax1*double(generator())/double(generator.max());
      }
      else if( x_CZ > CZ_low2 && x_CZ < CZ_high2 )
      {
        y_CZ = CZ_ymax2*double(generator())/double(generator.max());
      }
      else
      {
        y_CZ = CZ_ymax3*double(generator())/double(generator.max());
      }
    }
    if(!set_E)
    {
      x_E = (E_high-E_low)*double(generator())/double(generator.max()) + E_low;
      y_E = E_ymax*double(generator())/double(generator.max());
    }
    if(!set_T)
    {
      x_T = (T_high-T_low)*double(generator())/double(generator.max()) + T_low;
      y_T = T_ymax*double(generator())/double(generator.max());
    }

    set_X = true;
    set_Y = true;
    set_CX = true;
    set_CY = true;
    set_CZ = true;
    set_E = true;
    set_T = true;

    if( y_X > fit_X.Eval(x_X) )
    {
      set_X = false;
      fit_failed = true;
      //cerr << "X failed" << endl;
    }
    if( y_Y > fit_Y.Eval(x_Y) )
    {
      set_Y = false;
      fit_failed = true;
      //cerr << "Y failed" << endl;
    }
    if( y_E > fit_E.Eval(x_E) )
    {
      set_E = false;
      fit_failed = true;
      //cerr << "E failed" << endl;
    }
    if( y_T > fit_T.Eval(x_T) )
    {
      set_T = false;
      fit_failed = true;
      //cerr << "T failed" << endl;
    }
    if( x_CX > CX_low1 && x_CX < CX_high1 )
    {
      if( y_CX > fit_CX1.Eval(x_CX) )
      {
        set_CX = false;
        fit_failed = true;
      //cerr << "CX failed" << endl;
      }
    }
    else
    {
      if( y_CX > fit_CX2.Eval(x_CX) )
      {
        set_CX = false;
        fit_failed = true;
      //cerr << "CX failed" << endl;
      }
    }
    if( x_CY > CY_low1 && x_CY < CY_high1 )
    {
      if( y_CY > fit_CY1.Eval(x_CY) )
      {
        set_CY = false;
        fit_failed = true;
      //cerr << "CY failed" << endl;
      }
    }
    else
    {
      if( y_CY > fit_CY2.Eval(x_CY) )
      {
        set_CY = false;
        fit_failed = true;
      //cerr << "CY failed" << endl;
      }
    }
    if( x_CZ > CZ_low1 && x_CZ < CZ_high1 )
    {
      if( y_CZ > fit_CZ1.Eval(x_CZ) )
      {
        set_CZ = false;
        fit_failed = true;
      //cerr << "CZ failed in 1, " << y_CZ << " > " << fit_CZ1.Eval(x_CZ) << endl;
      }
    }
    else if(x_CZ > CZ_low2 && x_CZ < CZ_high2 )
    {
      if( y_CZ > double(fit_CZ2.Eval(x_CZ)) )
      {
        set_CZ = false;
        fit_failed = true;
      //cerr << "CZ failed in 2, " << y_CZ << " > " << double(fit_CZ2.Eval(x_CZ))  << endl;
      }
    }
    else
    {
      if( y_CZ > fit_CZ3.Eval(x_CZ) )
      {
        set_CZ = false;
        fit_failed = true;
      //cerr << "CZ failed in 3, " << y_CZ << " > " << fit_CZ3.Eval(x_CZ)  << endl;
      }
    }
    if(fit_failed)
    {
      i--;
      continue;
    }
    histo_X.Fill(x_X);
    histo_Y.Fill(x_Y);
    histo_CX.Fill(x_CX);
    histo_CY.Fill(x_CY);
    histo_CZ.Fill(x_CZ);
    histo_E.Fill(x_E);
    histo_T.Fill(x_T);
    fprintf(pfile,"%E %E %E %E %E %E %E\n",x_X,x_Y,x_CX,x_CY,x_CZ,x_E,x_T);
    set_X = false;
    set_Y = false;
    set_CX = false;
    set_CY = false;
    set_CZ = false;
    set_E = false;
    set_T = false;
  }
  fclose(pfile);
  TCanvas c;
  histo_X.Draw();
  c.Print("X.pdf");
  histo_Y.Draw();
  c.Print("Y.pdf");
  histo_CX.Draw();
  c.Print("CX.pdf");
  histo_CY.Draw();
  c.Print("CY.pdf");
  histo_CZ.Draw();
  c.Print("CZ.pdf");
  histo_E.Draw();
  c.Print("E.pdf");
  histo_T.Draw();
  c.Print("T.pdf");
  return 0;
}
