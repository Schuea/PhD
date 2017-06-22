#include <math.h>

//{
//#include <math.h>
//gROOT->Reset();
//
//TCanvas canvas;
//canvas.SetGridx();
//
//    // create a 2-d histogram to define the range
//   TH2F *hr = 
//    new TH2F("hr","Special labels on the axis",2,-0.4,1.2,2,0,12);
//   hr->GetYaxis()->SetLabelOffset(1);
//   hr->GetYaxis()->SetNdivisions(10);
//   hr->SetXTitle("X title");
//   hr->Draw();
//   canvas.GetFrame()->SetFillColor(21);
//   canvas.GetFrame()->SetBorderSize(12);
//
//int num_points = 5;
//
//float Time[] = {1.0, 2.0, 3.0, 4.0, 5.0};
////float Design1[] = {8,6,2,1,-2};
//float Design1[] = {1.0*std::pow(10,8),pow(10,6),pow(10,2),pow(10,1),pow(10,-2)};
//float Design2[] = {0.0,0.0,0.0,0.0,0.0};
//
//TText *t = new TText();
//t->SetTextAlign(32);
//t->SetTextSize(0.035);
//t->SetTextFont(72);
//char  *TimeLabels[5] = {"1 minute", "1 hour", "1 day", "1 month", "1 year"};
//
//TGraph g_1(num_points, Time, Design1);
//TGraph g_2(num_points, Time, Design2);
//
//g_1.SetTitle("Dose equivalent rate after a certain cooling time;Cooling time;Dose equivalent rate [mSv/s]");
//g_1.SetMarkerColor(kBlue);
//g_1.SetLineColor(kBlue);
//g_1.SetMarkerStyle(20);
//g_2.SetTitle("Dose equivalent rate after a certain cooling time;Cooling time;Dose equivalent rate [mSv/s]");
//g_2.SetMarkerColor(kMagenta);
//g_2.SetLineColor(kMagenta);
//g_2.SetMarkerStyle(20);
//
//g_1.Draw("APL");
//for(Int_t i=0; i<num_points; ++i){
//	t->DrawText(-0.42,Time[i],TimeLabels[i]);
//}
//canvas.Print("DoseEQ_Time_Design1.pdf");
////g_2.Draw("APL");
////canvas.Print("DoseEQ_Time_Design2.pdf");
////g_1.Draw("APL");
////g_2.Draw("PLSAME");
////TLegend leg(0.6, 0.7, 0.8, 0.89);
////leg.AddEntry(&g_1, "Design 1");
////leg.AddEntry(&g_2, "Design 2");
////leg.Draw();
////canvas.Print("DoseEQ_Time_Comparison.pdf");
//}

void Draw_Points_over_Time()
{
   const Int_t nx = 5;
   char *Time[nx] = {"One minute", "One hour", "One day", "One month", "One year"};
   //float Design1[nx] = {1.0, 2.0, 3.0, 4.0, 5.0};
   double Design1[nx] = {1.0*pow(10.,8.),pow(10.,6.),pow(10.,2.),pow(10.,1.),pow(10.,-2.)};
   TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,900,500);
   c1->SetGrid();
   c1->SetTopMargin(0.15);
   c1->SetLogy();
   gStyle->SetTitleW(0.5); //title width
   gStyle->SetTitleH(0.1); //title height
   TH1F *h1 = new TH1F("h1","test;;y axis [label]",5,0,5);
   h1->SetStats(0);
   h1->SetFillColor(38);
   //h1->SetBit(TH1::kCanRebin);
   for (Int_t i=0; i<nx; i++) {
      h1->Fill(Time[i],Design1[i]);
   }
   h1->LabelsDeflate();
   h1->GetXaxis()->SetLabelSize(0.08);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetTitleSize(0.06);
   h1->GetYaxis()->SetTitleOffset(0.65);
   h1->SetMarkerStyle(20);
   h1->SetMarkerColor(kRed);
   h1->Draw("hist p");
   
   float Design2[nx] = {5.0, 4.0, 3.0, 2.0, 1.0};
   TH1F *h2 = new TH1F("h2","test",5,0,5);
   h2->SetStats(0);
   h2->SetFillColor(38);
   //h2->SetBit(TH1::kCanRebin);
   for (Int_t i=0; i<nx; i++) {
      h2->Fill(Time[i],Design2[i]);
   }
   //h2->LabelsDeflate();
   h2->SetMarkerStyle(20);
   h2->SetMarkerColor(kBlue);
   h2->Draw("hist p same");
   TLegend leg(0.7,0.86,0.95,0.99);
   leg.AddEntry(h1,"Histogram 1","p");
   leg.AddEntry(h2,"Histogram 2","p");
   leg.SetTextSize(0.06);
   leg.Draw();
   c1->Print("Test.pdf");
}

