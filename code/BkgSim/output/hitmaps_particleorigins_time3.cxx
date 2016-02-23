void hitmaps_particleorigins_time3()
{
//=========Macro generated from canvas: canvas3/canvas
//=========  (Tue Feb 23 21:54:40 2016) by ROOT version6.04/14
   TCanvas *canvas3 = new TCanvas("canvas3", "canvas",0,0,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canvas3->SetHighLightColor(2);
   canvas3->Range(-4917.721,-429.6345,3943.038,2255.581);
   canvas3->SetFillColor(0);
   canvas3->SetBorderMode(0);
   canvas3->SetBorderSize(2);
   canvas3->SetTickx(1);
   canvas3->SetTicky(1);
   canvas3->SetLeftMargin(0.16);
   canvas3->SetRightMargin(0.05);
   canvas3->SetTopMargin(0.05);
   canvas3->SetBottomMargin(0.16);
   canvas3->SetFrameBorderMode(0);
   canvas3->SetFrameBorderMode(0);
   
   TH2D *histo3 = new TH2D("histo3","20ns <= Time < 50ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap",702,-3500,3500,106,0,2121.32);
   histo3->SetBinContent(1036,3);
   histo3->SetBinContent(1041,1);
   histo3->SetBinContent(1056,39);
   histo3->SetBinContent(1057,1);
   histo3->SetBinContent(1065,1);
   histo3->SetBinContent(1066,1);
   histo3->SetBinContent(1722,3);
   histo3->SetBinContent(1813,1);
   histo3->SetBinContent(2549,1);
   histo3->SetBinContent(7400,1);
   histo3->SetBinContent(8203,1);
   histo3->SetBinContent(8692,2);
   histo3->SetBinContent(8780,1);
   histo3->SetBinContent(11754,2);
   histo3->SetBinContent(36413,3);
   histo3->SetEntries(61);
   histo3->SetStats(0);
   histo3->SetContour(20);
   histo3->SetContourLevel(0,0);
   histo3->SetContourLevel(1,1.95);
   histo3->SetContourLevel(2,3.9);
   histo3->SetContourLevel(3,5.85);
   histo3->SetContourLevel(4,7.8);
   histo3->SetContourLevel(5,9.75);
   histo3->SetContourLevel(6,11.7);
   histo3->SetContourLevel(7,13.65);
   histo3->SetContourLevel(8,15.6);
   histo3->SetContourLevel(9,17.55);
   histo3->SetContourLevel(10,19.5);
   histo3->SetContourLevel(11,21.45);
   histo3->SetContourLevel(12,23.4);
   histo3->SetContourLevel(13,25.35);
   histo3->SetContourLevel(14,27.3);
   histo3->SetContourLevel(15,29.25);
   histo3->SetContourLevel(16,31.2);
   histo3->SetContourLevel(17,33.15);
   histo3->SetContourLevel(18,35.1);
   histo3->SetContourLevel(19,37.05);
   
   TPaletteAxis *palette = new TPaletteAxis(3544.304,0,3938.608,2121.32,histo3);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.05);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.05);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#f1f2e1");
   palette->SetFillColor(ci);
   palette->SetFillStyle(1001);
   histo3->GetListOfFunctions()->Add(palette,"br");
   histo3->SetLineWidth(2);
   histo3->SetMarkerStyle(20);
   histo3->SetMarkerSize(1.2);
   histo3->GetXaxis()->SetTitle("z [mm]");
   histo3->GetXaxis()->SetLabelFont(42);
   histo3->GetXaxis()->SetLabelSize(0.05);
   histo3->GetXaxis()->SetTitleSize(0.05);
   histo3->GetXaxis()->SetTitleOffset(1.4);
   histo3->GetXaxis()->SetTitleFont(42);
   histo3->GetYaxis()->SetTitle("r [mm]");
   histo3->GetYaxis()->SetLabelFont(42);
   histo3->GetYaxis()->SetLabelSize(0.05);
   histo3->GetYaxis()->SetTitleSize(0.05);
   histo3->GetYaxis()->SetTitleOffset(1.4);
   histo3->GetYaxis()->SetTitleFont(42);
   histo3->GetZaxis()->SetTitle("# of origins");
   histo3->GetZaxis()->SetLabelFont(42);
   histo3->GetZaxis()->SetLabelSize(0.05);
   histo3->GetZaxis()->SetTitleSize(0.05);
   histo3->GetZaxis()->SetTitleFont(42);
   histo3->Draw("colz");
   canvas3->Modified();
   canvas3->cd();
   canvas3->SetSelected(canvas3);
}
