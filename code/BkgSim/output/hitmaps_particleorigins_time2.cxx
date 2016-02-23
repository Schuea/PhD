void hitmaps_particleorigins_time2()
{
//=========Macro generated from canvas: canvas2/canvas
//=========  (Tue Feb 23 21:54:40 2016) by ROOT version6.04/14
   TCanvas *canvas2 = new TCanvas("canvas2", "canvas",0,0,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canvas2->SetHighLightColor(2);
   canvas2->Range(-4917.721,-429.6345,3943.038,2255.581);
   canvas2->SetFillColor(0);
   canvas2->SetBorderMode(0);
   canvas2->SetBorderSize(2);
   canvas2->SetTickx(1);
   canvas2->SetTicky(1);
   canvas2->SetLeftMargin(0.16);
   canvas2->SetRightMargin(0.05);
   canvas2->SetTopMargin(0.05);
   canvas2->SetBottomMargin(0.16);
   canvas2->SetFrameBorderMode(0);
   canvas2->SetFrameBorderMode(0);
   
   TH2D *histo2 = new TH2D("histo2","10ns <= Time < 20ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap",702,-3500,3500,106,0,2121.32);
   histo2->SetBinContent(1037,1);
   histo2->SetBinContent(1045,1);
   histo2->SetBinContent(1048,2);
   histo2->SetBinContent(1062,1);
   histo2->SetBinContent(1811,3);
   histo2->SetEntries(8);
   histo2->SetStats(0);
   histo2->SetContour(20);
   histo2->SetContourLevel(0,0);
   histo2->SetContourLevel(1,0.15);
   histo2->SetContourLevel(2,0.3);
   histo2->SetContourLevel(3,0.45);
   histo2->SetContourLevel(4,0.6);
   histo2->SetContourLevel(5,0.75);
   histo2->SetContourLevel(6,0.9);
   histo2->SetContourLevel(7,1.05);
   histo2->SetContourLevel(8,1.2);
   histo2->SetContourLevel(9,1.35);
   histo2->SetContourLevel(10,1.5);
   histo2->SetContourLevel(11,1.65);
   histo2->SetContourLevel(12,1.8);
   histo2->SetContourLevel(13,1.95);
   histo2->SetContourLevel(14,2.1);
   histo2->SetContourLevel(15,2.25);
   histo2->SetContourLevel(16,2.4);
   histo2->SetContourLevel(17,2.55);
   histo2->SetContourLevel(18,2.7);
   histo2->SetContourLevel(19,2.85);
   
   TPaletteAxis *palette = new TPaletteAxis(3544.304,0,3938.608,2121.32,histo2);
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
   histo2->GetListOfFunctions()->Add(palette,"br");
   histo2->SetLineWidth(2);
   histo2->SetMarkerStyle(20);
   histo2->SetMarkerSize(1.2);
   histo2->GetXaxis()->SetTitle("z [mm]");
   histo2->GetXaxis()->SetLabelFont(42);
   histo2->GetXaxis()->SetLabelSize(0.05);
   histo2->GetXaxis()->SetTitleSize(0.05);
   histo2->GetXaxis()->SetTitleOffset(1.4);
   histo2->GetXaxis()->SetTitleFont(42);
   histo2->GetYaxis()->SetTitle("r [mm]");
   histo2->GetYaxis()->SetLabelFont(42);
   histo2->GetYaxis()->SetLabelSize(0.05);
   histo2->GetYaxis()->SetTitleSize(0.05);
   histo2->GetYaxis()->SetTitleOffset(1.4);
   histo2->GetYaxis()->SetTitleFont(42);
   histo2->GetZaxis()->SetTitle("# of origins");
   histo2->GetZaxis()->SetLabelFont(42);
   histo2->GetZaxis()->SetLabelSize(0.05);
   histo2->GetZaxis()->SetTitleSize(0.05);
   histo2->GetZaxis()->SetTitleFont(42);
   histo2->Draw("colz");
   canvas2->Modified();
   canvas2->cd();
   canvas2->SetSelected(canvas2);
}
