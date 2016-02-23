{
//=========Macro generated from canvas: canvas/canvas
//=========  (Tue Feb 23 18:01:08 2016) by ROOT version5.34/05
   TCanvas *canvas = new TCanvas("canvas", "canvas",0,0,800,600);
   gStyle->SetOptStat(0);
   canvas->SetHighLightColor(2);
   canvas->Range(-1875,-1875,1875,1875);
   canvas->SetFillColor(0);
   canvas->SetBorderMode(0);
   canvas->SetBorderSize(2);
   canvas->SetFrameBorderMode(0);
   canvas->SetFrameBorderMode(0);
   
   TH2D *histo = new TH2D("histo","Hit maps of cells with full buffer for subdetector EcalEndcap",70,-1500,1500,70,-1500,1500);
   histo->SetBinContent(2123,0.125);
   histo->SetBinContent(2196,0.25);
   histo->SetBinContent(2345,0.125);
   histo->SetBinContent(2407,0.125);
   histo->SetBinContent(2417,0.125);
   histo->SetBinContent(2418,0.125);
   histo->SetBinContent(2840,0.125);
   histo->SetEntries(8);
   histo->SetContour(20);
   histo->SetContourLevel(0,0);
   histo->SetContourLevel(1,0.0125);
   histo->SetContourLevel(2,0.025);
   histo->SetContourLevel(3,0.0375);
   histo->SetContourLevel(4,0.05);
   histo->SetContourLevel(5,0.0625);
   histo->SetContourLevel(6,0.075);
   histo->SetContourLevel(7,0.0875);
   histo->SetContourLevel(8,0.1);
   histo->SetContourLevel(9,0.1125);
   histo->SetContourLevel(10,0.125);
   histo->SetContourLevel(11,0.1375);
   histo->SetContourLevel(12,0.15);
   histo->SetContourLevel(13,0.1625);
   histo->SetContourLevel(14,0.175);
   histo->SetContourLevel(15,0.1875);
   histo->SetContourLevel(16,0.2);
   histo->SetContourLevel(17,0.2125);
   histo->SetContourLevel(18,0.225);
   histo->SetContourLevel(19,0.2375);
   
   TPaletteAxis *palette = new TPaletteAxis(1518.75,-1500,1687.5,1500,histo);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.035);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.035);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   histo->GetListOfFunctions()->Add(palette,"br");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   histo->SetLineColor(ci);
   histo->GetXaxis()->SetTitle("x [mm]");
   histo->GetXaxis()->SetLabelFont(42);
   histo->GetXaxis()->SetLabelSize(0.035);
   histo->GetXaxis()->SetTitleSize(0.035);
   histo->GetXaxis()->SetTitleFont(42);
   histo->GetYaxis()->SetTitle(" y [mm]");
   histo->GetYaxis()->SetLabelFont(42);
   histo->GetYaxis()->SetLabelSize(0.035);
   histo->GetYaxis()->SetTitleSize(0.035);
   histo->GetYaxis()->SetTitleFont(42);
   histo->GetZaxis()->SetTitle("# of hits");
   histo->GetZaxis()->SetLabelFont(42);
   histo->GetZaxis()->SetLabelSize(0.035);
   histo->GetZaxis()->SetTitleSize(0.035);
   histo->GetZaxis()->SetTitleFont(42);
   histo->Draw("colz");
   
   TPaveText *pt = new TPaveText(0.01,0.9404546,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   TText *text = pt->AddText("Hit maps of cells with full buffer for subdetector EcalEndcap");
   pt->Draw();
   canvas->Modified();
   canvas->cd();
   canvas->SetSelected(canvas);
}
