{
//=========Macro generated from canvas: canvas/canvas
//=========  (Tue Feb 23 14:43:55 2016) by ROOT version5.34/05
   TCanvas *canvas = new TCanvas("canvas", "canvas",0,0,800,600);
   canvas->SetHighLightColor(2);
   canvas->Range(-4375,-50,4375,450);
   canvas->SetFillColor(0);
   canvas->SetBorderMode(0);
   canvas->SetBorderSize(2);
   canvas->SetFrameBorderMode(0);
   canvas->SetFrameBorderMode(0);
   
   TH2D *histo2 = new TH2D("histo2","10ns <= Time < 20ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap",702,-3500,3500,106,0,400);
   histo2->SetBinContent(1056,1);
   histo2->SetEntries(1);
   histo2->SetContour(20);
   histo2->SetContourLevel(0,0);
   histo2->SetContourLevel(1,0.05);
   histo2->SetContourLevel(2,0.1);
   histo2->SetContourLevel(3,0.15);
   histo2->SetContourLevel(4,0.2);
   histo2->SetContourLevel(5,0.25);
   histo2->SetContourLevel(6,0.3);
   histo2->SetContourLevel(7,0.35);
   histo2->SetContourLevel(8,0.4);
   histo2->SetContourLevel(9,0.45);
   histo2->SetContourLevel(10,0.5);
   histo2->SetContourLevel(11,0.55);
   histo2->SetContourLevel(12,0.6);
   histo2->SetContourLevel(13,0.65);
   histo2->SetContourLevel(14,0.7);
   histo2->SetContourLevel(15,0.75);
   histo2->SetContourLevel(16,0.8);
   histo2->SetContourLevel(17,0.85);
   histo2->SetContourLevel(18,0.9);
   histo2->SetContourLevel(19,0.95);
   
   TPaletteAxis *palette = new TPaletteAxis(3543.75,0,3937.5,400,histo2);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.035);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.035);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   histo2->GetListOfFunctions()->Add(palette,"br");
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.575,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("histo2");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 1      ");
   text = ptstats->AddText("Mean x =      0");
   text = ptstats->AddText("Mean y =      0");
   text = ptstats->AddText("RMS x =      0");
   text = ptstats->AddText("RMS y =      0");
   text = ptstats->AddText("       0|      0|      0
");
   text = ptstats->AddText("       0|      1|      0
");
   text = ptstats->AddText("       0|      0|      0
");
   ptstats->SetOptStat(11111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   histo2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(histo2);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   histo2->SetLineColor(ci);
   histo2->GetXaxis()->SetTitle("z [mm]");
   histo2->GetXaxis()->SetLabelFont(42);
   histo2->GetXaxis()->SetLabelSize(0.035);
   histo2->GetXaxis()->SetTitleSize(0.035);
   histo2->GetXaxis()->SetTitleFont(42);
   histo2->GetYaxis()->SetTitle(" r [mm]");
   histo2->GetYaxis()->SetLabelFont(42);
   histo2->GetYaxis()->SetLabelSize(0.035);
   histo2->GetYaxis()->SetTitleSize(0.035);
   histo2->GetYaxis()->SetTitleFont(42);
   histo2->GetZaxis()->SetLabelFont(42);
   histo2->GetZaxis()->SetLabelSize(0.035);
   histo2->GetZaxis()->SetTitleSize(0.035);
   histo2->GetZaxis()->SetTitleFont(42);
   histo2->Draw("colz");
   
   TPaveText *pt = new TPaveText(0.15,0.9341608,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("10ns <= Time < 20ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap");
   pt->Draw();
   canvas->Modified();
   canvas->cd();
   canvas->SetSelected(canvas);
}
