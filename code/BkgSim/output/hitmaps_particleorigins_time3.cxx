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
   
   TH2D *histo3 = new TH2D("histo3","20ns <= Time < 50ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap",702,-3500,3500,106,0,400);
   histo3->SetBinContent(1056,23);
   histo3->SetBinContent(3175,1);
   histo3->SetBinContent(3888,1);
   histo3->SetBinContent(4239,1);
   histo3->SetBinContent(64568,3);
   histo3->SetEntries(29);
   histo3->SetContour(20);
   histo3->SetContourLevel(0,0);
   histo3->SetContourLevel(1,1.15);
   histo3->SetContourLevel(2,2.3);
   histo3->SetContourLevel(3,3.45);
   histo3->SetContourLevel(4,4.6);
   histo3->SetContourLevel(5,5.75);
   histo3->SetContourLevel(6,6.9);
   histo3->SetContourLevel(7,8.05);
   histo3->SetContourLevel(8,9.2);
   histo3->SetContourLevel(9,10.35);
   histo3->SetContourLevel(10,11.5);
   histo3->SetContourLevel(11,12.65);
   histo3->SetContourLevel(12,13.8);
   histo3->SetContourLevel(13,14.95);
   histo3->SetContourLevel(14,16.1);
   histo3->SetContourLevel(15,17.25);
   histo3->SetContourLevel(16,18.4);
   histo3->SetContourLevel(17,19.55);
   histo3->SetContourLevel(18,20.7);
   histo3->SetContourLevel(19,21.85);
   
   TPaletteAxis *palette = new TPaletteAxis(3543.75,0,3937.5,400,histo3);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.035);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.035);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   histo3->GetListOfFunctions()->Add(palette,"br");
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.575,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("histo3");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 29     ");
   text = ptstats->AddText("Mean x =  49.98");
   text = ptstats->AddText("Mean y =  37.23");
   text = ptstats->AddText("RMS x =  792.6");
   text = ptstats->AddText("RMS y =    104");
   text = ptstats->AddText("       0|      0|      0
");
   text = ptstats->AddText("       0|     29|      0
");
   text = ptstats->AddText("       0|      0|      0
");
   ptstats->SetOptStat(11111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   histo3->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(histo3);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   histo3->SetLineColor(ci);
   histo3->GetXaxis()->SetTitle("z [mm]");
   histo3->GetXaxis()->SetLabelFont(42);
   histo3->GetXaxis()->SetLabelSize(0.035);
   histo3->GetXaxis()->SetTitleSize(0.035);
   histo3->GetXaxis()->SetTitleFont(42);
   histo3->GetYaxis()->SetTitle(" r [mm]");
   histo3->GetYaxis()->SetLabelFont(42);
   histo3->GetYaxis()->SetLabelSize(0.035);
   histo3->GetYaxis()->SetTitleSize(0.035);
   histo3->GetYaxis()->SetTitleFont(42);
   histo3->GetZaxis()->SetLabelFont(42);
   histo3->GetZaxis()->SetLabelSize(0.035);
   histo3->GetZaxis()->SetTitleSize(0.035);
   histo3->GetZaxis()->SetTitleFont(42);
   histo3->Draw("colz");
   
   TPaveText *pt = new TPaveText(0.15,0.9341608,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("20ns <= Time < 50ns, Hit maps of particle origins for those particles hitting subdetector EcalEndcap");
   pt->Draw();
   canvas->Modified();
   canvas->cd();
   canvas->SetSelected(canvas);
}
