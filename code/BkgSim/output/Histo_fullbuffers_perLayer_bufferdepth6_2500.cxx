{
//=========Macro generated from canvas: canvas/canvas
//=========  (Wed Feb 24 16:57:30 2016) by ROOT version5.34/05
   TCanvas *canvas = new TCanvas("canvas", "canvas",0,0,800,600);
   gStyle->SetOptFit(1);
   canvas->SetHighLightColor(2);
   canvas->Range(-7.272727,-0.9,38.18182,7.281818);
   canvas->SetFillColor(0);
   canvas->SetBorderMode(0);
   canvas->SetBorderSize(2);
   canvas->SetTickx(1);
   canvas->SetTicky(1);
   canvas->SetLeftMargin(0.16);
   canvas->SetRightMargin(0.18);
   canvas->SetTopMargin(0.12);
   canvas->SetFrameBorderMode(0);
   canvas->SetFrameBorderMode(0);
   
   TH1D *histo = new TH1D("histo","Number of cells with full buffer per layer for subdetector EcalEndcap",30,0,30);
   histo->SetBinContent(2,4);
   histo->SetBinContent(3,6);
   histo->SetBinContent(4,6);
   histo->SetBinContent(5,3);
   histo->SetBinContent(7,1);
   histo->SetBinContent(14,1);
   histo->SetBinContent(17,1);
   histo->SetBinContent(25,1);
   histo->SetEntries(23);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.715,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("histo");
   text->SetTextSize(0.06439999);
   text = ptstats->AddText("Entries = 23     ");
   text = ptstats->AddText("Mean  =  4.565");
   text = ptstats->AddText("RMS   =  5.468");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(1111);
   ptstats->Draw();
   histo->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(histo);
   histo->SetLineWidth(2);
   histo->SetMarkerStyle(20);
   histo->SetMarkerSize(1.2);
   histo->GetXaxis()->SetTitle("Layer number");
   histo->GetXaxis()->SetLabelFont(42);
   histo->GetXaxis()->SetLabelSize(0.05);
   histo->GetXaxis()->SetTitleSize(0.05);
   histo->GetXaxis()->SetTitleOffset(1.05);
   histo->GetXaxis()->SetTitleFont(42);
   histo->GetYaxis()->SetTitle("Number of dead cells");
   histo->GetYaxis()->SetLabelFont(42);
   histo->GetYaxis()->SetLabelSize(0.05);
   histo->GetYaxis()->SetTitleSize(0.05);
   histo->GetYaxis()->SetTitleOffset(1.4);
   histo->GetYaxis()->SetTitleFont(42);
   histo->GetZaxis()->SetLabelFont(42);
   histo->GetZaxis()->SetLabelSize(0.05);
   histo->GetZaxis()->SetTitleSize(0.05);
   histo->GetZaxis()->SetTitleOffset(1.3);
   histo->GetZaxis()->SetTitleFont(42);
   histo->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9404546,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   text = pt->AddText("Number of cells with full buffer per layer for subdetector EcalEndcap");
   pt->Draw();
   canvas->Modified();
   canvas->cd();
   canvas->SetSelected(canvas);
}
