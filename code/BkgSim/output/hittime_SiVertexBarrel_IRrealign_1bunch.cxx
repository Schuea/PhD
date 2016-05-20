{
//=========Macro generated from canvas: canvas1/canvas
//=========  (Fri May 20 10:01:52 2016) by ROOT version5.34/05
   TCanvas *canvas1 = new TCanvas("canvas1", "canvas",0,0,800,600);
   gStyle->SetOptFit(1);
   canvas1->SetHighLightColor(2);
   canvas1->Range(-16.9697,-14.14214,89.09091,114.4227);
   canvas1->SetFillColor(0);
   canvas1->SetBorderMode(0);
   canvas1->SetBorderSize(2);
   canvas1->SetTickx(1);
   canvas1->SetTicky(1);
   canvas1->SetLeftMargin(0.16);
   canvas1->SetRightMargin(0.18);
   canvas1->SetTopMargin(0.12);
   canvas1->SetFrameBorderMode(0);
   canvas1->SetFrameBorderMode(0);
   
   TH2D *HitTime_SiVertexBarrel_bunch#1 = new TH2D("HitTime_SiVertexBarrel_bunch#1","Radial position of hits over hit time for SiVertexBarrel",7,0,70,9,0,98.99495);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(19,889);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(21,262);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(22,6);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(28,274);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(30,25);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(37,66);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(39,78);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(46,21);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(55,6);
   HitTime_SiVertexBarrel_bunch#1->SetBinContent(57,1);
   HitTime_SiVertexBarrel_bunch#1->SetEntries(1628);
   
   TPaveStats *ptstats = new TPaveStats(0.65,0.83,0.85,0.9,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("HitTime_SiVertexBarrel_bunch#1");
   text->SetTextSize(0.0322);
   text = ptstats->AddText("Entries = 1628   ");
   ptstats->SetOptStat(11);
   ptstats->SetOptFit(1111);
   ptstats->Draw();
   HitTime_SiVertexBarrel_bunch#1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(HitTime_SiVertexBarrel_bunch#1);
   HitTime_SiVertexBarrel_bunch#1->SetLineWidth(2);
   HitTime_SiVertexBarrel_bunch#1->SetMarkerStyle(20);
   HitTime_SiVertexBarrel_bunch#1->SetMarkerSize(1.2);
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetTitle("Hit time [ns]");
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetLabelFont(42);
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetLabelSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetTitleSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetTitleOffset(1.05);
   HitTime_SiVertexBarrel_bunch#1->GetXaxis()->SetTitleFont(42);
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetTitle("r [mm]");
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetLabelFont(42);
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetLabelSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetTitleSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetTitleOffset(1.4);
   HitTime_SiVertexBarrel_bunch#1->GetYaxis()->SetTitleFont(42);
   HitTime_SiVertexBarrel_bunch#1->GetZaxis()->SetLabelFont(42);
   HitTime_SiVertexBarrel_bunch#1->GetZaxis()->SetLabelSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetZaxis()->SetTitleSize(0.05);
   HitTime_SiVertexBarrel_bunch#1->GetZaxis()->SetTitleOffset(1.3);
   HitTime_SiVertexBarrel_bunch#1->GetZaxis()->SetTitleFont(42);
   HitTime_SiVertexBarrel_bunch#1->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9383566,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   text = pt->AddText("Radial position of hits over hit time for SiVertexBarrel");
   pt->Draw();
   canvas1->Modified();
   canvas1->cd();
   canvas1->SetSelected(canvas1);
}
