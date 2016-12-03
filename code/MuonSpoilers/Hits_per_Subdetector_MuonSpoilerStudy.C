{				
				int n = 12;
				const char *detectors[12] = {"EcalBarrel","EcalEndcap","HcalBarrel","HcalEndcap","MuonBarrel","MuonEndcap","LumiCal","BeamCal","SiVertexBarrel","SiVertexEndcap","SiTrackerBarrel","SiTrackerEndcap"};

				TH1F *h1 = new TH1F("Hits","Number of hits in SiD per train - 5 Spoilers vs. 5 Spoilers+Wall ",n,0,n);
				h1->GetYaxis()->SetTitle("Number of hits");
				h1->SetMarkerColor(kBlue);
				h1->SetMarkerStyle(20);
				h1->SetFillColor(38);

				h1->SetBinContent(1,16997);
				h1->SetBinContent(2,31927);
				h1->SetBinContent(3,53408);
				h1->SetBinContent(4,53394);
				h1->SetBinContent(5,12505);
				h1->SetBinContent(6,145212);
				h1->SetBinContent(7,1038);
				h1->SetBinContent(8,965);
				h1->SetBinContent(9,2);
				h1->SetBinContent(10,0);
				h1->SetBinContent(11,4012);
				h1->SetBinContent(12,9196);

				h1->SetStats(0);

        TH1F *h2 = new TH1F("Hits","Number of hits in SiD per train - 5 Spoilers vs. 5 Spoilers+Wall",n,0,n);
				h2->GetYaxis()->SetTitle("Number of hits");
				h2->SetMarkerColor(kRed);
				h2->SetMarkerStyle(20);
				h2->SetFillColor(46);

				h2->SetBinContent(1,2218);
				h2->SetBinContent(2,2337);
				h2->SetBinContent(3,6982);
				h2->SetBinContent(4,3829);
				h2->SetBinContent(5,4344);
				h2->SetBinContent(6,33179);
				h2->SetBinContent(7,69);
				h2->SetBinContent(8,69);
				h2->SetBinContent(9,0);
				h2->SetBinContent(10,0);
				h2->SetBinContent(11,486);
				h2->SetBinContent(12,867);

				h2->SetStats(0);

				for (int i =1; i <=n; ++i) {
					h1->GetXaxis()->SetBinLabel(i,detectors[i-1]);
					h2->GetXaxis()->SetBinLabel(i,detectors[i-1]);
				}
				h1->GetXaxis()->SetLabelSize(0.045);
				h1->GetYaxis()->SetLabelSize(0.05);
				h1->GetYaxis()->SetTitleSize(0.05);
				h1->GetYaxis()->SetTitleOffset(0.95);
				h2->GetXaxis()->SetLabelSize(0.045);
				h2->GetYaxis()->SetLabelSize(0.05);
				h2->GetYaxis()->SetTitleSize(0.05);
				h2->GetYaxis()->SetTitleOffset(0.95);
				
				TCanvas *canvas = new TCanvas();
				canvas->SetGrid();
				canvas->SetLogy();

				h1->Draw();
				h2->Draw("SAME");
				canvas->Print("Hits_in_SiD_subdetectors_MuonSpoilerStudy.pdf");
				canvas->Print("Hits_in_SiD_subdetectors_MuonSpoilerStudy.cxx");
}

