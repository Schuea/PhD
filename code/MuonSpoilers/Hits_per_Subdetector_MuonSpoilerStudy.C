{				
				int n = 12;
				const char *detectors[12] = {"EcalBarrel","EcalEndcap","HcalBarrel","HcalEndcap","MuonBarrel","MuonEndcap","LumiCal","BeamCal","SiVertexBarrel","SiVertexEndcap","SiTrackerBarrel","SiTrackerEndcap"};

				TH1F *h1 = new TH1F("Hits","Number of hits in SiD per train - 5 Spoilers vs. 5 Spoilers+Wall ",n,0,n);
				h1->GetYaxis()->SetTitle("Number of hits");
				h1->SetMarkerColor(kBlue);
				h1->SetMarkerStyle(20);
				h1->SetFillColor(38);

				h1->SetBinContent(1,17033);
				h1->SetBinContent(2,31946);
				h1->SetBinContent(3,53566);
				h1->SetBinContent(4,53517);
				h1->SetBinContent(5,12548);
				h1->SetBinContent(6,145256);
				h1->SetBinContent(7,1042);
				h1->SetBinContent(8,968);
				h1->SetBinContent(9,5);
				h1->SetBinContent(10,27);
				h1->SetBinContent(11,4003);
				h1->SetBinContent(12,9188);

				h1->SetStats(0);

        TH1F *h2 = new TH1F("Hits","Number of hits in SiD per train - 5 Spoilers vs. 5 Spoilers+Wall",n,0,n);
				h2->GetYaxis()->SetTitle("Number of hits");
				h2->SetMarkerColor(kRed);
				h2->SetMarkerStyle(20);
				h2->SetFillColor(46);

				h2->SetBinContent(1,2239);
				h2->SetBinContent(2,2340);
				h2->SetBinContent(3,7007);
				h2->SetBinContent(4,3845);
				h2->SetBinContent(5,4415);
				h2->SetBinContent(6,33308);
				h2->SetBinContent(7,71);
				h2->SetBinContent(8,71);
				h2->SetBinContent(9,0);
				h2->SetBinContent(10,1);
				h2->SetBinContent(11,482);
				h2->SetBinContent(12,863);

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

