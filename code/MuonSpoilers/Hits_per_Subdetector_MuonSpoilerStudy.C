{				
				int n = 12;
				const char *detectors[12] = {"EcalBarrel","EcalEndcap","HcalEndcap","HcalBarrel","MuonBarrel","MuonEndcap","LumiCal","BeamCal","SiVertexEndcap","SiVertexBarrel","SiTrackerEndcap","SiTrackerBarrel"};

				TH1F *h = new TH1F("Hits","Number of hits in the SiD per train - Spoiler + Wall",n,0,n);
				h->GetYaxis()->SetTitle("Number of hits");
				h->SetMarkerColor(kBlue);
				h->SetMarkerStyle(20);
				h->SetFillColor(38);

				h->SetBinContent(1,9423/5);
				h->SetBinContent(2,12745/5);
				h->SetBinContent(3,21024/5);
				h->SetBinContent(4,39856/5);
				h->SetBinContent(5,27569/5);
				h->SetBinContent(6,198342/5);
				h->SetBinContent(7,212/5);
				h->SetBinContent(8,68/5);
				h->SetBinContent(9,0);
				h->SetBinContent(10,0);
				h->SetBinContent(11,3849/5);
				h->SetBinContent(12,1892/5);

				h->SetStats(0);

				for (int i =1; i <=n; ++i) {
					h->GetXaxis()->SetBinLabel(i,detectors[i-1]);
				}
				
				TCanvas *canvas = new TCanvas();
				canvas->SetGrid();
				canvas->SetLogy();

				h->Draw();
				canvas->Print("Hits_in_SiD_subdetectors_Spoiler.pdf");
				canvas->Print("Hits_in_SiD_subdetectors_Spoiler.cxx");
}

