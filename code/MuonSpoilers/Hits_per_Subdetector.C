{				
				int n = 12;
				const char *detectors[n] = {"EcalBarrel","EcalEndcap","HcalBarrel","HcalEndcap","MuonBarrel","MuonEndcap","LumiCal","BeamCal","SiVertexBarrel","SiVertexEndcap","SiTrackerBarrel","SiTrackerEndcap"};

				TH1F *h = new TH1F("Hits","Number of hits in the different SiD subdetectors",n,0,n);
				h->GetYaxis()->SetTitle("Number of hits");
				h->SetMarkerColor(kBlue);
				h->SetMarkerStyle(20);
				h->SetFillColor(38);

				h->SetBinContent(1,1);
				h->SetBinContent(2,2);
				h->SetBinContent(3,3);
				h->SetBinContent(4,4);
				h->SetBinContent(5,5);
				h->SetBinContent(6,6);
				h->SetBinContent(7,7);
				h->SetBinContent(8,5);
				h->SetBinContent(9,4);
				h->SetBinContent(10,3);
				h->SetBinContent(11,2);
				h->SetBinContent(12,1);

				h->SetStats(0);

				for (int i =1; i <=n; ++i) {
					h->GetXaxis()->SetBinLabel(i,detectors[i-1]);
				}
				
				TCanvas *canvas = new TCanvas();
				canvas->SetGrid();

				h->Draw();
				canvas.Print("Hits_in_SiD_subdetectors.pdf");
				canvas.Print("Hits_in_SiD_subdetectors.cxx");
}

