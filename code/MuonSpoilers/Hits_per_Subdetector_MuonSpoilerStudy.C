{				
				int n = 12;
				const char *detectors[12] = {"EcalBarrel","EcalEndcap","HcalBarrel","HcalEndcap","MuonBarrel","MuonEndcap","LumiCal","BeamCal","SiVertexBarrel","SiVertexEndcap","SiTrackerBarrel","SiTrackerEndcap"};

        //ILC500: 5spoilers
				TH1F *h1 = new TH1F("Hits","Number of hits in SiD per train - 5 spoilers vs. 5 spoilers+wall ",n,0,n);
				h1->GetYaxis()->SetTitle("Number of hits");
				h1->SetLineColor(kPink-1);
				h1->SetMarkerColor(kPink-1);
				h1->SetMarkerStyle(20);
				h1->SetFillColor(kPink-1);
				//h1->SetFillStyle(3003);

        h1->SetBinContent(1,   15468.9);
        h1->SetBinContent(2,   25258.6);
        h1->SetBinContent(3,   46458.3);
        h1->SetBinContent(4,   43568.6);
        h1->SetBinContent(5,   12551.1);
        h1->SetBinContent(6,   138431.3);
        h1->SetBinContent(7,   990.6);
        h1->SetBinContent(8,   1089.1);
        h1->SetBinContent(9,   0.5);
        h1->SetBinContent(10,  24.9);
        h1->SetBinContent(11,  3104.9);
        h1->SetBinContent(12,  6900.5);

				h1->SetStats(0);

        //ILC500: 5spoilers + wall
        TH1F *h2 = new TH1F("Hits2","Number of hits in SiD per train - 5 spoilers vs. 5 spoilers+wall",n,0,n);
				h2->GetYaxis()->SetTitle("Number of hits");
				h2->SetLineColor(kPink-1);
				h2->SetMarkerColor(kPink-1);
				h2->SetMarkerStyle(20);
				h2->SetFillColor(10);
				//h2->SetFillColorAlpha(10, 0.3);

        h2->SetBinContent(1,  2036.5);
        h2->SetBinContent(2,  1982.3);
        h2->SetBinContent(3,  5540.4);
        h2->SetBinContent(4,  3251.8);
        h2->SetBinContent(5,  3717.5);
        h2->SetBinContent(6,  28431.0);
        h2->SetBinContent(7,  64.8);
        h2->SetBinContent(8,  57.6);
        h2->SetBinContent(9,  0.0);
        h2->SetBinContent(10, 0.9);
        h2->SetBinContent(11, 529.6);
        h2->SetBinContent(12, 768.8);

				h2->SetStats(0);

        //ILC250: 5spoilers
        TH1F *h3 = new TH1F("Hits3","Number of hits in SiD per train - 5 spoilers vs. 5 spoilers+wall ",n,0,n);
				h3->GetYaxis()->SetTitle("Number of hits");
				h3->SetLineColor(kCyan+3);
				h3->SetMarkerColor(kCyan+3);
				h3->SetMarkerStyle(20);
				h3->SetFillColor(kCyan+3);

        h3->SetBinContent(1,  3634.5);
        h3->SetBinContent(2,  6843.8);
        h3->SetBinContent(3,  11294.5);
        h3->SetBinContent(4,  12167.9);
        h3->SetBinContent(5,  3158.7);
        h3->SetBinContent(6,  38184.4);
        h3->SetBinContent(7,  308.2);
        h3->SetBinContent(8,  264.5);
        h3->SetBinContent(9,  0.1);
        h3->SetBinContent(10, 6.4);
        h3->SetBinContent(11, 728.4);
        h3->SetBinContent(12, 1694.9);

				h3->SetStats(0);

        //ILC250: 5spoilers + wall
        TH1F *h4 = new TH1F("Hits4","Number of hits in SiD per train - 5 spoilers vs. 5 spoilers+wall",n,0,n);
				h4->GetYaxis()->SetTitle("Number of hits");
				h4->SetLineColor(kCyan+3);
				h4->SetMarkerColor(kCyan+3);
				h4->SetMarkerStyle(20);
				h4->SetFillColor(10);

        h4->SetBinContent(1,  138.7);
        h4->SetBinContent(2,  132.7);
        h4->SetBinContent(3,  365.5);
        h4->SetBinContent(4,  247.7);
        h4->SetBinContent(5,  178.3);
        h4->SetBinContent(6,  1087.3);
        h4->SetBinContent(7,  2.5);
        h4->SetBinContent(8,  3.2);
        h4->SetBinContent(9,  0.0);
        h4->SetBinContent(10, 0.0);
        h4->SetBinContent(11, 15.0);
        h4->SetBinContent(12, 37.0);
				
				h4->SetStats(0);


				for (int i =1; i <=n; ++i) {
					h1->GetXaxis()->SetBinLabel(i,detectors[i-1]);
					h2->GetXaxis()->SetBinLabel(i,detectors[i-1]);
					h3->GetXaxis()->SetBinLabel(i,detectors[i-1]);
					h4->GetXaxis()->SetBinLabel(i,detectors[i-1]);
				}
				h1->GetXaxis()->SetLabelSize(0.045);
				h1->GetYaxis()->SetLabelSize(0.05);
				h1->GetYaxis()->SetTitleSize(0.05);
				h1->GetYaxis()->SetTitleOffset(0.95);
				
        h2->GetXaxis()->SetLabelSize(0.045);
				h2->GetYaxis()->SetLabelSize(0.05);
				h2->GetYaxis()->SetTitleSize(0.05);
				h2->GetYaxis()->SetTitleOffset(0.95);
				
        h3->GetXaxis()->SetLabelSize(0.045);
				h3->GetYaxis()->SetLabelSize(0.05);
				h3->GetYaxis()->SetTitleSize(0.05);
				h3->GetYaxis()->SetTitleOffset(0.95);
				
        h4->GetXaxis()->SetLabelSize(0.045);
				h4->GetYaxis()->SetLabelSize(0.05);
				h4->GetYaxis()->SetTitleSize(0.05);
				h4->GetYaxis()->SetTitleOffset(0.95);
				
				TCanvas *canvas = new TCanvas();
				canvas->SetGrid();
				canvas->SetLogy();

				h1->Draw();
				h3->Draw("SAME");
				h2->Draw("SAME");
        TH1F *h2_ = (TH1F*)h2->Clone();
				h2_->SetFillStyle(3004);
				h2_->SetFillColor(kPink-1);
				h2_->Draw("SAME");
				h4->Draw("SAME");
        TH1F *h3_ = (TH1F*)h3->Clone();
				h3_->SetFillColorAlpha(kCyan+3,0.0);
				h3_->Draw("SAME");
				h4->Draw("SAME");
        TH1F *h4_ = (TH1F*)h4->Clone();
				h4_->SetFillStyle(3005);
				h4_->SetFillColor(kCyan+3);
				h4_->Draw("SAME");


        TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
        leg->AddEntry(h1,"ILC500, 5 spoilers","f");
        leg->AddEntry(h2_,"ILC500, 5 spoilers + wall","f");
        leg->AddEntry(h3,"ILC250, 5 spoilers","f");
        leg->AddEntry(h4_,"ILC250, 5 spoilers + wall","f");
        leg->Draw();

				canvas->Print("Hits_in_SiD_subdetectors_MuonSpoilerStudy.pdf");
				canvas->Print("Hits_in_SiD_subdetectors_MuonSpoilerStudy.cxx");
}

