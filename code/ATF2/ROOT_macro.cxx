{
  std::string Sampler_Tree_name = "Sampler_window";
  std::string inputfile_name = "/media/anne/USB\ DISK/Halo_1000000_Collimatoraperture12mm.root";
  std::cout << Sampler_Tree_name << std::endl;
  std::cout << inputfile_name << std::endl;


  //Make canvas for drawing histograms on
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  //Make histogram for storing the information
  std::string const title1 = "xy distribution of primaries;x [m];y [m]";
  TH2D *histo_xy_Primaries = new TH2D("xy_Primaries", title1.c_str(), 1000, -0.5, 0.5, 50, -0.05, 0.05);
  std::string const title2 = "xy distribution of secondaries;x [m];y [m]";
  TH2D *histo_xy_Secondaries = new TH2D("xy_Secondaries", title2.c_str(), 1000, -0.5, 0.5, 50, -0.05, 0.05);

  //Open root file and get the TTree
  TFile *file = TFile::Open(inputfile_name.c_str());
  TTree *tree = (TTree*)file->Get(Sampler_Tree_name.c_str());

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("parentID", 1);
  tree->SetBranchStatus("x", 1);
  tree->SetBranchStatus("y", 1);

  int parentID(0);
  float x(0.0), y(0.0);

  tree->SetBranchAddress("parentID", &parentID);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);

  long long int const entries = tree->GetEntries();
  for (long long int i = 0; i < entries; ++i) {
    tree->GetEntry(i);
    if (parentID == 0) {
      histo_xy_Primaries->Fill(x,y);
    } else {
      histo_xy_Secondaries->Fill(x,y);
    }
  } 

  //Draw first histogram
  histo_xy_Primaries->Draw("colz");
  canvas->Print("Histo_xy_Primaries_Halo.pdf");
  canvas->Print("Histo_xy_Primaries_Halo.cxx");


  //Draw second histogram
  histo_xy_Secondaries->Draw("colz");
  canvas->Print("Histo_xy_Secondaries_Halo.pdf");
  canvas->Print("Histo_xy_Secondaries_Halo.cxx");

  file->Close();
}
