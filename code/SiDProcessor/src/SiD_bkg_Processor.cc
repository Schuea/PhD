#include "SiD_bkg_Processor.h"
#include "Collection_Processor.h"
#include "OutputCollection_Processor.h"

#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/LCGenericObject.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>

#include "TPaveStats.h"
#include "TExec.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCut.h"
#ifdef USE_MARLIN
// streamlog include
#include "streamlog/streamlog.h"
#endif

using namespace lcio;
using namespace marlin;
using namespace std;


//SID_BKG_PROCESSOR -------------------------------------------------------------------------------------------------

SiD_bkg_Processor aSiD_bkg_Processor;

SiD_bkg_Processor::SiD_bkg_Processor() : Processor("SiD_bkg_Processor") {
	streamlog_out(DEBUG4) << "Running constructor" << endl;
	// Processor description 
  	// _description = "...";
	gStyle->SetNumberContours(999);
}

void SiD_bkg_Processor::init() { 
	streamlog_out(DEBUG4) << "Running init" << endl;

	File = new TFile("Sim.root","RECREATE"/*"UPDATE"*/,"TestBeam Simulation");
	Tree_MCParticle = new TTree("Tree_MCP","TTree for MCParticles");
	Tree_SiVertexEndcap = new TTree("Tree_SiVertexEndcap","TTree for SimTrackerHits");
	Tree_EcalBarrel = new TTree("Tree_EcalBarrel","TTree for SimCalorimeterHits");
	Tree_MuonEndcap = new TTree("Tree_MuonEndcap","TTree for SimCalorimeterHits");
	Tree_HcalEndcap = new TTree("Tree_HcalEndcap","TTree for SimCalorimeterHits");
	Tree_SiTrackerEndcap = new TTree("Tree_SiTrackerEndcap","TTree for SimTrackerHits");
	Tree_MuonBarrel = new TTree("Tree_MuonBarrel","TTree for SimCalorimeterHits");
	Tree_SiVertexBarrel = new TTree("Tree_SiVertexBarrel","TTree for SimTrackerHits");
	Tree_LumiCal = new TTree("Tree_LumiCal","TTree for SimCalorimeterHits");
	Tree_SiTrackerForward = new TTree("Tree_SiTrackerForward","TTree for SimTrackerHits");
	Tree_HcalBarrel = new TTree("Tree_HcalBarrel","TTree for SimCalorimeterHits");
	Tree_SiTrackerBarrel = new TTree("Tree_SiTrackerBarrel","TTree for SimTrackerHits");
	Tree_EcalEndcap = new TTree("Tree_EcalEndcap","TTree for SimCalorimeterHits");
	Tree_BeamCal = new TTree("Tree_BeamCal","TTree for SimCalorimeterHits");

	registerInputCollectionProcessor(new InputCollectionProcessor_MCParticle_collection(Tree_MCParticle,"InputCollectionName","MCParticle"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimTrackerHit_collection(Tree_SiVertexEndcap,"InputCollectionName2","SiVertexEndcapHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimTrackerHit_collection(Tree_SiTrackerEndcap,"InputCollectionName3","SiTrackerEndcapHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimTrackerHit_collection(Tree_SiVertexBarrel,"InputCollectionName4","SiVertexBarrelHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimTrackerHit_collection(Tree_SiTrackerForward,"InputCollectionName5","SiTrackerForwardHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimTrackerHit_collection(Tree_SiTrackerBarrel,"InputCollectionName6","SiTrackerBarrelHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_EcalBarrel,"InputCollectionName7","EcalBarrelHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_MuonEndcap,"InputCollectionName8","MuonEndcapHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_HcalEndcap,"InputCollectionName9","HcalEndcapHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_MuonBarrel,"InputCollectionName10","MuonBarrelHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_LumiCal,"InputCollectionName11","LumiCalHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_HcalBarrel,"InputCollectionName12","HcalBarrelHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_EcalEndcap,"InputCollectionName13","EcalEndcapHits"));
	registerInputCollectionProcessor(new InputCollectionProcessor_SimCalorimeterHit_collection(Tree_BeamCal,"InputCollectionName14","BeamCalHits"));
/*	registerInputCollectionProcessor(new InputCollectionProcessor_LCGenericObject_collection(Tree,"InputCollectionName7","MCParticleEndPointEnergy"));

	registerOutputCollectionProcessor(new OutputCollectionProcessor_MCParticle_collection("MCParticlePhotonSource"));
*/
  	printParameters() ;
    	_iRun = 0 ;
    	_iEvt = 0 ;

}

void SiD_bkg_Processor::registerInputCollectionProcessor(Collection_Processor_Interface* p){
streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
	registerInputCollection(p->getTypeName(), 
				p->getName(), 
				p->getDescription(),
				p->getCol(),
				p->getDefaultName());

	input_col_processors.push_back(p);
	streamlog_out(DEBUG4)<<"col_processors.back()->getCol()= "<<input_col_processors.back()->getCol()<<std::endl;
}

void SiD_bkg_Processor::registerOutputCollectionProcessor(OutputCollection_Processor_Interface* p){
streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
	registerOutputCollection(p->getTypeName(), 
				p->getName(), 
				p->getDescription(),
				p->getCol(),
				p->getDefaultName());

	output_col_processors.push_back(p);

}

void SiD_bkg_Processor::processRunHeader( LCRunHeader* run) { 
	streamlog_out(DEBUG) << "Running processRunHeader" << endl;

    	_iRun++ ;
} //end processRunHeader


void SiD_bkg_Processor::processEvent( LCEvent * evt ) {   
	streamlog_out(DEBUG0) << "Running processEvent" << endl;
	
	//Process single collections in specific functions
	for(size_t i=0; i<input_col_processors.size(); ++i){
		input_col_processors.at(i)->processCurrentEvent_InputCollection(evt); //in Collection_Processor.h
		streamlog_out(MESSAGE1) << "Process InputCollection #"<< i <<" for current event #" << evt->getEventNumber()  << endl;
	}

/*	for(size_t j=0; j<output_col_processors.size(); ++j){
streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;
		output_col_processors.at(j)->processCurrentEvent_OutputCollection(evt); //in Collection_Processor.h
		streamlog_out(DEBUG4) << "Process OutputCollection #"<< j <<" for current event #" << evt->getEventNumber()  << endl;
	}
*/

    	//-- note: this will not be printed if compiled w/o MARLINstreamlog_out(DEBUG=1 !
    	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
      	  << "   in run:  " << evt->getRunNumber() << std::endl ;

    	_iEvt ++ ;
}


void SiD_bkg_Processor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}// end check()


void SiD_bkg_Processor::end(){ 
streamlog_out(DEBUG0) << "Line " << __LINE__ << endl;


	File->Write();
	File->Close();

    	std::cout << "SiD_bkg_Processor::end()  " << name() 
  	<< " processed " << _iEvt << " events in " << _iRun << " runs "
    	<< std::endl ;

	while(!input_col_processors.empty()){
		streamlog_out(DEBUG4) << "while(!input_col_processors.empty()){"<<std::endl;
		delete input_col_processors.back();
		input_col_processors.pop_back();
	}

}//end end()
//END SID_BKG_PROCESSOR ---------------------------------------------------------------------------------------------
