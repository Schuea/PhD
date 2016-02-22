#ifndef COLLECTION_PROCESSOR_H__
#define COLLECTION_PROCESSOR_H__

#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/LCGenericObject.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>

#ifdef USE_MARLIN
// streamlog include
#include "streamlog/streamlog.h"
#endif
#include "TTree.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#define _USE_MATH_DEFINES

using namespace std;

#define PDGID_GAMMA 22
#define PDGID_POSITRON -11
#define PDGID_ELECTRON 11


template<typename T>
struct TypeNameCollection;

#define REGISTER_COLLECTION(X) template <> struct TypeNameCollection<X> \
    { static const char* TypeName; } ; const char* TypeNameCollection<X>::TypeName = #X

REGISTER_COLLECTION(EVENT::MCParticle);
REGISTER_COLLECTION(EVENT::SimCalorimeterHit);
REGISTER_COLLECTION(EVENT::SimTrackerHit);
REGISTER_COLLECTION(EVENT::LCGenericObject);


class Collection_Processor_Interface{

	public:
	Collection_Processor_Interface (const char* name,
		const char* description,const char* defaultName):name_(name),description_(description),defaultName_(defaultName){}
	
	virtual ~Collection_Processor_Interface (){std::cout<<"~Collection_Processor_Interface"<<std::endl;}	

	virtual void processCurrentEvent_InputCollection (lcio::LCEvent * evt)=0;
	
	
	virtual std::string& getTypeName()=0;
	std::string& getName(){return name_;}
	std::string& getDescription(){return description_;}
	std::string& getCol(){return col_;}
	std::string& getDefaultName(){return defaultName_;}

	private:

	std::string name_,description_,col_,defaultName_;

}; 


template <typename T>
class Collection_Processor_Template : public Collection_Processor_Interface{

	public:

	Collection_Processor_Template(const char* name,
				const char* description,
				const char* defaultName):Collection_Processor_Interface(name,description,defaultName){
		Typename_=std::string(TypeNameCollection<T>::TypeName);
	}

	virtual ~Collection_Processor_Template(){
		std::cout<<"~Collection_Processor_Template()"<<Typename_<<std::endl;
	}

	virtual void processCurrentEvent_InputCollection( lcio::LCEvent * evt){
		EVENT::LCCollection* col=evt->getCollection(getCol());
		
		tree_reset();
			int nMCP = col->getNumberOfElements()  ;
			for(int i=0; i< nMCP ; i++){
				T* p = dynamic_cast<T*>( col->getElementAt(i));		
				fillInHist(p,evt);
			}
	}

	virtual std::string& getTypeName(){
		return  Typename_;
	}

	virtual void fillInHist(T *p, lcio::LCEvent * evt)=0;
	virtual void tree_reset()=0;

	std::string Typename_;
};


class InputCollectionProcessor_MCParticle_collection: public Collection_Processor_Template<EVENT::MCParticle>{

	public:
	
	InputCollectionProcessor_MCParticle_collection(TTree* tree,const char* Name,const char* defaultName):Collection_Processor_Template<EVENT::MCParticle>(Name,"Collection name for MC Particles",defaultName),tree_(tree){		
		z_axis_ = new double[3];
		z_axis_[0] = 0.0;
		z_axis_[1] = 0.0;
		z_axis_[2] = 1.0;
		registerTTree();
		tree_reset();
	}
	~InputCollectionProcessor_MCParticle_collection(){
		std::cout<<"~InputCollectionProcessor_MCParticle_collection"<<std::endl;
	}

	double scalar_product(double a[3], double b[3]){
		double result = 0;
		for (int d=0; d<3; d++) {
			result += a[d]*b[d];
		}
		return result;
	}

	double norm(double a[3]){
		double result = 0;
		double tmp=0;
		for (int d=0; d<3; d++) {
			tmp += a[d]*a[d];
		}
		result = sqrt(tmp);
		return result;
	}
	virtual void fillInHist(EVENT::MCParticle *p, lcio::LCEvent * evt){
		streamlog_out(DEBUG0)<< "Line " << __LINE__  << " File " << __FILE__ << endl;
		
		event_id_= evt->getEventNumber();
		particle_id_ = p->getPDG();
		numberOfParents_ = p->getParents().size();
		numberOfDaughters_ = p->getDaughters().size();
		
		streamlog_out(DEBUG0) << "event_id_: " << event_id_ << endl;	
		streamlog_out(DEBUG0) << "numberOfParents_: " << p->getParents().size() << endl;
		
		if (numberOfParents_==1){
			parent_id_ = p->getParents().at(0)->getPDG();
			parent2_id_=0;
		
			streamlog_out(DEBUG0) << "p->getParents().at(0) : " << p->getParents().at(0) << endl;
			streamlog_out(DEBUG0) << "parent_id_: " << parent_id_ << endl;
			streamlog_out(DEBUG0) << "parent2_id_: " << parent2_id_ << endl;
		}
		else if (numberOfParents_>1){
			parent_id_ = p->getParents().at(0)->getPDG();
		 	parent2_id_ = p->getParents().at(1)->getPDG();
			
			streamlog_out(DEBUG0) << "p->getParents().at(0) : " << p->getParents().at(0) << endl;
			streamlog_out(DEBUG0) << "p->getParents().at(1) : " << p->getParents().at(1) << endl;
			streamlog_out(DEBUG0) << "parent_id_: " << parent_id_ << endl;
			streamlog_out(DEBUG0) << "parent2_id_: " << parent2_id_ << endl;
		}	
		else{
			parent_id_=0;
			parent2_id_=0;
		}
		streamlog_out(DEBUG0) << "parent_id_: " << parent_id_ << endl;
		streamlog_out(DEBUG0) << "parent2_id_: " << parent2_id_ << endl;
	
		createdInSimulation_status_ = p->isCreatedInSimulation(); 
		decayedInTracker_status_ = p->isDecayedInTracker(); 
		decayedInCalo_status_ = p->isDecayedInCalorimeter(); 
		hasLeftDetector_status_ = p->hasLeftDetector();
		createdInContinuousProcess_status_ = p->vertexIsNotEndpointOfParent();
		stopped_status_ = p->isStopped();
		backscatter_status_ = p->isBackscatter();
		
		charge_ = p->getCharge();
		energy_ = p->getEnergy();
		
		momentumx_ = p->getMomentum()[0];
		momentumy_ = p->getMomentum()[1];
		momentumz_ = p->getMomentum()[2];
		momentum_ = sqrt( (p->getMomentum()[0])*(p->getMomentum()[0]) + (p->getMomentum()[1])*(p->getMomentum()[1]) + (p->getMomentum()[2])*(p->getMomentum()[2]) );
				
		reflectionx_ = p->getEndpoint()[0];
		reflectiony_ = p->getEndpoint()[1];
		reflectionz_ = p->getEndpoint()[2];

		vertexx_ = p->getVertex()[0];
		vertexy_ = p->getVertex()[1];
		vertexz_ = p->getVertex()[2];

		
		tree_fill();
		tree_reset();
	}

	virtual void tree_reset(){
		event_id_=0;
		particle_id_=0;
		parent_id_=0;
		parent2_id_=0;
		numberOfParents_=0;		
		numberOfDaughters_=0;		

		createdInSimulation_status_=-1;
		decayedInTracker_status_=-1;
		decayedInCalo_status_=-1;
		hasLeftDetector_status_=-1;
 		createdInContinuousProcess_status_=-1;
		stopped_status_=-1;
		backscatter_status_=-1;

		charge_=0;
		energy_=0;
		momentumx_=0;
		momentumy_=0;
		momentumz_=0;
		momentum_=0;
        	reflectionx_=0;
        	reflectiony_=0;
        	reflectionz_=0;
	       	vertexx_=0;
        	vertexy_=0;
        	vertexz_=0;


	}

	virtual void tree_fill(){
		tree_->Fill();
	}

 	void registerTTree(){

		tree_->Branch("Event_ID",&event_id_,"Event_ID/I");
		tree_->Branch("Particle_ID",&particle_id_,"Particle_ID/I");
		tree_->Branch("Parent_ID",&parent_id_,"Parent_ID/I");
		tree_->Branch("Parent2_ID",&parent2_id_,"Parent2_ID/I");
		tree_->Branch("NumberOfParents",&numberOfParents_,"NumberOfParents/I");
		tree_->Branch("NumberOfDaughters",&numberOfDaughters_,"NumberOfDaughters/I");
		
		tree_->Branch("CreatedInSimulation_Status",&createdInSimulation_status_,"CreatedInSimulation_Status/O");	
		tree_->Branch("DecayedInTracker_Status",&decayedInTracker_status_,"DecayedInTracker_Status/O");	
		tree_->Branch("DecayedInCalo_Status",&decayedInCalo_status_,"DecayedInCalo_Status/O");	
		tree_->Branch("hasLeftDetector_Status",&hasLeftDetector_status_,"hasLeftDetector_Status/O");	
		tree_->Branch("createdInContinuousProcess_Status",&createdInContinuousProcess_status_,"createdInContinuousProcess_Status/O");	
		tree_->Branch("Stopped_Status",&stopped_status_,"Stopped_Status/O");	
		tree_->Branch("Backscatter_Status",&backscatter_status_,"Backscatter_Status/O");	

		tree_->Branch("Charge",&charge_,"Charge/F");
		tree_->Branch("Energy",&energy_,"Energy/D");
		tree_->Branch("Momentumx",&momentumx_,"Momentum/D");
		tree_->Branch("Momentumy",&momentumy_,"Momentum/D");
		tree_->Branch("Momentumz",&momentumz_,"Momentum/D");
		tree_->Branch("Momentum",&momentum_,"Momentum/D");
		tree_->Branch("Reflectionx",&reflectionx_,"Reflectionx/D");
		tree_->Branch("Reflectiony",&reflectiony_,"Reflectiony/D");
		tree_->Branch("Reflectionz",&reflectionz_,"Reflectionz/D");
		tree_->Branch("Vertexx",&vertexx_,"Vertexx/D");
		tree_->Branch("Vertexy",&vertexy_,"Vertexy/D");
		tree_->Branch("Vertexz",&vertexz_,"Vertexz/D");


cout << "Line " << __LINE__  << " File " << __FILE__ << endl;
	}


	private:

	TTree* tree_;

	int event_id_;
	int particle_id_;
	int parent_id_;
	int parent2_id_;
	int numberOfParents_;
	int numberOfDaughters_;
	
	bool createdInSimulation_status_;
	bool decayedInTracker_status_;
	bool decayedInCalo_status_;
	bool hasLeftDetector_status_;
	bool createdInContinuousProcess_status_;
	bool stopped_status_;
	bool backscatter_status_;
        
	float charge_;
	double energy_;
	double momentumx_;
	double momentumy_;
	double momentumz_;
	double momentum_;
	double reflectionx_;
	double reflectiony_;
	double reflectionz_;
	double vertexx_;
	double vertexy_;
	double vertexz_;

	double *z_axis_;
	
	
};


class InputCollectionProcessor_SimCalorimeterHit_collection: public Collection_Processor_Template<EVENT::SimCalorimeterHit>{

	public:
	
	InputCollectionProcessor_SimCalorimeterHit_collection(TTree* tree,const char* Name,const char* defaultName):Collection_Processor_Template<EVENT::SimCalorimeterHit> (Name,"Collection name for SimCalorimeterHits",defaultName),tree_(tree){		
		registerTTree();
		tree_reset();
	}

//____________EDIT:

	virtual void fillInHist(EVENT::SimCalorimeterHit *p, lcio::LCEvent * evt){

					event_id_ = evt->getEventNumber();

					Number_Contr_ = p->getNMCContributions();//Number of contributions to the hit
					//Global coordinates of hit
					HitPosition_x_ = p->getPosition()[0];
					HitPosition_y_ = p->getPosition()[1];
					HitPosition_z_ = p->getPosition()[2];

					HitCellID0_ = p->getCellID0();
					HitCellID1_ = p->getCellID1();

					HitEnergy_ = p->getEnergy();//Total energy [GeV] of hit

					for(int i = 0; i < Number_Contr_; ++i){

									HitContrEnergy_ = p->getEnergyCont(i);//Energy [GeV] of i-th contribution to the hit
									HitContrPDG_ = p->getPDGCont(i);//PDG of shower particle of i-th contribution to the hit
									HitContrTime_ = p->getTimeCont(i);//time [ns] of i-th contribution to the hit

									MCParticle *motherparticle = p->getParticleCont(i);//Mother particle that caused shower responsible for i-th contribution to the hit
									HitMotherLCIO_id_ = motherparticle->id();
									HitMotherCreationTime_ = motherparticle->getTime();

									HitMotherVertex_x_ = motherparticle->getVertex()[0];
									HitMotherVertex_y_ = motherparticle->getVertex()[1];
									HitMotherVertex_z_ = motherparticle->getVertex()[2];
									HitMotherEndpoint_x_ = motherparticle->getEndpoint()[0];
									HitMotherEndpoint_y_ = motherparticle->getEndpoint()[1];
									HitMotherEndpoint_z_ = motherparticle->getEndpoint()[2];

									HitMotherMomentum_x_ = motherparticle->getMomentum()[0];
									HitMotherMomentum_y_ = motherparticle->getMomentum()[1];
									HitMotherMomentum_z_ = motherparticle->getMomentum()[2];

									HitMotherParticle_id_ = motherparticle->getPDG();
									HitMotherParticleEnergy_ = motherparticle->getEnergy();
									HitMotherParticleCharge_ = motherparticle->getCharge();

									tree_->Fill();
									tree_reset();
					}

	}

	virtual void tree_reset(){

					event_id_ = -999;

					HitCellID0_ = -999;
					HitCellID1_ = -999;

					HitEnergy_ = -999;
					HitPosition_x_ = -999;
					HitPosition_y_ = -999;
					HitPosition_z_ = -999;

					Number_Contr_ = -999;
					HitContrEnergy_ = -999; 
					HitContrPDG_ = -999;
					HitContrTime_ = -999;

					HitMotherLCIO_id_ = -999;
					HitMotherCreationTime_ = -999;

					HitMotherVertex_x_ = -999;
					HitMotherVertex_y_ = -999;
					HitMotherVertex_z_ = -999;

					HitMotherEndpoint_x_ = -999;
					HitMotherEndpoint_y_ = -999;
					HitMotherEndpoint_z_ = -999;

					HitMotherMomentum_x_ = -999;
					HitMotherMomentum_y_ = -999;
					HitMotherMomentum_z_ = -999;

					HitMotherParticle_id_ = -999;
					HitMotherParticleCharge_ = -999;
					HitMotherParticleEnergy_ = -999;

	}

 	void registerTTree(){

		tree_->Branch("event_id",&event_id_,"event_id/I");

		tree_->Branch("HitCellID0",&HitCellID0_,"HitCellID0/I");
		tree_->Branch("HitCellID1",&HitCellID1_,"HitCellID1/I");

		tree_->Branch("HitEnergy",&HitEnergy_,"HitEnergy/F");
	
		tree_->Branch("HitPosition_x",&HitPosition_x_,"HitPosition_x/F");
		tree_->Branch("HitPosition_y",&HitPosition_y_,"HitPosition_y/F");
		tree_->Branch("HitPosition_z",&HitPosition_z_,"HitPosition_z/F");

		tree_->Branch("Number_Contr",&Number_Contr_,"Number_Contr/I");
		tree_->Branch("HitContrEnergy",&HitContrEnergy_,"HitContrEnergy/F"); 
		tree_->Branch("HitContrPDG",&HitContrPDG_,"HitContrPDG/I");
		tree_->Branch("HitContrTime",&HitContrTime_,"HitContrTime/F");

		tree_->Branch("HitMotherLCIO_id",&HitMotherLCIO_id_,"HitMotherLCIO_id/I");
		tree_->Branch("HitMotherCreationTime",&HitMotherCreationTime_,"HitMotherCreationTime/F");

		tree_->Branch("HitMotherVertex_x",&HitMotherVertex_x_,"HitMotherVertex_x/D");
		tree_->Branch("HitMotherVertex_y",&HitMotherVertex_y_,"HitMotherVertex_y/D");
		tree_->Branch("HitMotherVertex_z",&HitMotherVertex_z_,"HitMotherVertex_z/D");

		tree_->Branch("HitMotherEndpoint_x",&HitMotherEndpoint_x_,"HitMotherEndpoint_x/D");
		tree_->Branch("HitMotherEndpoint_y",&HitMotherEndpoint_y_,"HitMotherEndpoint_y/D");
		tree_->Branch("HitMotherEndpoint_z",&HitMotherEndpoint_z_,"HitMotherEndpoint_z/D");

		tree_->Branch("HitMotherMomentum_x",&HitMotherMomentum_x_,"HitMotherMomentum_x/D");
		tree_->Branch("HitMotherMomentum_y",&HitMotherMomentum_y_,"HitMotherMomentum_y/D");
		tree_->Branch("HitMotherMomentum_z",&HitMotherMomentum_z_,"HitMotherMomentum_z/D");

		tree_->Branch("HitMotherParticle_ID",&HitMotherParticle_id_,"HitMotherParticle_ID/I");
		tree_->Branch("HitMotherCharge",&HitMotherParticleCharge_,"HitMotherCharge/F");
		tree_->Branch("HitMotherParticleEnergy",&HitMotherParticleEnergy_,"HitMotherParticleEnergy/D");
		

	}

	private:

	TTree* tree_;

	int event_id_;

	int HitCellID0_;
	int HitCellID1_;
	
	float HitEnergy_;

	float HitPosition_x_;
	float HitPosition_y_;
	float HitPosition_z_;

	int   Number_Contr_;
	float HitContrEnergy_; 
	int   HitContrPDG_;
	float HitContrTime_;

	int		 HitMotherLCIO_id_;
	float	 HitMotherCreationTime_;

	double HitMotherVertex_x_;
	double HitMotherVertex_y_;
	double HitMotherVertex_z_;
	
	double HitMotherEndpoint_x_;
	double HitMotherEndpoint_y_;
	double HitMotherEndpoint_z_;
	
	double HitMotherMomentum_x_;
	double HitMotherMomentum_y_;
	double HitMotherMomentum_z_;

	int		 HitMotherParticle_id_;
	float  HitMotherParticleCharge_;
  double  HitMotherParticleEnergy_;
        
//_____________________________

};

class InputCollectionProcessor_SimTrackerHit_collection: public Collection_Processor_Template<EVENT::SimTrackerHit>{

	public:
	
	InputCollectionProcessor_SimTrackerHit_collection(TTree* tree,const char* Name,const char* defaultName):Collection_Processor_Template<EVENT::SimTrackerHit> (Name,"Collection name for SimTrackerHits",defaultName),tree_(tree){		
	//cout << Name << endl;	
		registerTTree();
		tree_reset();
	}

//____________EDIT:

	virtual void fillInHist(EVENT::SimTrackerHit *p, lcio::LCEvent * evt){
		event_id_= evt->getEventNumber();
		
		HitPosition_x_ = p->getPosition()[0];
		HitPosition_y_ = p->getPosition()[1];
		HitPosition_z_ = p->getPosition()[2];

		HitCellID_ = p->getCellID();
		HitdEdx_ = p->getdEdx();//dE/dx of the hit in GeV
		
		HitTime_ = p->getTime();//time of the hit in ns

		HitMomentum_x_ = p->getMomentum()[0];
		HitMomentum_y_ = p->getMomentum()[1];
		HitMomentum_z_ = p->getMomentum()[2];

		MCParticle *hitparticle = p->getMCParticle();
		HitParticleLCIO_id_ = hitparticle->id();
		HitParticleCreationTime_ = hitparticle->getTime();

		HitParticleVertex_x_ = hitparticle->getVertex()[0];
		HitParticleVertex_y_ = hitparticle->getVertex()[1];
		HitParticleVertex_z_ = hitparticle->getVertex()[2];
		HitParticleEndpoint_x_ = hitparticle->getEndpoint()[0];
		HitParticleEndpoint_y_ = hitparticle->getEndpoint()[1];
		HitParticleEndpoint_z_ = hitparticle->getEndpoint()[2];

		HitParticleMomentum_x_ = hitparticle->getMomentum()[0];
		HitParticleMomentum_y_ = hitparticle->getMomentum()[1];
		HitParticleMomentum_z_ = hitparticle->getMomentum()[2];

		HitParticle_id_ = hitparticle->getPDG();
		HitParticleEnergy_ = hitparticle->getEnergy();
		HitParticleCharge_ = hitparticle->getCharge();

		tree_->Fill();
		tree_reset();
	}

	virtual void tree_reset(){

					event_id_ = -999;

					HitCellID_ = -999;

					HitTime_ = -999;
					HitdEdx_ = -999;
					HitPosition_x_ = -999;
					HitPosition_y_ = -999;
					HitPosition_z_ = -999;

					HitMomentum_x_ = -999;
					HitMomentum_y_ = -999;
					HitMomentum_z_ = -999;

					HitParticleLCIO_id_ = -999;
					HitParticleCreationTime_ = -999;

					HitParticleVertex_x_ = -999;
					HitParticleVertex_y_ = -999;
					HitParticleVertex_z_ = -999;

					HitParticleEndpoint_x_ = -999;
					HitParticleEndpoint_y_ = -999;
					HitParticleEndpoint_z_ = -999;

					HitParticleMomentum_x_ = -999;
					HitParticleMomentum_y_ = -999;
					HitParticleMomentum_z_ = -999;

					HitParticle_id_ = -999;
					HitParticleCharge_ = -999;
					HitParticleEnergy_ = -999;

	}

 	void registerTTree(){
		tree_->Branch("event_id",&event_id_,"event_id/I");
		
		tree_->Branch("HitCellID",&HitCellID_,"HitCellID/I");
		tree_->Branch("HitdEdx",&HitdEdx_,"HitdEdx/F");
		tree_->Branch("HitTime",&HitTime_,"HitTime/F");
	
		tree_->Branch("HitPosition_x",&HitPosition_x_,"HitPosition_x/D");
		tree_->Branch("HitPosition_y",&HitPosition_y_,"HitPosition_y/D");
		tree_->Branch("HitPosition_z",&HitPosition_z_,"HitPosition_z/D");

		tree_->Branch("HitMomentum_x",&HitMomentum_x_,"HitMomentum_x/F");
		tree_->Branch("HitMomentum_y",&HitMomentum_y_,"HitMomentum_y/F");
		tree_->Branch("HitMomentum_z",&HitMomentum_z_,"HitMomentum_z/F");

		tree_->Branch("HitParticleLCIO_ID",&HitParticleLCIO_id_,"HitParticleLCIO_ID/I");
		tree_->Branch("HitParticleCreationTime",&HitParticleCreationTime_,"HitParticleCreationTime/F");
		tree_->Branch("HitParticleVertex_x",&HitParticleVertex_x_,"HitParicleVertex_x/D");
		tree_->Branch("HitParticleVertex_y",&HitParticleVertex_y_,"HitParicleVertex_y/D");
		tree_->Branch("HitParticleVertex_z",&HitParticleVertex_z_,"HitParicleVertex_z/D");

		tree_->Branch("HitParticleEndpoint_x",&HitParticleEndpoint_x_,"HitParticleEndpoint_x/D");
		tree_->Branch("HitParticleEndpoint_y",&HitParticleEndpoint_y_,"HitParticleEndpoint_y/D");
		tree_->Branch("HitParticleEndpoint_z",&HitParticleEndpoint_z_,"HitParticleEndpoint_z/D");

		tree_->Branch("HitParticleMomentum_x",&HitParticleMomentum_x_,"HitParticleMomentum_x/D");
		tree_->Branch("HitParticleMomentum_y",&HitParticleMomentum_y_,"HitParticleMomentum_y/D");
		tree_->Branch("HitParticleMomentum_z",&HitParticleMomentum_z_,"HitParticleMomentum_z/D");

		tree_->Branch("HitParticle_ID",&HitParticle_id_,"HitParticle_ID/I");
		tree_->Branch("HitParticleCharge",&HitParticleCharge_,"HitParticleCharge/F");
		tree_->Branch("HitParticleEnergy",&HitParticleEnergy_,"HitParticleEnergy/D");
	}

	private:

	TTree* tree_;

	int event_id_;

	int		HitCellID_;
	float HitdEdx_;
	float HitTime_;

	double HitPosition_x_;
	double HitPosition_y_;
	double HitPosition_z_;

	float HitMomentum_x_;
	float HitMomentum_y_;
	float HitMomentum_z_;

	int    HitParticleLCIO_id_;
	float	 HitParticleCreationTime_;

	double HitParticleVertex_x_;
	double HitParticleVertex_y_;
	double HitParticleVertex_z_;

	double HitParticleEndpoint_x_;
	double HitParticleEndpoint_y_;
	double HitParticleEndpoint_z_;

	double HitParticleMomentum_x_;
	double HitParticleMomentum_y_;
	double HitParticleMomentum_z_;

	int HitParticle_id_;
	float HitParticleCharge_;
	double HitParticleEnergy_;

	//_____________________________

};


class InputCollectionProcessor_LCGenericObject_collection: public Collection_Processor_Template<EVENT::LCGenericObject>{

	public:
	
	InputCollectionProcessor_LCGenericObject_collection(TTree* tree,const char* Name,const char* defaultName):Collection_Processor_Template<EVENT::LCGenericObject>(Name,"Collection name for LCGenericObject",defaultName),tree_(tree){		
		registerTTree();
		tree_reset();
	}

//____________EDIT:

	virtual void fillInHist(EVENT::LCGenericObject *p, lcio::LCEvent * evt){
/*
		event_id_= evt->getEventNumber();
		particle_id_ = p->getPDG();
		
		energy_ = p->getEnergy();
				
		reflectionx_ = p->getEndpoint()[0];
		reflectiony_ = p->getEndpoint()[1];
		reflectionz_ = p->getEndpoint()[2];
		tree_fill();
*/
	}

	virtual void tree_reset(){
/*
		event_id_=0;
		particle_id_=0;

		energy_=0;
        	reflectionx_=0;
        	reflectiony_=0;
        	reflectionz_=0;
*/
	}

	virtual void tree_fill(){
		tree_->Fill();
	}

 	void registerTTree(){
/*
		tree_->Branch("Event_ID",&event_id_,"Event_ID/I");
		tree_->Branch("Particle_ID",&particle_id_,"Particle_ID/I");
	
		tree_->Branch("Energy",&energy_,"Energy/F");
		tree_->Branch("Reflectionx",&reflectionx_,"Reflectionx/F");
		tree_->Branch("Reflectiony",&reflectiony_,"Reflectiony/F");
		tree_->Branch("Reflectionz",&reflectionz_,"Reflectionz/F");
*/
	}

	private:

	TTree* tree_;
/*
	int event_id_;
	int particle_id_;

        float energy_;
	float reflectionx_;
	float reflectiony_;
	float reflectionz_;
*/
//_____________________________

};

#endif
