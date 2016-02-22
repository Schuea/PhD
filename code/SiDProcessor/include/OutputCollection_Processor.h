#ifndef OUTPUTCOLLECTION_PROCESSOR_H__
#define OUTPUTCOLLECTION_PROCESSOR_H__

#include "marlin/Processor.h"
#include "marlin/Global.h"
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCGenericObject.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

#include <IMPL/LCCollectionVec.h>

#include "TTree.h"
#include <iostream>
#include <algorithm>

#ifdef USE_MARLIN
// streamlog include
#include "streamlog/streamlog.h"
#endif

#define PDGID_GAMMA 22
#define PDGID_POSITRON -11
#define PDGID_ELECTRON 11


template<typename T>
struct OutputTypeNameCollection;

#define REGISTER_OutputCOLLECTION(X) template <> struct OutputTypeNameCollection<X> \
    { static const char* OutputTypeName; } ; const char* OutputTypeNameCollection<X>::OutputTypeName = #X

REGISTER_OutputCOLLECTION(EVENT::MCParticle);
REGISTER_OutputCOLLECTION(EVENT::SimCalorimeterHit);
REGISTER_OutputCOLLECTION(EVENT::LCGenericObject);


class OutputCollection_Processor_Interface{

	public:

	OutputCollection_Processor_Interface (const char* name,
		const char* description,const char* defaultName):name_(name),description_(description),defaultName_(defaultName){}
	
	virtual ~OutputCollection_Processor_Interface (){std::cout<<"~OutputCollection_Processor_Interface"<<std::endl;}	

	virtual void processCurrentEvent_OutputCollection (lcio::LCEvent * evt)=0;
	
	virtual std::string& getTypeName()=0;
	std::string& getName(){return name_;}
	std::string& getDescription(){return description_;}
	std::string& getCol(){return col_;}
	std::string& getDefaultName(){return defaultName_;}

	private:

	std::string name_,description_,col_,defaultName_;

}; 


template <typename T>
class OutputCollection_Processor_Template : public OutputCollection_Processor_Interface{
	
	private:
		IMPL::LCCollectionVec* outputcol;

	public:

	OutputCollection_Processor_Template(const char* name,
				const char* description,
				const char* defaultName):OutputCollection_Processor_Interface(name,description,defaultName){
		Typename_=std::string(OutputTypeNameCollection<T>::OutputTypeName);
		outputcol = new IMPL::LCCollectionVec(Typename_);
	}

	virtual ~OutputCollection_Processor_Template(){
    		streamlog_out(DEBUG)<<"~OutputCollection_Processor_Template()"<<Typename_<<std::endl;
	}

	virtual void processCurrentEvent_OutputCollection( lcio::LCEvent * evt){

		//IMPL::LCCollectionVec* outputcol = new IMPL::LCCollectionVec(Typename_);
		std::string outputname=getCol();
		streamlog_out(DEBUG4) << "New OutputCollection created..." << std::endl;

                if (std::find(evt->getCollectionNames()->begin(),evt->getCollectionNames()->end(),outputname)==evt->getCollectionNames()->end()){
                        //Doesn't exist so add collection and add elements
			evt->addCollection(outputcol,outputname);
			streamlog_out(DEBUG4)<< "outputcol "<< outputcol << ", outputname " << outputname << std::endl;
                } else{ 
			std::cout<< "Collection " << outputname <<" already exists!"<<std::endl;
		}
		addElementtoOutputCollection(evt, outputname);
		streamlog_out(DEBUG4) << "All wanted elements of this event added." << std::endl;
	}

	virtual std::string& getTypeName(){
		return  Typename_;
	}

	virtual void addElementtoOutputCollection(lcio::LCEvent * evt, std::string colName)=0;

	std::string Typename_;

};


class OutputCollectionProcessor_MCParticle_collection: public OutputCollection_Processor_Template<EVENT::MCParticle>{
	private:
		EVENT::LCCollection* col;

	public:

		OutputCollectionProcessor_MCParticle_collection(const char* defaultName):OutputCollection_Processor_Template<EVENT::MCParticle>("OutputCollectionName", "OutputCollection name for MC Particles",defaultName){
		}

		~OutputCollectionProcessor_MCParticle_collection(){
			streamlog_out(DEBUG)<<"~OutputCollectionProcessor_MCParticle_collection"<<std::endl;
		}

		virtual void addElementtoOutputCollection(lcio::LCEvent * evt, std::string colName){

			for(StringVec::const_iterator it = evt->getCollectionNames()->begin(); it != evt->getCollectionNames()->end() ; ++it){
				//std::cout<< "StringVec::const_iterator it = evt->getCollectionNames()" << *it << std::endl;

				if(*it == "MCParticle"){
				streamlog_out(DEBUG4)<< "it = " << *it << std::endl;
				col=evt->getCollection(*it);
				int nMCP = col->getNumberOfElements();
				streamlog_out(DEBUG)<<nMCP<<std::endl;
			
				for(int i=0; i< nMCP ; i++){
					streamlog_out(DEBUG4) << "Going through elements in this collection..." << std::endl;

					MCParticle* p = dynamic_cast<MCParticle*>( col->getElementAt(i));
	
					if(p->getPDG()==PDGID_GAMMA) {
						streamlog_out(DEBUG4)<< "Photon was added!" <<std::endl;
						evt->getCollection(colName)->addElement(p);
						streamlog_out(DEBUG4)<<"evt->getCollection(colName)->getNumberOfElements()= ";
						streamlog_out(DEBUG4)<< evt->getCollection(colName)->getNumberOfElements() <<std::endl;
					}else continue;

				}streamlog_out(DEBUG4) << "Went through all elements in this collection!" << std::endl;

			}

		}streamlog_out(DEBUG4) << "Went through all collections in this event!" << std::endl;
	}

};

#endif
