#ifndef _R3BSTARTRACK_ORDERTS_
#define _R3BSTARTRACK_ORDERTS_

#include "FairTask.h"

using namespace std;

#include <vector>

class TClonesArray;  // just to indicate that this class exists and will be used
class TH1F;          // just to indicate that this class exists and will be used
class TH2F;          // just to indicate that this class exists and will be used
class TTree;

//class R3BTofCalPar;

class R3BStarTrackOrderTS : public FairTask {

 public:

  R3BStarTrackOrderTS();

  R3BStarTrackOrderTS(const char* taskName, Int_t verbose);
 
  virtual ~R3BStarTrackOrderTS();
   
  InitStatus Init(); // Initialisation
 
  void Exec(Option_t *option); // Implementation of event loop
  void FinishEvent(); // End of each event
  void FinishTask(); // End of Task
  void InsertionSort(vector<long long> & v_ts, vector<long long> & v_index);
  void InsertionSort2(vector<long long> & v_ts, vector<long long> & v_block_index, vector<long long> & v_hit_index);


 private:

  // Data members:

  // Input/Output
  TClonesArray* fRawData;
  // Additional data members
  //TH1F *fh_tdc;
  //TH1F *fh_tdc[16]; // we will collect the distribution for each of the 16 channels

// For Run 280-3364 (C target)
  TH1F *TS; 
  TH1F *TSext; 
  TH1F *TS_p; 
  TH1F *TSext_p; 
  TH1F *TS_n; 
  TH1F *TSext_n;
  TH1F *ADC;
  TH1F *TS_TSext_diff;
  TH2F *TS_TSext;
  TH2F *TS_event;
  TH2F *TSext_event;
  TH2F *ADC_TS;
  TH2F *ADC_TS_p;
  TH2F *ADC_TS_n;
  TH2F *Asic_Side;
 
	//TFile* outFile;
	TTree* output_Tree;

	struct struct_entry_sort{
		//long long tm_stp; //reconstructed timestamp (MSB+LSB)
		long long tm_stp; //reconstructed timestamp (MSB+LSB)
		// not used *R3B* ->    long long info; //MBS info data (external timestamp), anything else(?)
		long long tm_stp_ext; //reconstructed external timestamp trigger
		Int_t nhit;  // one hit is one strip hit (word 3)
		Int_t type;
		Int_t hit;
		Int_t ModuleId;
		Int_t Side;
		Int_t AsicId; // new *R3B*
		Int_t StripId;
    // not used *R3B* ->    int type; // QQQ: 0= 20 MeV or 1 GeV (decays), 1= 20 GeV (checked pulser data only in type 0)
              // type>=10: type = info_code+10 (i.e., PAUSE, RESUME, SYNC100, etc...)
		Int_t ADCdata;
		bool sync_flag; // check SYNC100 pulses received for this module
		bool pause_flag; // check Pause signals followed by proper Resume signal: true= SYNC100 paused...
	};
	
	
	  struct_entry_sort s_entry;


	
  Int_t fNevents;
  Int_t fTotalHits; //  

  // Added for handling calibration 
  //R3BTofCalPar* fCal_Par;


 public:
  ClassDef(R3BStarTrackOrderTS,1)

};

#endif
