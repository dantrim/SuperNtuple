// rjigsawAna.cxx


// std
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <math.h>

// ROOT
#include "TChain.h"
#include "TVectorD.h"
#include "TRandom.h"

// SusyNtuple
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/TriggerTools.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/PhysicsTools.h"

// RestFrames
#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace sflow;
using namespace RestFrames;

const string analysis_name = "rjigsawAna_WWBB";

////////////////////////////////////////
// Function Prototypes
////////////////////////////////////////
void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_, string& suffix_name_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys_, bool& split_sumw, bool dbg);


////////////////////////////////////////
// MAIN
////////////////////////////////////////
int main(int argc, char* argv[])
{
    ///////////////////////
    // READ IN
    ///////////////////////
    int n_skip_ = 0;
    int num_events_ = -1;
    bool m_dbg = false;
    string sample_;
    string suffix_name_ = "";
    bool do_sumw_split = false;
    SuperflowRunMode run_mode_ = SuperflowRunMode::nominal;
    SusyNtSys nt_sys_ = NtSys::NOM;

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    read_options(argc, argv, chain, n_skip_, num_events_, sample_, suffix_name_, run_mode_, nt_sys_, do_sumw_split, m_dbg);

    ////////////////////////////////////////////////////
    // Construct and configure the Superflow object
    ////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow();
    cutflow->setAnaName(analysis_name);
    cutflow->setAnaType(AnalysisType::Ana_WWBB);
    cutflow->setLumi(1000); // 1/fb
    cutflow->setSampleName(sample_);
    cutflow->setRunMode(run_mode_);
    cutflow->setCountWeights(true);
    cutflow->setChain(chain);
    if(suffix_name_!="")
        cutflow->setFileSuffix(suffix_name_);
    if(do_sumw_split) {
        string sumw_file = "./n0232val/sumw_file.txt";
        //string sumw_file = "/data/uclhc/uci/user/dantrim/n0229val/sumw_file.txt"; 
        cout << analysis_name << "    Reading sumw for sample from file: " << sumw_file << endl; 
        cutflow->setUseSumwFile(sumw_file);
    }
    cutflow->nttools().initTriggerTool(ChainHelper::firstFile(sample_, m_dbg));

    // print some useful
    cout << analysis_name << "    Total Entries    : " << chain->GetEntries() << endl;
    if(num_events_ > 0) {
    cout << analysis_name << "    Process Entries  : " << num_events_ << endl;
    } else {
    cout << analysis_name << "    Process Entries  : " << chain->GetEntries() << endl;
    }

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    //
    // Superflow methods [BEGIN]
    //
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    *cutflow << CutName("read in ") << [](Superlink* sl) -> bool { return true; };

    ////////////////////////////////////////////////////
    // Cleaning cuts
    ////////////////////////////////////////////////////
    int cutflags = 0;
    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];
        return (sl->tools->passGRL(cutflags));
    };
    *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passLarErr(cutflags));
    };
    *cutflow << CutName("Tile Error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTileErr(cutflags));
    };
    *cutflow << CutName("SCT error") << [&](Superlink* sl) -> bool {
        return (sl->tools->passSCTErr(cutflags));
    };
    *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passTTC(cutflags));
    };
    *cutflow << CutName("pass Good Vertex") << [&](Superlink * sl) -> bool {
        return (sl->tools->passGoodVtx(cutflags));
    };
    *cutflow << CutName("pass bad muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passBadMuon(sl->preMuons));
    };
    *cutflow << CutName("pass cosmic muon veto") << [&](Superlink* sl) -> bool {
        return (sl->tools->passCosmicMuon(sl->baseMuons));
    };
    *cutflow << CutName("pass jet cleaning") << [&](Superlink* sl) -> bool {
        return (sl->tools->passJetCleaning(sl->baseJets));
    };

//    *cutflow << CutName("pass met cleaning") << [&](Superlink* sl) -> bool {
//        return (sl->tools->passMetCleaning(sl->met));
//    };

    ///////////////////////////////////////////////////
    // Analysis Cuts
    ///////////////////////////////////////////////////
    *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    ///*
    *cutflow << CutName("lepton pTs > 20 GeV") << [](Superlink* sl) -> bool {
        return ( (sl->leptons->at(0)->Pt()>20) && (sl->leptons->at(1)->Pt()>20) );
    };
    //*/

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
    };

    ///*
    *cutflow << CutName("mll > 20 GeV") << [](Superlink* sl) -> bool {
        return ( (*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 20. );
    };
    //*/

    *cutflow << CutName("veto SF Z-window (within 10 GeV)") << [](Superlink* sl) -> bool {
        bool pass = true;
        bool isSF = false;
        if((sl->leptons->size()==2 && (sl->electrons->size()==2 || sl->muons->size()==2))) isSF = true;
        if(isSF) {
            double mll = (*sl->leptons->at(0) + *sl->leptons->at(1)).M();
            if( fabs(mll-91.2) < 20. ) pass = false;
        }
        return pass;
    };


    /*
    *cutflow << CutName("pass trigger requirement") << [&](Superlink* sl) -> bool {
        int year = sl->nt->evt()->treatAsYear;
        bool pass = false;

        bool passes_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
        bool passes_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
        bool passes_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
        bool passes_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
        bool passes_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
        bool passes_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14");

        bool trig_pass_2015 = (passes_mu18_mu8noL1 || passes_e17_lhloose_mu14 || passes_2e12_lhloose_L12EM10VH);
        bool trig_pass_2016 = (passes_2e17_lhvloose_nod0 || passes_mu22_mu8noL1 || passes_e17_lhloose_nod0_mu14);
        if( (year==2015 && trig_pass_2015==1) || (year==2016 && trig_pass_2016==1) ) {
            pass = true;
        }
        return pass;
    };
    */
    
    //*cutflow << CutName("isSF") << [&](Superlink* sl) -> bool {
    //    bool isSF = false;
    //    bool zveto = false;
    //    if((sl->leptons->size()==2 && (sl->electrons->size()==2 || sl->muons->size()==2))) isSF = true;
    //    if(isSF) {
    //        double mll = (*sl->leptons->at(0) + *sl->leptons->at(1)).M();
    //        if(fabs(mll-91.2)>10.) zveto = true;
    //    }
    //    if(isSF && zveto) return true;
    //    else { return false; }
    //};
    //*cutflow << CutName("isDF") << [&](Superlink* sl) -> bool {
    //    bool isDF = false;
    //    if((sl->leptons->size()==2 && sl->electrons->size()==1 && sl->muons->size()==1)) isDF = true;
    //    if(isDF) return true;
    //    else { return false; }
    //};

    ///////////////////////////////////////////////////
    // Ntuple Setup
    ///////////////////////////////////////////////////


    // TRIGGERS
    bool pass_mu18_mu8noL1;
    bool pass_mu20_mu8noL1;
    bool pass_e17_lhloose_mu14;
    bool pass_2e12_lhloose_L12EM10VH;
    bool pass_2e15_lhvloose_L12EM13VH;
    // updated
    bool pass_2e17_lhvloose_nod0;
    bool pass_mu22_mu8noL1;
    bool pass_e17_lhloose_nod0_mu14;

    // for eff
    bool p_mu20;
    bool p_mu20_iloose_L1MU15;
    bool p_mu22;
    bool p_mu24;
    bool p_mu24_iloose;
    bool p_mu24_imedium;
    bool p_mu24_iloose_L1MU15;
    bool p_mu26_imedium;
    bool p_mu40;
    bool p_mu50;
    bool p_mu60_0eta105_msonly;
    bool p_mu18_mu8noL1;
    bool p_mu20_mu8noL1;
    bool p_mu22_mu8noL1;
    bool p_mu24_mu8noL1;
    bool p_e12_lhloose;
    bool p_e12_lhloose_L1EM10VH;
    bool p_e15_lhloose_L1EM13VH;
    bool p_e17_lhloose;
    bool p_e17_lhmedium;
    bool p_e24_lhmedium_iloose_L1EM18VH;
    bool p_e24_lhmedium_iloose_L1EM20VH;
    bool p_e24_lhmedium_L1EM20VHI;
    bool p_e24_lhtight_iloose;
    bool p_e24_lhtight_nod0_iloose;
    bool p_e24_lhmedium_nod0_iloose_L1EM20VH;
    bool p_e24_lhtight_ivarloose;
    bool p_e24_lhtight_nod0_ivarloose;
    bool p_e26_lhmedium_L1EM22VHI;
    bool p_e26_lhtight_iloose;
    bool p_e26_lhtight_nod0_iloose;
    bool p_e26_lhtight_ivarloose;
    bool p_e26_lhtight_nod0_ivarloose;
    bool p_e60_lhmedium;
    bool p_e60_lhmedium_nod0;
    bool p_2e12_lhloose_L12EM10VH;
    bool p_2e15_lhloose_L12EM13VH;
    bool p_2e15_lhvloose_L12EM13VH;
    bool p_2e15_lhvloose_nod0_L12EM13VH;
    bool p_2e17_lhloose;
    bool p_2e17_lhvloose;
    bool p_2e17_lhvloose_nod0;
    bool p_e17_lhloose_mu14;
    bool p_e17_lhloose_nod0_mu14;
    bool p_e7_lhmedium_mu24;
    bool p_e7_lhmedium_nod0_mu24;
    bool p_e26_lhmedium_L1EM22VHI_mu8noL1;
    bool p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;
    bool p_e24_lhmedium_L1EM20VHI_mu8noL1;


    *cutflow << [&](Superlink* sl, var_void*) {

        // S2L

        pass_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
        pass_mu20_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_mu8noL1");
        pass_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
        pass_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
        pass_2e15_lhvloose_L12EM13VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e15_lhvloose_L12EM13VH");
        pass_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
        pass_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
        pass_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14"); 

        // eff
        p_mu20 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20");
        p_mu20_iloose_L1MU15 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_iloose_L1MU15");
        p_mu22 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22");
        p_mu24 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu24");
        p_mu24_iloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu24_iloose");
        p_mu24_imedium = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu24_imedium");
        p_mu24_iloose_L1MU15 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu24_iloose_L1MU15");
        p_mu26_imedium = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu26_imedium"); 
        p_mu40 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu40");
        p_mu50 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu50");
        p_mu60_0eta105_msonly = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu60_0eta105_msonly");
        p_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
        p_mu20_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_mu8noL1");
        p_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
        p_mu24_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu24_mu8noL1");
        p_e12_lhloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e12_lhloose");
        p_e12_lhloose_L1EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e12_lhloose_L1EM10VH");
        p_e15_lhloose_L1EM13VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e15_lhloose_L1EM13VH");
        p_e17_lhloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose");
        p_e17_lhmedium = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhmedium");
        p_e24_lhmedium_iloose_L1EM18VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_iloose_L1EM18VH");
        p_e24_lhmedium_iloose_L1EM20VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_iloose_L1EM20VH");
        p_e24_lhmedium_L1EM20VHI = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_L1EM20VHI");
        p_e24_lhtight_iloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhtight_iloose");
        p_e24_lhtight_nod0_iloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhtight_nod0_iloose");
        p_e24_lhmedium_nod0_iloose_L1EM20VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_nod0_iloose_L1EM20VH");
        p_e24_lhtight_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhtight_ivarloose");
        //p_e24_lhtight_nod0_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhtight_nod0_ivarloose");
        p_e24_lhtight_nod0_ivarloose = false;

        p_e26_lhmedium_L1EM22VHI = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhmedium_L1EM22VHI");
        p_e26_lhtight_iloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhtight_iloose");
        p_e26_lhtight_nod0_iloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhtight_nod0_iloose");
        p_e26_lhtight_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhtight_ivarloose");
        //p_e26_lhtight_nod0_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhtight_nod0_ivarloose");
        p_e26_lhtight_nod0_ivarloose = false;
    
        p_e60_lhmedium = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e60_lhmedium");
        p_e60_lhmedium_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e60_lhmedium_nod0");
        p_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
        p_2e15_lhloose_L12EM13VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e15_lhloose_L12EM13VH");
        p_2e17_lhloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhloose");
        p_2e17_lhvloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose");
        p_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
        p_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
        p_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14");
        p_e7_lhmedium_mu24 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e7_lhmedium_mu24");
        p_e7_lhmedium_nod0_mu24 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e7_lhmedium_nod0_mu24");
        p_e26_lhmedium_L1EM22VHI_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhmedium_L1EM22VHI_mu8noL1");
        p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1");
        p_e24_lhmedium_L1EM20VHI_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_L1EM20VHI_mu8noL1");
        
    };
    *cutflow << NewVar("pass mu18_mu8noL1"); {
        *cutflow << HFTname("trig_mu18_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu18_mu8noL1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu20_mu8noL1"); {
        *cutflow << HFTname("trig_mu20_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu20_mu8noL1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e17_lhloose_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_e17_lhloose_mu14;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e12_lhloose_L12EM10VH"); {
        *cutflow << HFTname("trig_2e12_lhloose_L12EM10VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e12_lhloose_L12EM10VH;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e15_lhvloose_L12EM13VH"); {
        *cutflow << HFTname("trig_2e15_lhvloose_L12EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e15_lhvloose_L12EM13VH;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass OR 2015 trigger set"); {
        *cutflow << HFTname("trig_pass2015");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_mu18_mu8noL1 || pass_e17_lhloose_mu14 || pass_2e12_lhloose_L12EM10VH);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass OR 2016 trigger set"); {
        *cutflow << HFTname("trig_pass2016");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_mu20_mu8noL1 || pass_e17_lhloose_mu14 || pass_2e15_lhvloose_L12EM13VH);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass HLT_2e17_lhvloose_nod0"); {
        *cutflow << HFTname("trig_2e17_lhvloose_nod0");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_2e17_lhvloose_nod0;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass mu22_mu8noL1"); {
        *cutflow << HFTname("trig_mu22_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_mu22_mu8noL1;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass HLT_e17_lhloose_nod0_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_nod0_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return pass_e17_lhloose_nod0_mu14;
        };
        *cutflow << SaveVar();
    } 
    *cutflow << NewVar("pass OR UPDATED trigger set"); {
        *cutflow << HFTname("trig_pass2016update");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return (pass_2e17_lhvloose_nod0 || pass_mu22_mu8noL1 || pass_e17_lhloose_nod0_mu14); 
        };
        *cutflow << SaveVar();
    }


    // eff trig cuts
    *cutflow << NewVar("pass_mu20"); {
        *cutflow << HFTname("trig_mu20");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu20;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu20_iloose_L1MU15"); {
        *cutflow << HFTname("trig_mu20_iloose_L1MU15");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu20_iloose_L1MU15;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu22"); {
        *cutflow << HFTname("trig_mu22");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu22;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu24"); {
        *cutflow << HFTname("trig_mu24");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu24;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu24_iloose"); {
        *cutflow << HFTname("trig_mu24_iloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu24_iloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu24_imedium"); {
        *cutflow << HFTname("trig_mu24_imedium");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu24_imedium;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu24_iloose_L1MU15"); {
        *cutflow << HFTname("trig_mu24_iloose_L1MU15");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu24_iloose_L1MU15;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu26_imedium"); {
        *cutflow << HFTname("trig_mu26_imedium");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu26_imedium;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu40"); {
        *cutflow << HFTname("trig_mu40");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu40;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu50"); {
        *cutflow << HFTname("trig_mu50");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu50;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu60_0eta105_msonly"); {
        *cutflow << HFTname("trig_mu60_0eta105_msonly");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu60_0eta105_msonly;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu18_mu8noL1"); {
        *cutflow << HFTname("trig_mu18_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu18_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu20_mu8noL1"); {
        *cutflow << HFTname("trig_mu20_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu20_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu22_mu8noL1"); {
        *cutflow << HFTname("trig_mu22_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu22_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_mu24_mu8noL1"); {
        *cutflow << HFTname("trig_mu24_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_mu24_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e12_lhloose"); {
        *cutflow << HFTname("trig_e12_lhloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e12_lhloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e12_lhloose_L1EM10VH"); {
        *cutflow << HFTname("trig_e12_lhloose_L1EM10VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e12_lhloose_L1EM10VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e15_lhloose_L1EM13VH"); {
        *cutflow << HFTname("trig_e15_lhloose_L1EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e15_lhloose_L1EM13VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e17_lhloose"); {
        *cutflow << HFTname("trig_e17_lhloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e17_lhloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e17_lhmedium"); {
        *cutflow << HFTname("trig_e17_lhmedium");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e17_lhmedium;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhmedium_iloose_L1EM18VH"); {
        *cutflow << HFTname("trig_e24_lhmedium_iloose_L1EM18VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhmedium_iloose_L1EM18VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhmedium_iloose_L1EM20VH"); {
        *cutflow << HFTname("trig_e24_lhmedium_iloose_L1EM20VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhmedium_iloose_L1EM20VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhmedium_L1EM20VHI"); {
        *cutflow << HFTname("trig_e24_lhmedium_L1EM20VHI");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhmedium_L1EM20VHI;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhtight_iloose"); {
        *cutflow << HFTname("trig_e24_lhtight_iloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhtight_iloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhtight_nod0_iloose"); {
        *cutflow << HFTname("trig_e24_lhtight_nod0_iloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhtight_nod0_iloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhmedium_nod0_iloose_L1EM20VH"); {
        *cutflow << HFTname("trig_e24_lhmedium_nod0_iloose_L1EM20VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhmedium_nod0_iloose_L1EM20VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhtight_ivarloose"); {
        *cutflow << HFTname("trig_e24_lhtight_ivarloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhtight_ivarloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhtight_nod0_ivarloose"); {
        *cutflow << HFTname("trig_e24_lhtight_nod0_ivarloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhtight_nod0_ivarloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhmedium_L1EM22VHI"); {
        *cutflow << HFTname("trig_e26_lhmedium_L1EM22VHI");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhmedium_L1EM22VHI;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhtight_iloose"); {
        *cutflow << HFTname("trig_e26_lhtight_iloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhtight_iloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhtight_nod0_iloose"); {
        *cutflow << HFTname("trig_e26_lhtight_iloose_nod0");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhtight_nod0_iloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhtight_ivarloose"); {
        *cutflow << HFTname("trig_e26_lhtight_ivarloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhtight_ivarloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhtight_nod0_ivarloose"); {
        *cutflow << HFTname("trig_e26_lhtight_nod0_ivarloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhtight_nod0_ivarloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e60_lhmedium"); {
        *cutflow << HFTname("trig_e60_lhmedium");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e60_lhmedium;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e60_lhmedium_nod0"); {
        *cutflow << HFTname("trig_e60_lhmedium_nod0");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e60_lhmedium_nod0;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e12_lhloose_L12EM10VH"); {
        *cutflow << HFTname("trig_2e12_lhloose_L12EM10VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e12_lhloose_L12EM10VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e15_lhloose_L12EM13VH"); {
        *cutflow << HFTname("trig_2e15_lhloose_L12EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e15_lhloose_L12EM13VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e15_lhvloose_L12EM13VH"); {
        *cutflow << HFTname("trig_2e15_lhvloose_L12EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e15_lhvloose_L12EM13VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e15_lhvloose_nod0_L12EM13VH"); {
        *cutflow << HFTname("trig_2e15_lhvloose_nod0_L12EM13VH");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e15_lhvloose_nod0_L12EM13VH;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e17_lhloose"); {
        *cutflow << HFTname("trig_2e17_lhloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e17_lhloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e17_lhvloose"); {
        *cutflow << HFTname("trig_2e17_lhvloose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e17_lhvloose;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_2e17_lhvloose_nod0"); {
        *cutflow << HFTname("trig_2e17_lhvloose_nod0");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_2e17_lhvloose_nod0;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e17_lhloose_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e17_lhloose_mu14;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e17_lhloose_nod0_mu14"); {
        *cutflow << HFTname("trig_e17_lhloose_nod0_mu14");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e17_lhloose_nod0_mu14;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e7_lhmedium_mu24"); {
        *cutflow << HFTname("trig_e7_lhmedium_mu24");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e7_lhmedium_mu24;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e7_lhmedium_nod0_mu24"); {
        *cutflow << HFTname("trig_e7_lhmedium_nod0_mu24");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e7_lhmedium_nod0_mu24;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhmedium_L1EM22VHI_mu8noL1"); {
        *cutflow << HFTname("trig_e26_lhmedium_L1EM22VHI_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhmedium_L1EM22VHI_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e26_lhmedium_nod0_L1EM22VHI_mu8noL1"); {
        *cutflow << HFTname("trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("pass_e24_lhmedium_L1EM20VHI_mu8noL1"); {
        *cutflow << HFTname("trig_e24_lhmedium_L1EM20VHI_mu8noL1");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return p_e24_lhmedium_L1EM20VHI_mu8noL1;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("mcid"); {
        *cutflow << HFTname("mcid");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->mcChannel;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("year (either 2015 or 2016 for data)"); {
        *cutflow << HFTname("year");
        *cutflow << [&](Superlink* sl, var_int*) -> int { 
            return sl->nt->evt()->treatAsYear;
        };
        *cutflow << SaveVar();
    }

    // standard variables
    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("event weight with Sherpa V+Jets weight"); {
        *cutflow << HFTname("eventweightVJets");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            double weight = sl->weights->product() * sl->nt->evt()->wPileup;
            if(sl->nt->evt()->isSherpaVjetsSample)
                weight *= sl->nt->evt()->sherpa22VjetsWeight;
            return weight;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("event weight (no pileup)"); {
        *cutflow << HFTname("eventweightNOPUPW");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product();
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("Pile-up weight"); {
        *cutflow << HFTname("pupw");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pile-up weight (up variation)"); {
        *cutflow << HFTname("pupw_up");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup_up;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Pile-up weight (down variation"); {
        *cutflow << HFTname("pupw_down");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup_dn;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Left polarized stop"); {
        *cutflow << HFTname("susy3BodyLeftPol");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyLeftPol;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Right polarized (70%) stop"); {
        *cutflow << HFTname("susy3BodyRightPol");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyRightPol;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Three body re-weighting -- Matrix element only"); {
        *cutflow << HFTname("susy3BodyOnlyMass");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->susy3BodyOnlyMass;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is MC"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of primary vertices"); {
        *cutflow << HFTname("nVtx");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->nVtx; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("average interactions per b.c."); {
        *cutflow << HFTname("avgMu");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->avgMu; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("average interactions per b.c. with data scale factor applied"); {
        *cutflow << HFTname("avgMuDataSF");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->avgMuDataSF; };
        *cutflow << SaveVar();
    } 

    // lepton variables
    // lepton variables
    // lepton variables

    LeptonVector leptons;
    ElectronVector electrons;
    MuonVector muons;
    *cutflow << [&](Superlink* sl, var_void*) { leptons = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { electrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { muons = *sl->muons; };

    *cutflow << NewVar("number of leptons"); {
        *cutflow << HFTname("nLeptons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return leptons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of electrons"); {
        *cutflow << HFTname("nElectrons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return electrons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of muons"); {
        *cutflow << HFTname("nMuons");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return muons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton charge"); {
        *cutflow << HFTname("l_q");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->q);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0"); {
        *cutflow << HFTname("l_d0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->d0);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0sig"); {
        *cutflow << HFTname("l_d0sig");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->d0sigBSCorr);
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0sinTheta"); {
        *cutflow << HFTname("l_z0sinTheta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->z0SinTheta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton pt"); {
        *cutflow << HFTname("l_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Pt());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("etconetopo20"); {
        *cutflow << HFTname("l_etconetopo20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("etconetopo30"); {
        *cutflow << HFTname("l_etconetopo30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->etconetopo30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptcone20"); {
        *cutflow << HFTname("l_ptcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptcone20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptcone30"); {
        *cutflow << HFTname("l_ptcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptcone30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptvarcone20"); {
        *cutflow << HFTname("l_ptvarcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("ptvarcone30"); {
        *cutflow << HFTname("l_ptvarcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->ptvarcone30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton eta"); {
        *cutflow << HFTname("l_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Eta());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton phi"); {
        *cutflow << HFTname("l_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < leptons.size(); i++) {
                out.push_back(leptons.at(i)->Phi());
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mll leptons"); {
        *cutflow << HFTname("mll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mll = -999.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                mll = (*l0 + *l1).M();
            }
            return mll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton pT"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double pTll = -999.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                pTll = (*l0 + *l1).Pt();
            }
            return pTll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between to leptons"); {
        *cutflow << HFTname("dphi_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double dphi = -999.0;
            if(leptons.size() == 2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                dphi = l0.DeltaPhi(l1);
            }
            return dphi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta eta between two leptons"); {
        *cutflow << HFTname("deta_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double deta = -999.0;
            if(leptons.size() == 2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                deta = l0.Eta() - l1.Eta();
            }
            return deta;
        };
        *cutflow << SaveVar();
    }

    // jet variables
    // jet variables
    // jet variables

    JetVector jets;
    JetVector bjets;
    JetVector sjets;

    JetVector bjets30;
    JetVector sjets30;
    JetVector bjets40;
    JetVector sjets40;
    JetVector bjets50;
    JetVector sjets50;
    JetVector bjets60;
    JetVector sjets60;

    *cutflow << [&](Superlink* sl, var_void*) { jets = *sl->jets; };
    *cutflow << [&](Superlink* sl, var_void*) {
        for(int i = 0; i < jets.size(); i++) {
            Jet* j = jets[i];
            if(sl->tools->jetSelector().isB(j))  bjets.push_back(j);
            if(!sl->tools->jetSelector().isB(j)) sjets.push_back(j);

            // 30
            if(sl->tools->jetSelector().isBMod(j, 85, 30)) bjets30.push_back(j);
            if(!sl->tools->jetSelector().isBMod(j, 85, 30)) sjets30.push_back(j);
            // 40
            if(sl->tools->jetSelector().isBMod(j, 85, 40)) bjets40.push_back(j);
            if(!sl->tools->jetSelector().isBMod(j, 85, 40)) sjets40.push_back(j);
            // 50
            if(sl->tools->jetSelector().isBMod(j, 85, 50)) bjets50.push_back(j);
            if(!sl->tools->jetSelector().isBMod(j, 85, 50)) sjets50.push_back(j);
            // 60
            if(sl->tools->jetSelector().isBMod(j, 85, 60)) bjets60.push_back(j);
            if(!sl->tools->jetSelector().isBMod(j, 85, 60)) sjets60.push_back(j);

        }// i
    };


    *cutflow << NewVar("number of jets"); {
        *cutflow << HFTname("nJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return jets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets"); {
        *cutflow << HFTname("nSJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets (b>30)"); {
        *cutflow << HFTname("nSJets30");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets30.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets (b>40)"); {
        *cutflow << HFTname("nSJets40");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets40.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets (b>50)"); {
        *cutflow << HFTname("nSJets50");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets50.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets (b>60)"); {
        *cutflow << HFTname("nSJets60");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets60.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of bjets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of bjets pt>30"); {
        *cutflow << HFTname("nBJets30");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets30.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of bjets pt>40"); {
        *cutflow << HFTname("nBJets40");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets40.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of bjets pt>50"); {
        *cutflow << HFTname("nBJets50");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets50.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of bjets pt>60"); {
        *cutflow << HFTname("nBJets60");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets60.size();
        };
        *cutflow << SaveVar();
    }




    // add cut
    //*cutflow << CutName("at least 2 b-tagged jets") << [&](Superlink* sl) -> bool {
    //    return (bjets.size()>=2);
    //};

    *cutflow << NewVar("jet pt"); {
        *cutflow << HFTname("j_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet pt"); {
        *cutflow << HFTname("sj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt"); {
        *cutflow << HFTname("bj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt 30"); {
        *cutflow << HFTname("bj_pt30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets30.size(); i++) {
                out.push_back(bjets30.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt 40"); {
        *cutflow << HFTname("bj_pt40");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets40.size(); i++) {
                out.push_back(bjets40.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt 50"); {
        *cutflow << HFTname("bj_pt50");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets50.size(); i++) {
                out.push_back(bjets50.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt 60"); {
        *cutflow << HFTname("bj_pt60");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets60.size(); i++) {
                out.push_back(bjets60.at(i)->Pt());
            }
            return out;
            };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("jet eta"); {
        *cutflow << HFTname("j_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet eta"); {
        *cutflow << HFTname("sj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet eta"); {
        *cutflow << HFTname("bj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Eta());
            }
            return out;
            };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet phi"); {
        *cutflow << HFTname("j_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < jets.size(); i++) {
                out.push_back(jets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet phi"); {
        *cutflow << HFTname("sj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet phi"); {
        *cutflow << HFTname("bj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Phi());
            }
            return out;
            };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between dilepton system and leading jet"); {
        *cutflow << HFTname("dphi_j0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(jets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = jets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading sjet"); {
        *cutflow << HFTname("dphi_j0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(jets.size()>0) {
                out = jets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading sjet"); {
        *cutflow << HFTname("dphi_sj0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(sjets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = sjets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading sjet"); {
        *cutflow << HFTname("dphi_sj0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(sjets.size()>0) {
                out = sjets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(bjets.size()>0) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = bjets.at(0)->DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_l0");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double out = -999;
            if(bjets.size()>0) {
                out = bjets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    // met variables
    // met variables
    // met variables
    Met met;
    *cutflow << [&](Superlink* sl, var_void*) { met = *sl->met; };
    *cutflow << NewVar("transverse missing energy (Etmiss)"); {
        *cutflow << HFTname("met");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Pt(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("phi coord. of Etmiss"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Phi(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met TST"); {
        *cutflow << HFTname("metTST");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.softTerm_et; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and met"); {
        *cutflow << HFTname("dphi_met_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0, l1, ll;
            l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
            ll = l0 + l1;
            return met.lv().DeltaPhi(ll);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mt2"); {
        *cutflow << HFTname("mt2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mt2 = -999.0;
            if(leptons.size() == 2) {
                mt2 = kin::getMT2(*sl->leptons, *sl->met);
            }
            return mt2;
        };
        *cutflow << SaveVar();
    }

    double meff;
    *cutflow << NewVar("meff : scalar sum pt of all jets, leptons, and met"); {
        *cutflow << HFTname("meff");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            meff = 0.0;
            // met
            meff += met.lv().Pt();
            // jets
            for(unsigned int ij = 0; ij < jets.size(); ij++){
                meff += jets.at(ij)->Pt();
            }
            // leptons
            for(unsigned int il=0; il < leptons.size(); il++){
                meff += leptons.at(il)->Pt();
            }
            return meff;
        };
        *cutflow << SaveVar();
    }
    double meff_S2L;
    *cutflow << NewVar("meff S2L : scalar sum pt of leptons, met, and up to two jets"); {
        *cutflow << HFTname("meff_S2L");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            meff_S2L = 0.0;
            // met
            meff_S2L += met.lv().Pt();
            // leptons
            for(int il=0; il < (int)leptons.size(); il++){
                meff_S2L += leptons.at(il)->Pt();
            }
            // jets
            int n_j = 0;
            for(int ij = 0; ij < (int)jets.size(); ij++){
                if(n_j < 2) {
                    meff_S2L += jets.at(ij)->Pt();
                    n_j++;
                }
            }
            return meff_S2L;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R1 : met / meff"); {
        *cutflow << HFTname("R1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R1 = -999.0;
            if(meff>0.0) {
                R1 = met.lv().Pt() / meff * 1.0;
            }
            return R1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R1 S2L : met / meff_S2L"); {
        *cutflow << HFTname("R1_S2L");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R1_S2L = -999.0;
            if(meff_S2L>0.0) {
                R1_S2L = met.lv().Pt() / (meff_S2L * 1.0);
            }
            return R1_S2L;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("R2 : met / (met + l0pt + l1pt)"); {
        *cutflow << HFTname("R2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double R2 = -999.0;
            if(leptons.size() == 2) {
                double denom = met.lv().Pt() + leptons.at(0)->Pt() + leptons.at(1)->Pt();
                R2 = met.lv().Pt() / denom * 1.0;
            }
            return R2;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("cosThetaB (WW-like)"); {
        *cutflow << HFTname("cosThetaB");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double cosThetaB = -999;
            if(leptons.size()==2) {
                TLorentzVector lp, lm;
                for(int il = 0; il < (int)leptons.size(); il++) {
                    if(leptons.at(il)->q < 0) lm = *leptons.at(il);
                    else if(leptons.at(il)->q > 0) lp = *leptons.at(il);
                } // il
                TLorentzVector ll = lp+lm;
                TVector3 boost = ll.BoostVector();
                lp.Boost(-boost);
                lm.Boost(-boost);
                cosThetaB = tanh((lp.Eta()-lm.Eta())/2.);
            }
            return cosThetaB;
        };
        *cutflow << SaveVar();
    }


    ///////////////////////////////////////////////
    // WWBB variables
    ///////////////////////////////////////////////

    // angles
    // dRll
    *cutflow << NewVar("delta R between two leptons"); {
        *cutflow << HFTname("dRll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return (leptons.at(0)->DeltaR(*leptons.at(1)));
        };
        *cutflow << SaveVar();
    }

    // M_bb
    *cutflow << NewVar("invariant mass of di-bjet system"); {
        *cutflow << HFTname("mbb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                return (*bjets.at(0) + *bjets.at(1)).M();
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dibjet mass vector"); {
        *cutflow << HFTname("mbb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            if(bjets.size()>=2) {
                out.push_back( (*bjets.at(0) + *bjets.at(1)).M());
            }
            else {
                out.push_back(-10.);
            }
            if(bjets30.size()>=2) {
                out.push_back( (*bjets30.at(0) + *bjets30.at(1)).M());
            }
            else {
                out.push_back(-10.);
            }
            if(bjets40.size()>=2) {
                out.push_back( (*bjets40.at(0) + *bjets40.at(1)).M());
            }
            else {
                out.push_back(-10.);
            }
            if(bjets50.size()>=2) {
                out.push_back( (*bjets50.at(0) + *bjets50.at(1)).M());
            }
            else {
                out.push_back(-10.);
            }
            if(bjets60.size()>=2) {
                out.push_back( (*bjets60.at(0) + *bjets60.at(1)).M());
            }
            else {
                out.push_back(-10.);
            }

            return out;
        };
        *cutflow << SaveVar();
    }

    // dRbb
    *cutflow << NewVar("delta R between two leading b-jets"); {
        *cutflow << HFTname("dRbb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                return (bjets.at(0)->DeltaR(*bjets.at(1)));
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta R between two leading b-jet vec"); {
        *cutflow << HFTname("dRbb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            if(bjets.size()>=2) {
                out.push_back( bjets.at(0)->DeltaR(*bjets.at(1)) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets30.size()>=2) {
                out.push_back( bjets30.at(0)->DeltaR(*bjets30.at(1)) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets40.size()>=2) {
                out.push_back( bjets40.at(0)->DeltaR(*bjets40.at(1)) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets50.size()>=2) {
                out.push_back( bjets50.at(0)->DeltaR(*bjets50.at(1)) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets60.size()>=2) {
                out.push_back( bjets60.at(0)->DeltaR(*bjets60.at(1)) );
            }
            else {
                out.push_back(-10.0);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    // dR_ll_bb
    *cutflow << NewVar("delta R between dilepton system and di-bjet system"); {
        *cutflow << HFTname("dR_ll_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));

                return ( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dR_ll_bb_vec"); {
        *cutflow << HFTname("dR_ll_bb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            TLorentzVector l0 = (*leptons.at(0));
            TLorentzVector l1 = (*leptons.at(1));
            if(bjets.size()>=2) {
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));

                out.push_back( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets30.size()>=2) {
                TLorentzVector b0 = (*bjets30.at(0));
                TLorentzVector b1 = (*bjets30.at(1));

                out.push_back( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets40.size()>=2) {
                TLorentzVector b0 = (*bjets40.at(0));
                TLorentzVector b1 = (*bjets40.at(1));

                out.push_back( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets50.size()>=2) {
                TLorentzVector b0 = (*bjets50.at(0));
                TLorentzVector b1 = (*bjets50.at(1));

                out.push_back( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets60.size()>=2) {
                TLorentzVector b0 = (*bjets60.at(0));
                TLorentzVector b1 = (*bjets60.at(1));

                out.push_back( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            else {
                out.push_back(-10.0);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    // dphi bb  ll
    *cutflow << NewVar("delta phi between bb and ll systems"); {
        *cutflow << HFTname("dphi_ll_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                return ( (*bjets.at(0) + *bjets.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1)) ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dphi_bb_ll_vec"); {
        *cutflow << HFTname("dphi_ll_bb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            if(bjets.size()>=2) {
                out.push_back((*bjets.at(0) + *bjets.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1))));
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets30.size()>=2) {
                out.push_back((*bjets30.at(0) + *bjets30.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1))));
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets40.size()>=2) {
                out.push_back((*bjets40.at(0) + *bjets40.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1))));
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets50.size()>=2) {
                out.push_back((*bjets50.at(0) + *bjets50.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1))));
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets60.size()>=2) {
                out.push_back((*bjets60.at(0) + *bjets60.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1))));
            }
            else {
                out.push_back(-10.0);
            }
            return out;
        };
        *cutflow << SaveVar();
    }


    // dphi WW bb
    *cutflow << NewVar("delta phi between WW and bb systems"); {
        *cutflow << HFTname("dphi_WW_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                return ( (met.lv() + *leptons.at(0) + *leptons.at(1)).DeltaPhi( (*bjets.at(0) + *bjets.at(1)) ) ); 
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    float delta_phi_x = -10.;
    float delta_phi_x_scaled = -10.;
    float delta_phi_x_scaled_ww = -10.;

    // mass_X
    *cutflow << NewVar("mass of X"); {
        *cutflow << HFTname("mass_X");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                delta_phi_x = (met.lv() + *leptons.at(0) + * leptons.at(1)).DeltaPhi((*bjets.at(0) + *bjets.at(1)));
                return ( (met.lv() + *leptons.at(0) + *leptons.at(1) + *bjets.at(0) + *bjets.at(1)).M() );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // mass_X
    *cutflow << NewVar("mass of X trans"); {
        *cutflow << HFTname("mass_X_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));
                l0.SetPz(0.0);
                l1.SetPz(0.0);
                b0.SetPz(0.0);
                b1.SetPz(0.0);
                return ( ( met.lv() + l0 + l1 + b0 + b1).M() );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of X trans"); {
        *cutflow << HFTname("mass_X_T_2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));
                return ( ( met.lv() + l0 + l1 + b0 + b1).Mt() );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }


    // mass_X scaled
    *cutflow << NewVar("mass of X with bb scaling"); {
        *cutflow << HFTname("mass_X_scaled");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                float m_bb = (*bjets.at(0) + *bjets.at(1)).M();
                float scaling = 125.09 / m_bb;
                TLorentzVector bjet_system = (*bjets.at(0) + *bjets.at(1));
                bjet_system.SetPtEtaPhiE( bjet_system.Pt() * scaling, bjet_system.Eta(),
                                            bjet_system.Phi(), bjet_system.E() * scaling );

                delta_phi_x_scaled = (met.lv() + *leptons.at(0) + *leptons.at(1)).DeltaPhi(bjet_system);

                return ( (met.lv() + *leptons.at(0) + *leptons.at(1) + bjet_system).M() );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // dphi_met_ll
    *cutflow << NewVar("delta phi between MET and dilepton system"); {
        *cutflow << HFTname("dphi_met_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return ( (met.lv().DeltaPhi( (*leptons.at(0) + *leptons.at(1)) )) );
        };
        *cutflow << SaveVar();
    }

    // mass_met_ll
    *cutflow << NewVar("mass of met and dilepton system"); {
        *cutflow << HFTname("mass_met_ll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return ( (met.lv() + *leptons.at(0) + *leptons.at(1)).M() );
        };
        *cutflow << SaveVar();
    }

    // mass_met_ll_T
    *cutflow << NewVar("mass of met and dilepton system transv"); {
        *cutflow << HFTname("mass_met_ll_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0 = (*leptons.at(0));
            TLorentzVector l1 = (*leptons.at(1));
            l0.SetPz(0.0);
            l1.SetPz(0.0);
            return ( ( met.lv() + l0 + l1).M() );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of met and dilepton system transv"); {
        *cutflow << HFTname("mass_met_ll_T_2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0 = (*leptons.at(0));
            TLorentzVector l1 = (*leptons.at(1));
            return ( ( met.lv() + l0 + l1).Mt() );
        };
        *cutflow << SaveVar();
    }

    // met_pTll
    *cutflow << NewVar("pT of met + dilepton ssytem"); {
        *cutflow << HFTname("met_pTll");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return ( (met.lv() + *leptons.at(0) + *leptons.at(1)).Pt() );
        };
        *cutflow << SaveVar();
    }

    // HT2
    *cutflow << NewVar("HT2"); {
        *cutflow << HFTname("HT2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                double HT2 = ( (*bjets.at(0) + *bjets.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );
                return HT2;
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // HT2Ratio
    *cutflow << NewVar("HT2Ratio"); {
        *cutflow << HFTname("HT2Ratio");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {

                double num = ( (*bjets.at(0) + *bjets.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets.at(0)).Pt());
                den += (*bjets.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();

                return (num / den);
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("HT2Ratio vec"); {
        *cutflow << HFTname("HT2Ratio_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;

            //////// pt 20
            if(bjets.size() >= 2) {
                double num = ( (*bjets.at(0) + *bjets.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets.at(0)).Pt());
                den += (*bjets.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                double ht2 = (num / den);
                out.push_back(ht2);
            }
            else {
                out.push_back(-10.0);
            }
            //////// pt 30
            if(bjets30.size() >= 2) {
                double num = ( (*bjets30.at(0) + *bjets30.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets30.at(0)).Pt());
                den += (*bjets30.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                double ht2 = (num / den);
                out.push_back(ht2);
            }
            else {
                out.push_back(-10.0);
            }
            //////// pt 40
            if(bjets40.size() >= 2) {
                double num = ( (*bjets40.at(0) + *bjets40.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets40.at(0)).Pt());
                den += (*bjets40.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                double ht2 = (num / den);
                out.push_back(ht2);
            }
            else {
                out.push_back(-10.0);
            }
            //////// pt 50
            if(bjets50.size() >= 2) {
                double num = ( (*bjets50.at(0) + *bjets50.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets50.at(0)).Pt());
                den += (*bjets50.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                double ht2 = (num / den);
                out.push_back(ht2);
            }
            else {
                out.push_back(-10.0);
            }
            //////// pt 60
            if(bjets60.size() >= 2) {
                double num = ( (*bjets60.at(0) + *bjets60.at(1)).Pt() + 
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((*bjets60.at(0)).Pt());
                den += (*bjets60.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                double ht2 = (num / den);
                out.push_back(ht2);
            }
            else {
                out.push_back(-10.0);
            }

            return out;

        };
        *cutflow << SaveVar();
    }

    // MT_HWW
    *cutflow << NewVar("MT_HWW"); {
        *cutflow << HFTname("MT_HWW");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double ptll_met = (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt(); 
            double ptll2 = (*leptons.at(0) + *leptons.at(1)).Pt();
            ptll2 = ptll2 * ptll2;
            double mll2 = (*leptons.at(0) + *leptons.at(1)).M();
            mll2 = mll2 * mll2;
            double ET_ll = sqrt(ptll2 + mll2);
            return sqrt( (ET_ll + met.lv().Pt())*(ET_ll + met.lv().Pt()) - (ptll_met * ptll_met) );
        };
        *cutflow << SaveVar();
    }

    // MT_1
    *cutflow << NewVar("MT_1"); {
        *cutflow << HFTname("MT_1");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1)
                            + *bjets.at(0) + *bjets.at(1));
                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);
                return sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ( (vis + met.lv()).Pt() * (vis + met.lv()).Pt() ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // MT_1_scaled
    *cutflow << NewVar("MT_1_scaled"); {
        *cutflow << HFTname("MT_1_scaled");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector bjet_system = (*bjets.at(0) + *bjets.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                return ( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) );
            }
            return -10.;

        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("MT_1_test"); {
        *cutflow << HFTname("MT_1_test");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector bjet_system = (*bjets.at(0) + * bjets.at(1));
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);
                return (vis + met.lv()).Mt();
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("MT_1_test_scaled"); {
        *cutflow << HFTname("MT_1_test_scaled");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector bjet_system = (*bjets.at(0) + * bjets.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);
                return (vis + met.lv()).Mt();
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
        
        

    // MT_1_scaled with h->WW scaled in addition to h->bb
    *cutflow << NewVar("mass_X_scaled_ww"); {
        *cutflow << HFTname("mass_X_scaled_ww");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {

                TLorentzVector bjet_system = (*bjets.at(0) + *bjets.at(1));
                double mbb = bjet_system.M();
                double scaling = 125.09/mbb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
    
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets.at(0) + *bjets.at(1) );
                float mt2_llbb = kin::getMT2( v0, v1, met );

                TLorentzVector ww_system = (*leptons.at(0) + *leptons.at(1) + met.lv());
                //double mww = ww_system.M();
                double ww_scaling = 125.09/mt2_llbb;
                ww_system.SetPtEtaPhiE(ww_system.Pt() * ww_scaling, ww_system.Eta(), ww_system.Phi(), ww_system.E() * ww_scaling);

                delta_phi_x_scaled_ww = ww_system.DeltaPhi(bjet_system);

                return (bjet_system + ww_system).M();
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta_phi_x"); {
        *cutflow << HFTname("delta_phi_x");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return delta_phi_x;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta_phi_x_scaled"); {
        *cutflow << HFTname("delta_phi_x_scaled");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return delta_phi_x_scaled;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta_phi_x_scaled_ww"); {
        *cutflow << HFTname("delta_phi_x_scaled_ww");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return delta_phi_x_scaled_ww;
        };
        *cutflow << SaveVar();
    }



    *cutflow << NewVar("MT1_scaled_vec"); {
        *cutflow << HFTname("MT_1_scaled_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;

            ////////////////pt 20
            if(bjets.size()>=2) {
                TLorentzVector bjet_system = (*bjets.at(0) + *bjets.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                out.push_back(( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) ));
            }
            else {
                out.push_back(-10.);
            }
            ////////////////pt 30
            if(bjets30.size()>=2) {
                TLorentzVector bjet_system = (*bjets30.at(0) + *bjets30.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                out.push_back(( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) ));
            }
            else {
                out.push_back(-10.);
            }
            ////////////////pt 40
            if(bjets40.size()>=2) {
                TLorentzVector bjet_system = (*bjets40.at(0) + *bjets40.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                out.push_back(( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) ));
            }
            else {
                out.push_back(-10.);
            }
            ////////////////pt 50
            if(bjets50.size()>=2) {
                TLorentzVector bjet_system = (*bjets50.at(0) + *bjets50.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                out.push_back(( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) ));
            }
            else {
                out.push_back(-10.);
            }
            ////////////////pt 60
            if(bjets60.size()>=2) {
                TLorentzVector bjet_system = (*bjets60.at(0) + *bjets60.at(1));
                double m_bb = bjet_system.M();
                double scaling = 125.09/m_bb;
                bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(), bjet_system.Phi(), bjet_system.E() * scaling);
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1) + bjet_system);

                double pt_vis = vis.Pt();
                double m_vis = vis.M();
                double et_vis = sqrt(pt_vis * pt_vis + m_vis * m_vis);

                out.push_back(( sqrt( (et_vis + met.lv().Pt()) * (et_vis + met.lv().Pt()) -
                            ((vis + met.lv()).Pt() * (vis + met.lv()).Pt()) ) ));
            }
            else {
                out.push_back(-10.);
            }

            return out;
        };
        *cutflow << SaveVar();
    }

    // mt1 (l0, b0), (l1, b1)
    *cutflow << NewVar("mt2_00"); {
        *cutflow << HFTname("mt2_00");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *bjets.at(0) );
                const TLorentzVector v1 = ( *leptons.at(1) + *bjets.at(1) );
                return kin::getMT2( v0, v1, met );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // mt1 (l0, b1), (l1, b0)
    *cutflow << NewVar("mt2_01"); {
        *cutflow << HFTname("mt2_01");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *bjets.at(1) );
                const TLorentzVector v1 = ( *leptons.at(1) + *bjets.at(0) );
                return kin::getMT2( v0, v1, met );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // mt1 (l1, b0), (l0, b1)
    *cutflow << NewVar("mt2_10"); {
        *cutflow << HFTname("mt2_10");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(1) + *bjets.at(0) );
                const TLorentzVector v1 = ( *leptons.at(0) + *bjets.at(1) );
                return kin::getMT2( v0, v1, met );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // mt1 (l1, b1), (l0, b0)
    *cutflow << NewVar("mt2_11"); {
        *cutflow << HFTname("mt2_11");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(1) + *bjets.at(1) );
                const TLorentzVector v1 = ( *leptons.at(0) + *bjets.at(0) );
                return kin::getMT2( v0, v1, met );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // mt1 (l0, l1), (b0, b1)
    *cutflow << NewVar("mt2_llbb"); {
        *cutflow << HFTname("mt2_llbb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets.at(0) + *bjets.at(1) );
                return kin::getMT2( v0, v1, met );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mt2_llbb_vec"); {
        *cutflow << HFTname("mt2_llbb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;

            ///////// pt 20
            if(bjets.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets.at(0) + *bjets.at(1) );
                out.push_back(kin::getMT2(v0, v1, met));
            } 
            else {
                out.push_back(-10.0);
            }
            ///////// pt 30
            if(bjets30.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets30.at(0) + *bjets30.at(1) );
                out.push_back(kin::getMT2(v0, v1, met));
            } 
            else {
                out.push_back(-10.0);
            }
            ///////// pt 40
            if(bjets40.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets40.at(0) + *bjets40.at(1) );
                out.push_back(kin::getMT2(v0, v1, met));
            } 
            else {
                out.push_back(-10.0);
            }
            ///////// pt 50
            if(bjets50.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets50.at(0) + *bjets50.at(1) );
                out.push_back(kin::getMT2(v0, v1, met));
            } 
            else {
                out.push_back(-10.0);
            }
            ///////// pt 60
            if(bjets60.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( *bjets60.at(0) + *bjets60.at(1) );
                out.push_back(kin::getMT2(v0, v1, met));
            } 
            else {
                out.push_back(-10.0);
            }

            return out;
        };
        *cutflow << SaveVar();
    }


    // mt2_bb
    *cutflow << NewVar("mt2_bb"); {
        *cutflow << HFTname("mt2_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector b0 = (*bjets.at(0));
                const TLorentzVector b1 = (*bjets.at(1));
                return kin::getMT2(b0,b1,met);
            }
            return -10.;

        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mt2_bb_vec"); {
        *cutflow << HFTname("mt2_bb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;

            ////////////pt 20
            if(bjets.size()>=2) {
                const TLorentzVector b0 = (*bjets.at(0));
                const TLorentzVector b1 = (*bjets.at(1));
                out.push_back(kin::getMT2(b0,b1,met));
            }
            else {
                out.push_back(-10.0);
            }
            ////////////pt 30
            if(bjets30.size()>=2) {
                const TLorentzVector b0 = (*bjets30.at(0));
                const TLorentzVector b1 = (*bjets30.at(1));
                out.push_back(kin::getMT2(b0,b1,met));
            }
            else {
                out.push_back(-10.0);
            }
            ////////////pt 40
            if(bjets40.size()>=2) {
                const TLorentzVector b0 = (*bjets40.at(0));
                const TLorentzVector b1 = (*bjets40.at(1));
                out.push_back(kin::getMT2(b0,b1,met));
            }
            else {
                out.push_back(-10.0);
            }
            ////////////pt 50
            if(bjets50.size()>=2) {
                const TLorentzVector b0 = (*bjets50.at(0));
                const TLorentzVector b1 = (*bjets50.at(1));
                out.push_back(kin::getMT2(b0,b1,met));
            }
            else {
                out.push_back(-10.0);
            }
            ////////////pt 60
            if(bjets60.size()>=2) {
                const TLorentzVector b0 = (*bjets60.at(0));
                const TLorentzVector b1 = (*bjets60.at(1));
                out.push_back(kin::getMT2(b0,b1,met));
            }
            else {
                out.push_back(-10.0);
            }

            return out;
        };
        *cutflow << SaveVar();
    }


    // mt2_lvis
    *cutflow << NewVar("mt2_lvis"); {
        *cutflow << HFTname("mt2_lvis");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));

                b0.SetPz(0.0);
                b1.SetPz(0.0);

                Susy::Met met_t = met;
                double phi = (met.lv() + b0 + b1).Phi();
                double et = (met.lv() + b0 + b1).E();

                met_t.Et = et;
                met_t.phi = phi;

                return kin::getMT2(l0, l1, met_t);
            };
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // mt2_bvis
    *cutflow << NewVar("mt2_bvis"); {
        *cutflow << HFTname("mt2_bvis");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                const TLorentzVector b0 = (*bjets.at(0));
                const TLorentzVector b1 = (*bjets.at(1));

                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                l0.SetPz(0.0);
                l1.SetPz(0.0);
                const TLorentzVector l0t = l0;
                const TLorentzVector l1t = l1;

                Susy::Met met_t = met;
                double phi = (met.lv() + l0 + l1).Phi();
                double et = (met.lv() + l0 + l1).E();

                met_t.Et = et;
                met_t.phi = phi;

                return kin::getMT2(b0, b1, met_t);
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }


    

    // cosTheta2
    *cutflow << NewVar("cosTheta2"); {
        *cutflow << HFTname("cosTheta2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector llsystem = (*leptons.at(0) + *leptons.at(1));
            TVector3 boost = (llsystem + met.lv()).BoostVector();
            TLorentzVector lp;
            TLorentzVector lm;
            for(int i = 0; i < 2; i++) {
                if(leptons.at(i)->q < 0) lm = (*leptons.at(i));
                else if(leptons.at(i)->q > 0) lp = (*leptons.at(i));
            }

            lp.Boost(-boost);
            lm.Boost(-boost);
            return tanh((lp.Eta()-lm.Eta())/2.);
        };
        *cutflow << SaveVar();
    }

    // dphi_bb
    *cutflow << NewVar("dphi_bb"); {
        *cutflow << HFTname("dphi_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                return (bjets.at(0)->DeltaPhi(*bjets.at(1)));
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("dphi_bb_vec"); {
        *cutflow << HFTname("dphi_bb_vec");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            if(bjets.size()>=2) {
                out.push_back( (bjets.at(0)->DeltaPhi(*bjets.at(1))) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets30.size()>=2) {
                out.push_back( (bjets30.at(0)->DeltaPhi(*bjets30.at(1))) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets40.size()>=2) {
                out.push_back( (bjets40.at(0)->DeltaPhi(*bjets40.at(1))) );
            }
            else {
                out.push_back(-10.0);
            }
            if(bjets50.size()>=2) {
                out.push_back( (bjets50.at(0)->DeltaPhi(*bjets50.at(1))) );
            }
            else {
                out.push_back(-10.0);
            }

            if(bjets60.size()>=2) {
                out.push_back( (bjets60.at(0)->DeltaPhi(*bjets60.at(1))) );
            }
            else {
                out.push_back(-10.0);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    // dphi_boost_CM_ll_met
    *cutflow << NewVar("dphi_boost_CM_ll_met"); {
        *cutflow << HFTname("dhpi_boost_CM_ll_met");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));
                TLorentzVector met_tlv = met.lv();

                TVector3 boost_to_cm = (l0 + l1 + met.lv() + b0 + b1).BoostVector();

                l0.Boost(-boost_to_cm);
                l1.Boost(-boost_to_cm);
                met_tlv.Boost(-boost_to_cm);

                TVector3 boost = (l0 + l1 + met_tlv).BoostVector();
                l0.Boost(-boost);
                l1.Boost(-boost);
                return (l0 + l1).Vect().DeltaPhi(boost);
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // dphi_boost_CM_ll_met_T
    *cutflow << NewVar("dphi_boost_CM_ll_met_T"); {
        *cutflow << HFTname("dhpi_boost_CM_ll_met_T");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));
                TLorentzVector met_tlv = met.lv();

                l0.SetPz(0.0);
                l1.SetPz(0.0);
                b0.SetPz(0.0);
                b1.SetPz(0.0);

                TVector3 boost_to_cm = (l0 + l1 + met.lv() + b0 + b1).BoostVector();

                l0.Boost(-boost_to_cm);
                l1.Boost(-boost_to_cm);
                met_tlv.Boost(-boost_to_cm);

                TVector3 boost = (l0 + l1 + met_tlv).BoostVector();
                l0.Boost(-boost);
                l1.Boost(-boost);
                return (l0 + l1).Vect().DeltaPhi(boost);
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // mT_llmet
    *cutflow << NewVar("mT_llmet"); {
        *cutflow << HFTname("mT_llmet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            TLorentzVector l0t = (*leptons.at(0));
            TLorentzVector l1t = (*leptons.at(1));
            l0t.SetPz(0.0);
            l1t.SetPz(0.0);
            return ( (l0t + l1t).M() );
        };
        *cutflow << SaveVar();
    }

    // mT_bb
    *cutflow << NewVar("mT_bb"); {
        *cutflow << HFTname("mT_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size() >= 2) {
                TLorentzVector b0t = (*bjets.at(0));
                TLorentzVector b1t = (*bjets.at(1));
                return (b0t + b1t).M();
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // min/max
    *cutflow << NewVar("max_mT_llmet_bb"); {
        *cutflow << HFTname("max_mT_llmet_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));

                l0.SetPz(0.0);
                l1.SetPz(0.0);
                b0.SetPz(0.0);
                b1.SetPz(0.0);

                double mb = (b0 + b1).M();
                double ml = (l0 + l1).M();

                double out = 0.0;
                if(mb > ml) out = mb;
                else { out = ml; }

                return out;
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }
    // min/max
    *cutflow << NewVar("min_mT_llmet_bb"); {
        *cutflow << HFTname("min_mT_llmet_bb");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(bjets.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (*bjets.at(0));
                TLorentzVector b1 = (*bjets.at(1));

                l0.SetPz(0.0);
                l1.SetPz(0.0);
                b0.SetPz(0.0);
                b1.SetPz(0.0);

                double mb = (b0 + b1).M();
                double ml = (l0 + l1).M();

                double out = 0.0;
                if(mb > ml) out = ml;
                else { out = mb; }

                return out;
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    vector<float> event_shapes;
    vector<float> event_shapes_b20;
    vector<float> event_shapes_b30;
    vector<float> event_shapes_b40;
    vector<float> event_shapes_b50;
    vector<float> event_shapes_b60;

    *cutflow << [&](Superlink* sl, var_void*) {

        vector<float> px;
        vector<float> py;
        vector<float> pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        for(auto x : jets) px.push_back(x->Px());
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        for(auto x : jets) py.push_back(x->Py());
        py.push_back(met.lv().Py());

        // px vector
        for(auto x : leptons) pz.push_back(x->Pz());
        for(auto x : jets) pz.push_back(x->Pz());

        event_shapes = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        vector<float> px;
        vector<float> py;
        vector<float>pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        if(bjets.size()>=2) {
            px.push_back(bjets.at(0)->Px());
            px.push_back(bjets.at(1)->Px());
        } 
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        if(bjets.size()>=2) {
            py.push_back(bjets.at(0)->Py());
            py.push_back(bjets.at(1)->Py());
        }
        py.push_back(met.lv().Py());

        // pz vector
        for(auto x : leptons) pz.push_back(x->Pz());
        if(bjets.size()>=2) {
            pz.push_back(bjets.at(0)->Pz());
            pz.push_back(bjets.at(1)->Pz());
        }

        event_shapes_b20 = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        vector<float> px;
        vector<float> py;
        vector<float>pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        if(bjets30.size()>=2) {
            px.push_back(bjets30.at(0)->Px());
            px.push_back(bjets30.at(1)->Px());
        } 
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        if(bjets30.size()>=2) {
            py.push_back(bjets30.at(0)->Py());
            py.push_back(bjets30.at(1)->Py());
        }
        py.push_back(met.lv().Py());

        // pz vector
        for(auto x : leptons) pz.push_back(x->Pz());
        if(bjets30.size()>=2) {
            pz.push_back(bjets30.at(0)->Pz());
            pz.push_back(bjets30.at(1)->Pz());
        }

        event_shapes_b30 = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        vector<float> px;
        vector<float> py;
        vector<float>pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        if(bjets40.size()>=2) {
            px.push_back(bjets40.at(0)->Px());
            px.push_back(bjets40.at(1)->Px());
        } 
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        if(bjets40.size()>=2) {
            py.push_back(bjets40.at(0)->Py());
            py.push_back(bjets40.at(1)->Py());
        }
        py.push_back(met.lv().Py());

        // pz vector
        for(auto x : leptons) pz.push_back(x->Pz());
        if(bjets40.size()>=2) {
            pz.push_back(bjets40.at(0)->Pz());
            pz.push_back(bjets40.at(1)->Pz());
        }

        event_shapes_b40 = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        vector<float> px;
        vector<float> py;
        vector<float>pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        if(bjets50.size()>=2) {
            px.push_back(bjets50.at(0)->Px());
            px.push_back(bjets50.at(1)->Px());
        } 
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        if(bjets50.size()>=2) {
            py.push_back(bjets50.at(0)->Py());
            py.push_back(bjets50.at(1)->Py());
        }
        py.push_back(met.lv().Py());

        // pz vector
        for(auto x : leptons) pz.push_back(x->Pz());
        if(bjets50.size()>=2) {
            pz.push_back(bjets50.at(0)->Pz());
            pz.push_back(bjets50.at(1)->Pz());
        }

        event_shapes_b50 = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << [&](Superlink* sl, var_void*) {
        vector<float> px;
        vector<float> py;
        vector<float>pz;

        // px vector
        for(auto x : leptons) px.push_back(x->Px());
        if(bjets60.size()>=2) {
            px.push_back(bjets60.at(0)->Px());
            px.push_back(bjets60.at(1)->Px());
        } 
        px.push_back(met.lv().Px());

        // py vector
        for(auto x : leptons) py.push_back(x->Py());
        if(bjets60.size()>=2) {
            py.push_back(bjets60.at(0)->Py());
            py.push_back(bjets60.at(1)->Py());
        }
        py.push_back(met.lv().Py());

        // pz vector
        for(auto x : leptons) pz.push_back(x->Pz());
        if(bjets60.size()>=2) {
            pz.push_back(bjets60.at(0)->Pz());
            pz.push_back(bjets60.at(1)->Pz());
        }

        event_shapes_b60 = PhysicsTools::event_shapes(px, py, pz);
    };

    *cutflow << NewVar("sphericity"); {
        *cutflow << HFTname("sphericity");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity"); {
        *cutflow << HFTname("transverse_sphericity");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity"); {
        *cutflow << HFTname("aplanarity");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity"); {
        *cutflow << HFTname("mod_aplanarity");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes.at(2) ); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sphericity_b20"); {
        *cutflow << HFTname("sphericity_b20");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b20.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity_b20"); {
        *cutflow << HFTname("transverse_sphericity_b20");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b20.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity_b20"); {
        *cutflow << HFTname("aplanarity_b20");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b20.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity_b20"); {
        *cutflow << HFTname("mod_aplanarity_b20");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes_b20.at(2) ); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sphericity_b30"); {
        *cutflow << HFTname("sphericity_b30");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b30.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity_b30"); {
        *cutflow << HFTname("transverse_sphericity_b30");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b30.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity_b30"); {
        *cutflow << HFTname("aplanarity_b30");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b30.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity_b30"); {
        *cutflow << HFTname("mod_aplanarity_b30");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes_b30.at(2) ); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sphericity_b40"); {
        *cutflow << HFTname("sphericity_b40");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b40.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity_b40"); {
        *cutflow << HFTname("transverse_sphericity_b40");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b40.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity_b40"); {
        *cutflow << HFTname("aplanarity_b40");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b40.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity_b40"); {
        *cutflow << HFTname("mod_aplanarity_b40");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes_b40.at(2) ); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sphericity_b50"); {
        *cutflow << HFTname("sphericity_b50");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b50.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity_b50"); {
        *cutflow << HFTname("transverse_sphericity_b50");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b50.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity_b50"); {
        *cutflow << HFTname("aplanarity_b50");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b50.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity_b50"); {
        *cutflow << HFTname("mod_aplanarity_b50");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes_b50.at(2) ); 
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sphericity_b60"); {
        *cutflow << HFTname("sphericity_b60");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b60.at(0);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("transverse_sphericity_b60"); {
        *cutflow << HFTname("transverse_sphericity_b60");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b60.at(1);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("aplanarity_b60"); {
        *cutflow << HFTname("aplanarity_b60");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return event_shapes_b60.at(2);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mod_aplanarity_b60"); {
        *cutflow << HFTname("mod_aplanarity_b60");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return exp( -8. * event_shapes_b60.at(2) ); 
        };
        *cutflow << SaveVar();
    }




    *cutflow << [&](Superlink* sl, var_void*) { event_shapes.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { event_shapes_b20.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { event_shapes_b30.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { event_shapes_b40.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { event_shapes_b50.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { event_shapes_b60.clear(); };

    *cutflow << [&](Superlink* sl, var_void*) { delta_phi_x = -10.; delta_phi_x_scaled = -10.; delta_phi_x_scaled_ww = -10.; };


    // clear the wectors
    *cutflow << [&](Superlink* sl, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets30.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets40.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets50.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets60.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets30.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets40.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets50.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets60.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Sysystematics [BEGIN]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////
    // weight systematics
    ////////////////////////////////////

    // electron eff
    *cutflow << NewSystematic("shift in electron ID efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_UP, SupersysWeight::EL_EFF_ID_DN);
        *cutflow << TreeName("EL_EFF_ID");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron ISO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ISO_UP, SupersysWeight::EL_EFF_ISO_DN);
        *cutflow << TreeName("EL_EFF_Iso");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in electron RECO efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_RECO_UP, SupersysWeight::EL_EFF_RECO_DN);
        *cutflow << TreeName("EL_EFF_Reco");
        *cutflow << SaveSystematic();
    }

    // muon eff
    *cutflow << NewSystematic("muon eff stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_UP, SupersysWeight::MUON_EFF_STAT_DN);
        *cutflow << TreeName("MUON_EFF_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff stat uncertainty (low pt)"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_LOWPT_UP, SupersysWeight::MUON_EFF_STAT_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_STAT_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_UP, SupersysWeight::MUON_EFF_SYS_DN);
        *cutflow << TreeName("MUON_EFF_SYS");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff syst uncertainty (low pt"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_LOWPT_UP, SupersysWeight::MUON_EFF_SYS_LOWPT_DN);
        *cutflow << TreeName("MUON_EFF_SYS_LOWPT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso stat uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_STAT_UP, SupersysWeight::MUON_EFF_ISO_STAT_DN);
        *cutflow << TreeName("MUON_ISO_STAT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon eff iso syst uncertainty"); {
        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_SYS_UP, SupersysWeight::MUON_EFF_ISO_SYS_DN);
        *cutflow << TreeName("MUON_ISO_SYS");
        *cutflow << SaveSystematic();
    }

    // flavor tagging eff
    *cutflow << NewSystematic("shift in b-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_B_UP, SupersysWeight::FT_EFF_B_DN);
        *cutflow << TreeName("FT_EFF_B");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in c-tag efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_C_UP, SupersysWeight::FT_EFF_C_DN);
        *cutflow << TreeName("FT_EFF_C");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in light tag (i.e. mis-tag) efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_LT_UP, SupersysWeight::FT_EFF_LT_DN);
        *cutflow << TreeName("FT_EFF_Light");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAP_UP, SupersysWeight::FT_EFF_EXTRAP_DN);
        *cutflow << TreeName("FT_EFF_extrapolation");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift flavor tagging extrapolation - charm ?"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAPC_UP, SupersysWeight::FT_EFF_EXTRAPC_DN);
        *cutflow << TreeName("FT_EFF_extrapolation_charm");
        *cutflow << SaveSystematic();
    }

    // jvt eff
    *cutflow << NewSystematic("shift in JVT efficiency"); {
        *cutflow << WeightSystematic(SupersysWeight::JVT_EFF_UP, SupersysWeight::JVT_EFF_DN);
        *cutflow << TreeName("JET_JVTEff");
        *cutflow << SaveSystematic();
    }

    // pileup
    *cutflow << NewSystematic("shift in data mu (pile-up)"); {
        *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DN);
        *cutflow << TreeName("PILEUP");
        *cutflow << SaveSystematic();
    }


    ////////////////////////////////////
    // shape systematics
    ////////////////////////////////////

    // egamma
    *cutflow << NewSystematic("shift in e-gamma resolution (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
        *cutflow << TreeName("EG_RESOLUTION_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma resolution (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
        *cutflow << TreeName("EG_RESOLUTION_ALL_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (UP)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
        *cutflow << TreeName("EG_SCALE_ALL_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("shift in e-gamma scale (DOWN)"); {
        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
        *cutflow << TreeName("EG_SCALE_ALL_DN");
        *cutflow << SaveSystematic();
    }
    // muon
    *cutflow << NewSystematic("muon ID (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_ID_UP);
        *cutflow << TreeName("MUON_ID_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon ID (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUON_ID_DN);
        *cutflow << TreeName("MUON_ID_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_MS_UP);
        *cutflow << TreeName("MUON_MS_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon MS (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MUON_MS_DN);
        *cutflow << TreeName("MUON_MS_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (UP)"); {
        *cutflow << EventSystematic(NtSys::MUON_SCALE_UP);
        *cutflow << TreeName("MUON_SCALE_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("muon scale shift (DN)"); {
        *cutflow << EventSystematic(NtSys::MUON_SCALE_DN);
        *cutflow << TreeName("MUON_SCALE_DN");
        *cutflow << SaveSystematic();
    }

    // jet
    *cutflow << NewSystematic("JER"); {
        *cutflow << EventSystematic(NtSys::JER);
        *cutflow << TreeName("JER");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
        *cutflow << TreeName("JET_GroupedNP_1_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 1 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
        *cutflow << TreeName("JET_GroupedNP_1_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 2 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_UP);
        *cutflow << TreeName("JET_GroupedNP_2_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 2 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_DN);
        *cutflow << TreeName("JET_GroupedNP_2_DN");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 3 (up)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_UP);
        *cutflow << TreeName("JET_GroupedNP_3_UP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JES NP set 3 (down)"); {
        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_DN);
        *cutflow << TreeName("JET_GroupedNP_3_DN");
        *cutflow << SaveSystematic();
    }
    //#warning only setting MET systematics
    // met
    *cutflow << NewSystematic("MET TST Soft-Term resolution (parallel)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
        *cutflow << TreeName("MET_SoftTrk_ResoPara");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term resolution (perpendicular)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
        *cutflow << TreeName("MET_SoftTrk_ResoPerp");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term shift in scale (UP)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
        *cutflow << TreeName("MET_SoftTrk_ScaleUp");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET TST Soft-Term shift in scale (DOWN)"); {
        *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
        *cutflow << TreeName("MET_SoftTrk_ScaleDown");
        *cutflow << SaveSystematic();
    }


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Superflow methods [END]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // initialize the cutflow and start the event loop
    chain->Process(cutflow, sample_.c_str(), num_events_, n_skip_);
    delete cutflow;
    delete chain;
    cout << "La Fin." << endl;
    exit(0);

}

void read_options(int argc, char* argv[], TChain* chain, int& n_skip_, int& num_events_, string& sample_,
    string& suffix_name_, SuperflowRunMode& run_mode_, SusyNtSys& nt_sys, bool& do_sumw_split, bool dbg)
{
    bool nominal_ = false;
    bool nominal_and_weight_sys_ = false;
    bool all_sys_ = false;
    string input;

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-n") == 0)
            num_events_ = atoi(argv[++i]);
        else if(strcmp(argv[i], "-i") == 0) {
            input = argv[++i];
        }
        else if(strcmp(argv[i], "-c") == 0) {
            nominal_ = true;
        }
        else if(strcmp(argv[i], "-w") == 0) {
            nominal_and_weight_sys_ = true;
        }
        else if(strcmp(argv[i], "-a") == 0) {
            all_sys_ = true;
        }
        else if(strcmp(argv[i], "-d") == 0) {
            dbg = true;
        }
        else if(strcmp(argv[i], "--suffix") == 0) {
            suffix_name_ = argv[++i];
        }
        else if(strcmp(argv[i], "--sumw") == 0) {
            do_sumw_split = true;
        }
        else {
            cout << analysis_name << "    Error (fatal) : Bad arguments." << endl;
            exit(1);
        }
    } // i

    sample_ = input;

    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    cout << analysis_name << "    input: " << input << endl;
    ChainHelper::addInput(chain, input, dbg); 

    Long64_t tot_num_events = chain->GetEntries();
    num_events_ = (num_events_ < 0 ? tot_num_events : num_events_);
    if(nominal_) {
        run_mode_ = SuperflowRunMode::nominal;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal" << endl;
    }
    if(nominal_and_weight_sys_) {
        run_mode_ = SuperflowRunMode::nominal_and_weight_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if(all_sys_) {
        run_mode_ = SuperflowRunMode::all_syst;
        cout << analysis_name << "    run mode: SuperflowRunMode::all_syst" << endl;
    }
    if(suffix_name_!="") {
        cout << analysis_name << "    Setting output file suffix to: " << suffix_name_ << endl;
    }
}

