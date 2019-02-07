// std
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <math.h>
#include <random>

// ROOT
#include "TChain.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TF1.h"

// SusyNtuple
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/SusyNtSys.h"
#include "SusyNtuple/KinematicTools.h"

// Superflow
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "Superflow/Cut.h"
#include "Superflow/StringTools.h"
#include "Superflow/input_options.h"

// lwtnn
#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/parse_json.hh"

// DileptonTriggerTool
#include "DileptonTriggerTool/DileptonTriggerTool.h"

// HH Weight
#include "hhTruthWeightTools/hhWeightTool.h"

// LambdaTool
#include "HHLambdaTool/HHLambdaTool.h"

using namespace std;
using namespace sflow;

const string analysis_name = "ntupler_muj";

//string network_dir = "/data/uclhc/uci/user/dantrim/n0303val/susynt-read/data/";
string network_dir = "./susynt-read/data/";
string nn_file = network_dir + "nn_descriptor_nombbmt2_1k.json";
string nn_file_update = network_dir + "lwtnn_NN_jan15_mt2bb_3x_0p5_250node_adam.json";

int main(int argc, char* argv[])
{
    /////////////////////////////////////////////////////////////////////
    // Read in the command-line options (input file, num events, etc...)
    /////////////////////////////////////////////////////////////////////
    SFOptions options(argc, argv);
    options.ana_name = analysis_name;
    if(!read_options(options)) {
        exit(1);
    }


    // load the NN
    if(!(std::ifstream(nn_file).good())) {
        cout << options.ana_name << "    ERROR Input NN file (=" << nn_file << ") not found!" << endl;
        exit(1);
    }
    std::ifstream input_nn_file(nn_file);
    std::string output_layer_name = "OutputLayer";
    auto config = lwt::parse_json_graph(input_nn_file);
    lwt::LightweightGraph nn_graph(config, output_layer_name);
    std::map< std::string, double > nn_input;

    // load updated NN
    if(!(std::ifstream(nn_file_update).good()))
    {
        cout << options.ana_name << "    ERROR Input NN file (=" << nn_file_update << ") not found!" << endl;
        exit(1);
    }
    std::ifstream input_nn_file_update(nn_file_update);
    std::string output_layer_name_update = "OutputLayer";
    auto config_update = lwt::parse_json_graph(input_nn_file_update);
    lwt::LightweightGraph nn_graph_update(config_update, output_layer_name_update);

    std::map<std::string, float> var_means;
    var_means["j0_pt"] = 134.898;
    var_means["j1_pt"] = 69.5677;
    var_means["bj0_pt"] = 122.180;
    var_means["bj1_pt"] = 54.1126;
    var_means["j0_eta"] = 0.0035159;
    var_means["j1_eta"] = -0.0014209;
    var_means["bj0_eta"] = -0.005168;
    var_means["bj1_eta"] = 0.00638;
    var_means["j0_phi"] = 0.0146455;
    var_means["j1_phi"] = 0.0051678;
    var_means["bj0_phi"] = 0.013698;
    var_means["bj1_phi"] = 0.011199;
    var_means["dphi_j0_ll"] = 0.0058626;
    var_means["dphi_j0_l0"] = -0.0030659;
    var_means["dphi_bj0_ll"] = 0.0086884;
    var_means["dphi_bj0_l0"] = -0.0016912;
    var_means["mbb"] = 144.59769;
    var_means["dRbb"] = 2.130620;
    var_means["dR_ll_bb"] = 2.815526;
    var_means["dphi_ll_bb"] = 0.00045855;
    var_means["dphi_WW_bb"] = -0.0093672;
    var_means["HT2"] = 279.0936;
    var_means["HT2Ratio"] = 0.63980;
    var_means["MT_1"] = 478.057;
    var_means["MT_1_scaled"] = 470.3389;
    var_means["mt2_llbb"] = 172.9586;
    var_means["mt2_bb"] = 66.25853;
    var_means["dphi_bb"] = -0.003595;
    var_means["mT_bb"] = 144.5976;

    // HH reweighting
    //cout << analysis_name << "    Initializing HH weight tool" << endl;
    //xAOD::hhWeightTool hhWeightTool("hhWeights");
    //hhWeightTool.setProperty("ReweightFile", "hhTruthWeightTools/SMhh_mhh_ReWeight.root");
    //hhWeightTool.initialize();

    //// lambda scan
    //cout << analysis_name << "    Initializing lambda scan tool" << endl;
    //hh::HHLambdaTool lambdaTool;
    //// load the default, which are the SM histograms
    //lambdaTool.load();
    
    

    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    bool verbose = true;
    ChainHelper::addInput(chain, options.input, verbose);
    Long64_t tot_num_events = chain->GetEntries();
    options.n_events_to_process = (options.n_events_to_process < 0 ? tot_num_events : options.n_events_to_process);

    ////////////////////////////////////////////////////
    // Construct and configure the Superflow object
    ////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow();
    cutflow->setAnaName(options.ana_name);
    cutflow->setAnaType(AnalysisType::Ana_WWBB);

    float lumi_to_set_in_pb = 1000.;
    cutflow->setLumi(lumi_to_set_in_pb); // 1/fb
    cutflow->setSampleName(options.input);
    cutflow->setRunMode(options.run_mode);
    cutflow->setCountWeights(true);
    cutflow->setChain(chain);
    cutflow->setDebug(options.dbg);
    if(options.suffix_name != "") {
        cutflow->setFileSuffix(options.suffix_name);
    }
    if(options.sumw_file_name != "") {
        cout << options.ana_name << "    Reading sumw for sample from file: " << options.sumw_file_name << endl;
        cutflow->setUseSumwFile(options.sumw_file_name);
    }
    cutflow->nttools().initTriggerTool(ChainHelper::firstFile(options.input, options.dbg));

    // trigger tool
    dileptrig::DileptonTriggerTool* dilep_trig_tool1516 = new dileptrig::DileptonTriggerTool();
    dilep_trig_tool1516->initialize("WWBB1516", &cutflow->nttools().triggerTool());
    dileptrig::DileptonTriggerTool* dilep_trig_tool151617 = new dileptrig::DileptonTriggerTool();
    dilep_trig_tool151617->initialize("WWBB151617RAND", &cutflow->nttools().triggerTool());
    dileptrig::DileptonTriggerTool* dilep_trig_tool15161718 = new dileptrig::DileptonTriggerTool();
    dilep_trig_tool15161718->initialize("WWBB15161718RAND", &cutflow->nttools().triggerTool());

    // print some useful
    cout << analysis_name << "    Total Entries    : " << chain->GetEntries() << endl;
    if(options.n_events_to_process > 0) {
    cout << analysis_name << "    Process Entries  : " << options.n_events_to_process << endl;
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

    *cutflow << CutName("read in ") << [](Superlink* /* sl */) -> bool { return true; };

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
    *cutflow << CutName("==2 signal leptons") << [](Superlink* sl) -> bool {
        return sl->leptons->size() == 2;
    };

    //#warning REMOVING LEPTON PT CUT
  /*
    *cutflow << CutName("lepton pTs > (25,20) GeV") << [](Superlink* sl) -> bool {
        return ( (sl->leptons->at(0)->Pt()>25) && (sl->leptons->at(1)->Pt()>20) );
    };
  */

    *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
        return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
    };

  ///*
    *cutflow << CutName("mll > 20 GeV") << [](Superlink* sl) -> bool {
        return ( (*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 20. );
    };
    //*/

    //*cutflow << CutName("veto SF Z-window (within 20 GeV)") << [](Superlink* sl) -> bool {
    //    bool pass = true;
    //    bool isSF = false;
    //    if((sl->leptons->size()==2 && (sl->electrons->size()==2 || sl->muons->size()==2))) isSF = true;
    //    if(isSF) {
    //        double mll = (*sl->leptons->at(0) + *sl->leptons->at(1)).M();
    //        if( fabs(mll-91.2) < 20. ) pass = false;
    //    }
    //    return pass;
    //};

   *cutflow << CutName(">=2 b-tagged jets") << [](Superlink* sl) -> bool {
       int n_bjets = 0;
       for(auto j : (*sl->jets)) {
           if(sl->tools->jetSelector().isBJet(j)) { n_bjets++; }
       }
       if(n_bjets>=2) return true;
       else return false;
   };

    ///*
//    *cutflow << CutName("pass trigger requirement") << [&](Superlink* sl) -> bool {
//        int year = sl->nt->evt()->treatAsYear;
//        bool pass = false;
//
//        bool passes_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
//        bool passes_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
//        bool passes_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
//        bool passes_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
//        bool passes_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
//        bool passes_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14");
//
//        bool trig_pass_2015 = (passes_mu18_mu8noL1 || passes_e17_lhloose_mu14 || passes_2e12_lhloose_L12EM10VH);
//        bool trig_pass_2016 = (passes_2e17_lhvloose_nod0 || passes_mu22_mu8noL1 || passes_e17_lhloose_nod0_mu14);
//        if( (year==2015 && trig_pass_2015==1) || (year==2016 && trig_pass_2016==1) ) {
//            pass = true;
//        }
//        return pass;
//    };
    //*/

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


    // HH signal variables
    float hh_truth_mHH = -1.;
    //std::map<int, float> lambda_weights;
    //static bool print_lambdas = true;
    //*cutflow << [&](Superlink* sl, var_void*) {
    //    int dsid = sl->nt->evt()->mcChannel;
    //    if(dsid == 346218)
    //    {
    //        hh_truth_mHH = sl->nt->evt()->truth_mhh;
    //        lambda_weights = lambdaTool.get_weights(hh_truth_mHH);
    //        if(print_lambdas)
    //        {
    //            cout << analysis_name << "    Loading weights for lambda_hhh values: ";
    //            for(auto & l : lambda_weights)
    //            {
    //                cout << " " << l.first;
    //            }
    //            cout << endl;
    //            print_lambdas = false;
    //        }
    //    }
    //};


    *cutflow << NewVar("WWbb truth mHH"); {
        *cutflow << HFTname("truth_mHH");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            if(sl->nt->evt()->mcChannel == 346218)
            {
                return hh_truth_mHH;
            }
            else
            {
                return -1.0;
            }
        };
        *cutflow << SaveVar();
    }

    //*cutflow << NewVar("WWbb NLO reweight"); {
    //    *cutflow << HFTname("hh_NLO_weight");
    //    *cutflow << [&](Superlink* sl, var_float*) -> double {
    //        if(sl->nt->evt()->mcChannel == 346218)
    //        {
    //            return hhWeightTool.getWeight( hh_truth_mHH * 1000. ); // the tool requires values in MeV
    //        }
    //        return 1.0;
    //    };
    //    *cutflow << SaveVar();
    //}
//
//    *cutflow << NewVar("WWbb lambda weight (-20)"); {
//        *cutflow << HFTname("hh_lambda_w_neg20");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-20);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-19)"); {
//        *cutflow << HFTname("hh_lambda_w_neg19");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-19);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-18)"); {
//        *cutflow << HFTname("hh_lambda_w_neg18");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-18);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-17)"); {
//        *cutflow << HFTname("hh_lambda_w_neg17");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-17);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-16)"); {
//        *cutflow << HFTname("hh_lambda_w_neg16");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-16);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-15)"); {
//        *cutflow << HFTname("hh_lambda_w_neg15");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-15);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-14)"); {
//        *cutflow << HFTname("hh_lambda_w_neg14");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-14);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-13)"); {
//        *cutflow << HFTname("hh_lambda_w_neg13");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-13);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-12)"); {
//        *cutflow << HFTname("hh_lambda_w_neg12");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-12);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-11)"); {
//        *cutflow << HFTname("hh_lambda_w_neg11");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-11);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-10)"); {
//        *cutflow << HFTname("hh_lambda_w_neg10");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-10);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-9)"); {
//        *cutflow << HFTname("hh_lambda_w_neg9");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-9);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-8)"); {
//        *cutflow << HFTname("hh_lambda_w_neg8");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-8);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-7)"); {
//        *cutflow << HFTname("hh_lambda_w_neg7");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-7);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-6)"); {
//        *cutflow << HFTname("hh_lambda_w_neg6");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-6);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-5)"); {
//        *cutflow << HFTname("hh_lambda_w_neg5");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-5);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-4)"); {
//        *cutflow << HFTname("hh_lambda_w_neg4");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-4);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-3)"); {
//        *cutflow << HFTname("hh_lambda_w_neg3");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-3);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-2)"); {
//        *cutflow << HFTname("hh_lambda_w_neg2");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-2);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (-1)"); {
//        *cutflow << HFTname("hh_lambda_w_neg1");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(-1);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (0)"); {
//        *cutflow << HFTname("hh_lambda_w_pos0");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(0);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (1)"); {
//        *cutflow << HFTname("hh_lambda_w_pos1");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(1);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (2)"); {
//        *cutflow << HFTname("hh_lambda_w_pos2");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(2);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (3)"); {
//        *cutflow << HFTname("hh_lambda_w_pos3");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(3);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (4)"); {
//        *cutflow << HFTname("hh_lambda_w_pos4");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(4);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (5)"); {
//        *cutflow << HFTname("hh_lambda_w_pos5");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(5);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (6)"); {
//        *cutflow << HFTname("hh_lambda_w_pos6");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(6);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (7)"); {
//        *cutflow << HFTname("hh_lambda_w_pos7");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(7);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (8)"); {
//        *cutflow << HFTname("hh_lambda_w_pos8");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(8);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (9)"); {
//        *cutflow << HFTname("hh_lambda_w_pos9");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(9);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (10)"); {
//        *cutflow << HFTname("hh_lambda_w_pos10");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(10);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (11)"); {
//        *cutflow << HFTname("hh_lambda_w_pos11");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(11);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (12)"); {
//        *cutflow << HFTname("hh_lambda_w_pos12");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(12);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (13)"); {
//        *cutflow << HFTname("hh_lambda_w_pos13");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(13);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (14)"); {
//        *cutflow << HFTname("hh_lambda_w_pos14");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(14);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (15)"); {
//        *cutflow << HFTname("hh_lambda_w_pos15");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(15);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (16)"); {
//        *cutflow << HFTname("hh_lambda_w_pos16");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(16);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (17)"); {
//        *cutflow << HFTname("hh_lambda_w_pos17");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(17);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (18)"); {
//        *cutflow << HFTname("hh_lambda_w_pos18");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(18);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (19)"); {
//        *cutflow << HFTname("hh_lambda_w_pos19");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(19);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }
//    
//    
//    *cutflow << NewVar("WWbb lambda weight (20)"); {
//        *cutflow << HFTname("hh_lambda_w_pos20");
//        *cutflow << [&](Superlink* sl, var_float*) -> double {
//            if(sl->nt->evt()->mcChannel == 346218)
//            {
//                return lambda_weights.at(20);
//            }
//            return 1.0;
//        };
//        *cutflow << SaveVar();
//    }



    bool p_mu20_iloose_L1MU15;
    bool p_mu20_ivarloose_L1MU15;
    bool p_mu26_ivarmedium;
    bool p_mu18_mu8noL1;
    bool p_mu22_mu8noL1;
    bool p_e24_lhmedium_L1EM20VH;
    bool p_e24_lhmedium_L1EM20VHI;
    bool p_e24_lhtight_nod0_ivarloose;
    bool p_e26_lhtight_nod0_ivarloose;
    bool p_2e12_lhloose_L12EM10VH;
    bool p_2e17_lhvloose_nod0;
    bool p_2e17_lhvloose_nod0_L12EM15VHI;
    bool p_2e24_lhvloose_nod0;
    bool p_e7_lhmedium_nod0_mu24;
    bool p_e7_lhmedium_mu24;
    bool p_e17_lhloose_mu14;
    bool p_e17_lhloose_nod0_mu14;
    bool p_e24_lhmedium_nod0_L1EM20VHI_mu8noL1;
    bool p_e26_lhmedium_nod0_mu8noL1;
    bool p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;

    *cutflow << [&](Superlink* sl, var_void*) {
    	p_mu20_iloose_L1MU15 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_iloose_L1MU15");
    	p_mu20_ivarloose_L1MU15 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_ivarloose_L1MU15");
    	p_mu26_ivarmedium = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu26_ivarmedium");
    	p_mu18_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1");
    	p_mu22_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1");
    	p_e24_lhmedium_L1EM20VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_L1EM20VH");
    	p_e24_lhmedium_L1EM20VHI = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_L1EM20VHI");
    	p_e24_lhtight_nod0_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhtight_nod0_ivarloose");
    	p_e26_lhtight_nod0_ivarloose = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhtight_nod0_ivarloose");
    	p_2e12_lhloose_L12EM10VH = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH");
    	p_2e17_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0");
    	p_2e17_lhvloose_nod0_L12EM15VHI = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0_L12EM15VHI");
    	p_2e24_lhvloose_nod0 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e24_lhvloose_nod0");
    	p_e7_lhmedium_nod0_mu24 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e7_lhmedium_nod0_mu24");
    	p_e7_lhmedium_mu24 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e7_lhmedium_mu24");
    	p_e17_lhloose_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14");
    	p_e17_lhloose_nod0_mu14 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14");
    	p_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1");
    	p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1");
    	p_e26_lhmedium_nod0_mu8noL1 = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e26_lhmedium_nod0_mu8noL1");
    };
    *cutflow << NewVar("pass mu20_iloose_L1MU15"); {
    	*cutflow << HFTname("trig_mu20_iloose_L1MU15");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_mu20_iloose_L1MU15;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu20_ivarloose_L1MU15"); {
    	*cutflow << HFTname("trig_mu20_ivarloose_L1MU15");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_mu20_ivarloose_L1MU15;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu26_ivarmedium"); {
    	*cutflow << HFTname("trig_mu26_ivarmedium");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_mu26_ivarmedium;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu18_mu8noL1"); {
    	*cutflow << HFTname("trig_mu18_mu8noL1");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_mu18_mu8noL1;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass mu22_mu8noL1"); {
    	*cutflow << HFTname("trig_mu22_mu8noL1");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_mu22_mu8noL1;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e24_lhmedium_L1EM20VH"); {
    	*cutflow << HFTname("trig_e24_lhmedium_L1EM20VH");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e24_lhmedium_L1EM20VH;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e24_lhmedium_L1EM20VHI"); {
    	*cutflow << HFTname("trig_e24_lhmedium_L1EM20VHI");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e24_lhmedium_L1EM20VHI;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e24_lhtight_nod0_ivarloose"); {
    	*cutflow << HFTname("trig_e24_lhtight_nod0_ivarloose");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e24_lhtight_nod0_ivarloose;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e26_lhtight_nod0_ivarloose"); {
    	*cutflow << HFTname("trig_e26_lhtight_nod0_ivarloose");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e26_lhtight_nod0_ivarloose;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e12_lhloose_L12EM10VH"); {
    	*cutflow << HFTname("trig_2e12_lhloose_L12EM10VH");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_2e12_lhloose_L12EM10VH;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e17_lhvloose_nod0"); {
    	*cutflow << HFTname("trig_2e17_lhvloose_nod0");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_2e17_lhvloose_nod0;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e17_lhvloose_nod0_L12EM15VHI"); {
    	*cutflow << HFTname("trig_2e17_lhvloose_nod0_L12EM15VHI");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_2e17_lhvloose_nod0_L12EM15VHI;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass 2e24_lhvloose_nod0"); {
    	*cutflow << HFTname("trig_2e24_lhvloose_nod0");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_2e24_lhvloose_nod0;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e7_lhmedium_nod0_mu24"); {
    	*cutflow << HFTname("trig_e7_lhmedium_nod0_mu24");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e7_lhmedium_nod0_mu24;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e7_lhmedium_mu24"); {
    	*cutflow << HFTname("trig_e7_lhmedium_mu24");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e7_lhmedium_mu24;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e17_lhloose_mu14"); {
    	*cutflow << HFTname("trig_e17_lhloose_mu14");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e17_lhloose_mu14;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e17_lhloose_nod0_mu14"); {
    	*cutflow << HFTname("trig_e17_lhloose_nod0_mu14");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e17_lhloose_nod0_mu14;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e24_lhmedium_nod0_L1EM20VHI_mu8noL1"); {
    	*cutflow << HFTname("trig_e24_lhmedium_nod0_L1EM20VHI_mu8noL1");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e24_lhmedium_nod0_L1EM20VHI_mu8noL1;
    	};
    	*cutflow << SaveVar();
    }
    *cutflow << NewVar("pass e26_lhmedium_nod0_L1EM22VHI_mu8noL1"); {
    	*cutflow << HFTname("trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1");
    	*cutflow << [&](Superlink* /*sl*/, var_bool*) -> bool {
    		return p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;
    	};
    	*cutflow << SaveVar();
    }

    *cutflow << NewVar("pass_wwbb_trig_1516"); {
        *cutflow << HFTname("wwbb_trig_1516");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return dilep_trig_tool1516->pass_trigger_selection(sl->nt->evt(),
                    sl->leptons->at(0), sl->leptons->at(1), false);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass_wwbb_trig_151617"); {
        *cutflow << HFTname("wwbb_trig_151617Rand");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return dilep_trig_tool151617->pass_trigger_selection(sl->nt->evt(),
                    sl->leptons->at(0), sl->leptons->at(1), false);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("pass_wwbb_trig_151617"); {
        *cutflow << HFTname("wwbb_trig_15161718Rand");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            return dilep_trig_tool15161718->pass_trigger_selection(sl->nt->evt(),
                    sl->leptons->at(0), sl->leptons->at(1), false);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass_trig_tight_2015"); {

        *cutflow << HFTname("trig_tight_2015");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt =  sl->leptons->at(1)->Pt();

            if(isEE) {
                if(lead_pt>=26 && p_e24_lhmedium_L1EM20VH) {
                    return true;
                }
                else if( (lead_pt>=14 && sub_pt>=14) && p_2e12_lhloose_L12EM10VH) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isMM) {
                if(lead_pt>=22 && p_mu20_iloose_L1MU15) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=10) && p_mu18_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isEM) {
                if(lead_pt>=26 && p_e24_lhmedium_L1EM20VH) {
                    return true;
                }
                //else if( (lead_pt>=26 && sub_pt>=10) && p_e24_lhmedium_nod0_L1EM20VHI_mu8noL1) {
                //    return true;
                //}
                else if( (lead_pt>=19 && sub_pt>=16) && p_e17_lhloose_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isME) {
                if(lead_pt>=22 && p_mu20_iloose_L1MU15) {
                    return true;
                }
                else if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "ERROR in pass tight 2015!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass_trig_tight_2015dil"); {

        *cutflow << HFTname("trig_tight_2015dil");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt =  sl->leptons->at(1)->Pt();

            if(isEE) {
                if( (lead_pt>=14 && sub_pt>=14) && p_2e12_lhloose_L12EM10VH) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isMM) {
                if( (lead_pt>=20 && sub_pt>=10) && p_mu18_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isEM) {
                //if( (lead_pt>=26 && sub_pt>=10) && p_e24_lhmedium_nod0_L1EM20VHI_mu8noL1) {
                //    return true;
                //}
                if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isME) {
                if( (lead_pt>=26 && sub_pt>=10) && p_e7_lhmedium_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "ERROR in pass tight 2015!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass_trig_tight_2016"); {

        *cutflow << HFTname("trig_tight_2016");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt =  sl->leptons->at(1)->Pt();

            if(isEE) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isMM) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isEM) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isME) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2016 error!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass_trig_tight_2016dil"); {

        *cutflow << HFTname("trig_tight_2016dil");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() &&  sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt =  sl->leptons->at(1)->Pt();

            if(isEE) {
                if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isMM) {
                if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isEM) {
                if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_L1EM22VHI_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if(isME) {
                if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2016 error!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass trig tight 2017"); {
        *cutflow << HFTname("trig_tight_2017");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            /////////////////////////////////////
            if(isEE) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                    return true;
                }
                else {
                    return false;
                }
            } // EE
            /////////////////////////////////////
            else if(isMM) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            } // MM
            /////////////////////////////////////
            else if(isEM) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            /////////////////////////////////////
            else if(isME) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2017 error!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<float> uniform_distribution(0.0, 1.0);

    *cutflow << NewVar("pass trig tight 2017 with random"); {
        *cutflow << HFTname("trig_tight_2017rand");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            int run_number = sl->nt->evt()->run;
            bool isMC = sl->nt->evt()->isMC;

            /////////////////////////////////////
            if(isEE) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if(isMC) {
                    float random_number = uniform_distribution(generator);
                    if(random_number < (0.6 / 140.)) {
                    //if(random_number < (0.6 / 78.2)) {
                        if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                            return true;
                        }
                    }
                    else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
                else if(!isMC) {
                    if(run_number>=326834 && run_number<=328393) {
                        if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                            return true;
                        }
                    }
                    else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
            } // EE
            /////////////////////////////////////
            else if(isMM) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            } // MM
            /////////////////////////////////////
            else if(isEM) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            /////////////////////////////////////
            else if(isME) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2017rand error!" << endl;
                exit(1);
            }
            return false;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass trig tight 2017dil"); {
        *cutflow << HFTname("trig_tight_2017dil");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            /////////////////////////////////////
            if(isEE) {
                if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                    return true;
                }
                else {
                    return false;
                }
            } // EE
            /////////////////////////////////////
            else if(isMM) {
                if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            } // MM
            /////////////////////////////////////
            else if(isEM) {
                if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            /////////////////////////////////////
            else if(isME) {
                if( (lead_pt>=26 && sub_pt>=10) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2017dil error!" << endl;
                exit(1);
            }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass trig tight 2017 dilepton only with random"); {
        *cutflow << HFTname("trig_tight_2017dilrand");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());

            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            int run_number = sl->nt->evt()->run;
            bool isMC = sl->nt->evt()->isMC;

            /////////////////////////////////////
            if(isEE) {
                if(isMC) {
                    float random_number = uniform_distribution(generator);
                    if(random_number < (0.6 / 140.)) {
                        if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                            return true;
                        }
                    }
                    else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
                else if(!isMC) {
                    if(run_number>=326834 && run_number<=328393) {
                        if( (lead_pt>=26 && sub_pt>=26) && p_2e24_lhvloose_nod0) {
                            return true;
                        }
                    }
                    else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
            } // EE
            /////////////////////////////////////
            else if(isMM) {
                if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1) {
                    return true;
                }
                else {
                    return false;
                }
            } // MM
            /////////////////////////////////////
            else if(isEM) {
                if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            /////////////////////////////////////
            else if(isME) {
                if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else {
                cout << "trig pass tight 2017dilrand error!" << endl;
                exit(1);
            }
            return false;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass trig tight 2018"); {
        *cutflow << HFTname("trig_tight_2018");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());
            
            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            //////////////////////////////////////////////
            if(isEE) {
                if(lead_pt >= 28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                    return true;
                }
                else {
                    return false;
                }
            } // isEE
            //////////////////////////////////////////////
            else if(isMM) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1 ) {
                    return true;
                }
                else {
                    return false;
                }
            } // isMM
            //////////////////////////////////////////////
            else if(isEM) {
                if(lead_pt>=28 && p_e26_lhtight_nod0_ivarloose) {
                    return true;
                }
                else if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            } // isEM
            //////////////////////////////////////////////
            else if(isME) {
                if(lead_pt>=28 && p_mu26_ivarmedium) {
                    return true;
                }
                else if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24 ) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            } // isME
            else {
                cout << "tright pass tight 2018 error!" << endl;
                exit(1);
            }
            return false;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pass trig tight 2018 (dilepton only)"); {
        *cutflow << HFTname("trig_tight_2018dil");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool {
            if(sl->leptons->size()<2) return false;
            bool isEE = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
            bool isMM = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu());
            bool isEM = (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu());
            bool isME = (sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle());
            
            float lead_pt = sl->leptons->at(0)->Pt();
            float sub_pt = sl->leptons->at(1)->Pt();

            //////////////////////////////////////////////
            if(isEE) {
                if( (lead_pt>=19 && sub_pt>=19) && p_2e17_lhvloose_nod0_L12EM15VHI) {
                    return true;
                }
                else {
                    return false;
                }
            } // isEE
            //////////////////////////////////////////////
            else if(isMM) {
                if( (lead_pt>=24 && sub_pt>=10) && p_mu22_mu8noL1 ) {
                    return true;
                }
                else {
                    return false;
                }
            } // isMM
            //////////////////////////////////////////////
            else if(isEM) {
                if( (lead_pt>=28 && sub_pt>=10) && p_e26_lhmedium_nod0_mu8noL1) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            } // isEM
            //////////////////////////////////////////////
            else if(isME) {
                if( (lead_pt>=26 && sub_pt>=9) && p_e7_lhmedium_nod0_mu24 ) {
                    return true;
                }
                else if( (lead_pt>=20 && sub_pt>=17) && p_e17_lhloose_nod0_mu14) {
                    return true;
                }
                else {
                    return false;
                }
            } // isME
            else {
                cout << "tright pass tight 2018 error!" << endl;
                exit(1);
            }
            return false;
        };
        *cutflow << SaveVar();
    }
    
    *cutflow << NewVar("eventNumber"); {
        *cutflow << HFTname("eventNumber");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->eventNumber;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("run"); {
        *cutflow << HFTname("runNumber");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->run;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lumi block"); {
        *cutflow << HFTname("lumi_block");
        *cutflow << [](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->lb;
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

    *cutflow << NewVar("mc campaign (Susy::MCType)"); {
        *cutflow << HFTname("mcType");
        *cutflow << [](Superlink* sl, var_int*) -> int {
            return sl->nt->evt()->mcType;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("year"); {
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

    *cutflow << NewVar("event weight without pileup weight"); {
        *cutflow << HFTname("eventweightNoPRW");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF"); {
        *cutflow << HFTname("eventweightbtag");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->nt->evt()->wPileup * sl->weights->btagSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF NoPRW"); {
        *cutflow << HFTname("eventweightbtagNoPRW");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->weights->btagSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF x jvtSf"); {
        *cutflow << HFTname("eventweightBtagJvt");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->nt->evt()->wPileup * sl->weights->btagSf * sl->weights->jvtSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF x jvtSf NoPRW"); {
        *cutflow << HFTname("eventweightBtagJvtNoPRW");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product() * sl->weights->btagSf * sl->weights->jvtSf;
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

    *cutflow << NewVar("event weight (multi period)"); {
        *cutflow << HFTname("eventweight_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi() * sl->nt->evt()->wPileup;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight without pileup weight"); {
        *cutflow << HFTname("eventweightNoPRW_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF"); {
        *cutflow << HFTname("eventweightbtag_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi() * sl->nt->evt()->wPileup * sl->weights->btagSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF NoPRW"); {
        *cutflow << HFTname("eventweightbtagNoPRW_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi() * sl->weights->btagSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF x jvtSf"); {
        *cutflow << HFTname("eventweightBtagJvt_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi() * sl->nt->evt()->wPileup * sl->weights->btagSf * sl->weights->jvtSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight x btag SF x jvtSf NoPRW"); {
        *cutflow << HFTname("eventweightBtagJvtNoPRW_multi");
        *cutflow << [&](Superlink* sl, var_double*) -> double {
            return sl->weights->product_multi() * sl->weights->btagSf * sl->weights->jvtSf;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("pile-up weight with period weight divided out"); {
        *cutflow << HFTname("pupwNoPeriod");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return (sl->nt->evt()->wPileup / sl->nt->evt()->wPileup_period);
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

    *cutflow << NewVar("Pile-up weight period weight"); {
        *cutflow << HFTname("pupw_period");
        *cutflow << [](Superlink* sl, var_double*) -> double {
            return sl->nt->evt()->wPileup_period;
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

    *cutflow << NewVar("actual interactions per b.c."); {
        *cutflow << HFTname("actualMu");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->actualMu; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("actual interactions per b.c. with data scale factor applied"); {
        *cutflow << HFTname("actualMuDataSF");
        *cutflow << [](Superlink* sl, var_float*) -> double { return sl->nt->evt()->actualMuDataSF; };
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
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return leptons.size(); };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("number of electrons"); {
       *cutflow << HFTname("nElectrons");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return electrons.size(); };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("number of muons"); {
       *cutflow << HFTname("nMuons");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return muons.size(); };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("is an EE event"); {
       *cutflow << HFTname("isEE");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int  {
            if(leptons.size()<2) return 0;
            int val = 0;
            if(leptons.at(0)->isEle() && leptons.at(1)->isEle()) { val = 1; }
            else { val = 0; }
            nn_input["isEE"] = val;
            return val;
       };
       *cutflow << SaveVar();
   }

   *cutflow << NewVar("is an MM event"); {
       *cutflow << HFTname("isMM");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(leptons.size()<2) return 0;
            int val = 0;
            if(leptons.at(0)->isMu() && leptons.at(1)->isMu()) { val = 1; }
            else { val = 0; }
            nn_input["isMM"] = val;
            return val;
       };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("is an EM event"); {
       *cutflow << HFTname("isEM");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(leptons.size()<2) return 0;
            int val = 0;
            if(leptons.at(0)->isEle() && leptons.at(1)->isMu()) { val = 1; }
            else { val = 0; }
            nn_input["isEM"] = val;
            return val;
       };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("is an ME event"); {
       *cutflow << HFTname("isME");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(leptons.size()<2) return 0;
            int val = 0;
            if(leptons.at(0)->isMu() && leptons.at(1)->isEle()) { val = 1; }
            else { val = 0; }
            nn_input["isME"] = val;
            return val;
       };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("is a SF event"); {
       *cutflow << HFTname("isSF");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(leptons.size()<2) return 0;
             int val = 0;
            if( (leptons.at(0)->isEle() && leptons.at(1)->isEle()) ||
                (leptons.at(0)->isMu() && leptons.at(1)->isMu()) ) { val = 1; }
            else { val = 0; }
             nn_input["isSF"] = val;
             return val;
       };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("is a DF event"); {
       *cutflow << HFTname("isDF");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
           if(leptons.size()<2) return 0;
            int val = 0;
           if( (leptons.at(0)->isEle() && leptons.at(1)->isMu()) ||
               (leptons.at(0)->isMu() && leptons.at(1)->isEle()) )  { val = 1; }
           else { val = 0; }
            nn_input["isDF"] = val;
            return val;
       };
       *cutflow << SaveVar();
   }
   *cutflow << NewVar("lepton flavor [EE=0,MM=1,EM=2,ME=3]"); {
       *cutflow << HFTname("l_flav");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
           if(leptons.size()<2) return -1;
           bool e0 = leptons.at(0)->isEle();
           bool e1 = leptons.at(1)->isEle();

           if( e0 && e1 ) { return 0; }
           else if( !e0 && !e1 ) { return 1; }
           else if( e0 && !e1 ) { return 2; }
           else if( !e0 && e1 ) { return 3; }
           else { return -1; }
       };
       *cutflow << SaveVar();
   }

   *cutflow << NewVar("lead lepton flavor [E=0, M=1]"); {
       *cutflow << HFTname("l0_flav");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
           bool e = leptons.at(0)->isEle();
           bool m = leptons.at(0)->isMu();

           if(e && !m) return 0;
           if(!e && m) return 1;
           else { return -1; }
       };
       *cutflow << SaveVar();
   }

   *cutflow << NewVar("sub lead lepton flavor [E=0, M=1]"); {
       *cutflow << HFTname("l1_flav");
       *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
           if(leptons.size()<2) return -1;
           bool e = leptons.at(1)->isEle();
           bool m = leptons.at(1)->isMu();
           if(e && !m) return 0;
           if(!e && m) return 1;
           else { return -1; }
       };
       *cutflow << SaveVar();
   }

    *cutflow << NewVar("lead lepton q"); {
        *cutflow << HFTname("l0_q");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return leptons.at(0)->q; };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lepton q"); {
        *cutflow << HFTname("l1_q");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(leptons.size()<2) return 0;
            return leptons.at(1)->q;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lepton d0"); {
        *cutflow << HFTname("l0_d0");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return leptons.at(0)->d0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lepton d0"); {
        *cutflow << HFTname("l1_d0");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -10;
            return leptons.at(1)->d0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lepton d0sig"); {
        *cutflow << HFTname("l0_d0sig");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return leptons.at(0)->d0sigBSCorr; };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lepton d0sig"); {
        *cutflow << HFTname("l1_d0sig");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -10;
            return leptons.at(1)->d0sigBSCorr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lepton z0sinTheta"); {
        *cutflow << HFTname("l0_z0sinTheta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return leptons.at(0)->z0SinTheta(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lepton z0sinTheta"); {
        *cutflow << HFTname("l1_z0sinTheta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -10;
            return leptons.at(1)->z0SinTheta();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead lepton pt"); {
        *cutflow << HFTname("l0_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            nn_input["l0_pt"] = leptons.at(0)->Pt();
            return leptons.at(0)->Pt();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sublead lepton pt"); {
        *cutflow << HFTname("l1_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -1;
            nn_input["l1_pt"] = leptons.at(1)->Pt();
            return leptons.at(1)->Pt();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lep eta"); {
        *cutflow << HFTname("l0_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            nn_input["l0_eta"] = leptons.at(0)->Eta();
            return leptons.at(0)->Eta();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lep eta"); {
        *cutflow << HFTname("l1_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -5;
            nn_input["l1_eta"] = leptons.at(1)->Eta();
            return leptons.at(1)->Eta();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead lep phi"); {
        *cutflow << HFTname("l0_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            nn_input["l0_phi"] = leptons.at(0)->Phi();
            return leptons.at(0)->Phi();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sublead lep phi"); {
        *cutflow << HFTname("l1_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
           if(leptons.size()<2) return -5;
           nn_input["l1_phi"] = leptons.at(1)->Phi();
           return leptons.at(1)->Phi();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mll leptons"); {
        *cutflow << HFTname("mll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double mll = -10.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                mll = (*l0 + *l1).M();
            }
            nn_input["mll"] = mll;
            return mll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton pT"); {
        *cutflow << HFTname("pTll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double pTll = -10.0;
            if(leptons.size() == 2) {
                Lepton* l0 = leptons.at(0);
                Lepton* l1 = leptons.at(1);
                pTll = (*l0 + *l1).Pt();
            }
            nn_input["pTll"] = pTll;
            return pTll;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between to leptons"); {
        *cutflow << HFTname("dphi_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double dphi = -10.0;
            if(leptons.size() == 2) {
                Lepton l0 = *leptons.at(0);
                Lepton l1 = *leptons.at(1);
                dphi = l0.DeltaPhi(l1);
            }
            nn_input["dphi_ll"] = dphi;
            return dphi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta eta between two leptons"); {
        *cutflow << HFTname("deta_ll");
        *cutflow << [&](Superlink* /* sl */, var_float*) -> double {
            double deta = -10.0;
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
    //JetVector bjets;
    std::vector<Susy::Jet> bjets;
    JetVector sjets;

    *cutflow << [&](Superlink* sl, var_void*) { jets = *sl->jets; };
    *cutflow << [&](Superlink* sl, var_void*) {
        JetVector bjets_tmp;
        for(int i = 0; i < (int)jets.size(); i++) {
            Jet* j = jets[i];
            if(sl->tools->jetSelector().isBJet(j))  bjets_tmp.push_back(j);
            else { sjets.push_back(j); }
        }// i

        // correct the b-jets
        auto bjet_muons = sl->tools->muons_for_bjet_correction(*sl->preMuons);
        bjets = sl->tools->mu_jet_correct(bjets_tmp, bjet_muons);
    };

    int lead_bjet_flavor = -1;
    int sublead_bjet_flavor = -1;
    *cutflow << [&](Superlink* /*sl*/, var_void*) {
        if(bjets.size()>=2)
        {
            auto lead_b = bjets.at(0);
            auto sublead_b = bjets.at(1);
            lead_bjet_flavor = std::abs(lead_b.truthLabel);
            sublead_bjet_flavor = std::abs(sublead_b.truthLabel);
        }
        else
        {
            lead_bjet_flavor = -1;
            sublead_bjet_flavor = -1;
        }
    };

    *cutflow << NewVar("bjets are BB"); {
        *cutflow << HFTname("isBB");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            return (lead_bjet_flavor == 5  && sublead_bjet_flavor == 5);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjets are CC"); {
        *cutflow << HFTname("isCC");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            return (lead_bjet_flavor == 4  && sublead_bjet_flavor == 4);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjets are BC or CB"); {
        *cutflow << HFTname("isBC");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            return ((lead_bjet_flavor == 5  && sublead_bjet_flavor == 4) || (lead_bjet_flavor == 4 && sublead_bjet_flavor == 5));
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjets are LL"); {
        *cutflow << HFTname("isLL");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            bool lead_light = (lead_bjet_flavor != 5 && lead_bjet_flavor != 4);
            bool sub_light = (sublead_bjet_flavor != 5 && sublead_bjet_flavor != 4);
            return (lead_light && sub_light);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjets are BL or LB"); {
        *cutflow << HFTname("isBL");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            bool lead_light = (lead_bjet_flavor != 5 && lead_bjet_flavor != 4);
            bool sub_light = (sublead_bjet_flavor != 5 && sublead_bjet_flavor != 4);
            bool lead_b = (lead_bjet_flavor == 5);
            bool sub_b = (sublead_bjet_flavor == 5);
            bool is_bl = (lead_b && sub_light);
            bool is_lb = (lead_light && sub_b);
            return (is_bl || is_lb);
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjets are CL or LC"); {
        *cutflow << HFTname("isCL");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            bool lead_light = (lead_bjet_flavor != 5 && lead_bjet_flavor != 4);
            bool sub_light = (sublead_bjet_flavor != 5 && sublead_bjet_flavor != 4);
            bool lead_c = (lead_bjet_flavor == 4);
            bool sub_c = (sublead_bjet_flavor == 4);
            bool is_cl = (lead_c && sub_light);
            bool is_lc = (lead_light && sub_c);
            return (is_cl || is_lc);
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead jet jvt"); {
        *cutflow << HFTname("j0_jvt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>0) return jets.at(0)->jvt;
            else { return -10; }
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead sjet jvt"); {
        *cutflow << HFTname("sj0_jvt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->jvt;
            else { return -10; }
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead bjet jvt"); {
        *cutflow << HFTname("bj0_jvt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>0) return bjets.at(0).jvt;
            else { return -10; }
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet nTracks"); {
        *cutflow << HFTname("j0_nTracks");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(jets.size()>0) return jets.at(0)->nTracks;
            else { return -1; }
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet nTracks"); {
        *cutflow << HFTname("sj0_nTracks");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(sjets.size()>0) return sjets.at(0)->nTracks;
            else { return -1; }
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet nTracks"); {
        *cutflow << HFTname("bj0_nTracks");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            if(bjets.size()>0) return bjets.at(0).nTracks;
            else { return -1; }
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("jet sumTrkPt"); {
        *cutflow << HFTname("j0_sumTrkPt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>0) return jets.at(0)->sumTrkPt;
            else return -1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet sumTrkPt"); {
        *cutflow << HFTname("sj0_sumTrkPt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->sumTrkPt;
            else return -1;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet sumTrkPt"); {
        *cutflow << HFTname("bj0_sumTrkPt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>0) return bjets.at(0).sumTrkPt;
            else return -1;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet mv2c10"); {
        *cutflow << HFTname("j0_mv2c10");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>0) return jets.at(0)->mv2c10;
            else return -10;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet mv2c10"); {
        *cutflow << HFTname("sj0_mv2c10");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->mv2c10;
            else return -10;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet mv2c10"); {
        *cutflow << HFTname("bj0_mv2c10");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>0) return bjets.at(0).mv2c10;
            else return -10;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("jet dl1"); {
        *cutflow << HFTname("j0_dl1");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>0) return jets.at(0)->dl1;
            else return -10;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet dl1"); {
        *cutflow << HFTname("sj0_dl1");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->dl1;
            else return -10;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet dl1"); {
        *cutflow << HFTname("bj0_dl1");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>0) return bjets.at(0).dl1;
            else return -10;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of jets"); {
        *cutflow << HFTname("nJets");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            nn_input["nJets"] = jets.size();
            return jets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of sjets"); {
        *cutflow << HFTname("nSJets");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            nn_input["nSJets"] = sjets.size();
            return sjets.size();
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("number of bjets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
            nn_input["nBJets"] = bjets.size();
            return bjets.size();
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("lead jet pt"); {
        *cutflow << HFTname("j0_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j0_pt");
            if(jets.size()>0) val = jets.at(0)->Pt();
            nn_input["j0_pt"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead jet pt (no bdef)"); {
        *cutflow << HFTname("j0_pt_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val =  -10;
            if(jets.size()>0) val = jets.at(0)->Pt();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet pt"); {
        *cutflow << HFTname("j1_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j1_pt");
            if(jets.size()>1) val = jets.at(1)->Pt();
            nn_input["j1_pt"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet pt (no bdef)"); {
        *cutflow << HFTname("j1_pt_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(jets.size()>1) val = jets.at(1)->Pt();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead jet pt"); {
        *cutflow << HFTname("j2_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>2) return jets.at(2)->Pt();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead sjet pt"); {
        *cutflow << HFTname("sj0_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->Pt();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead sjet pt"); {
        *cutflow << HFTname("sj1_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>1) return sjets.at(1)->Pt();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead sjet pt"); {
        *cutflow << HFTname("sj2_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>2) return sjets.at(2)->Pt();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead bjet pt"); {
        *cutflow << HFTname("bj0_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj0_pt");
            if(bjets.size()>0) val = bjets.at(0).Pt();
            nn_input["bj0_pt"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead bjet pt (no bdef)"); {
        *cutflow << HFTname("bj0_pt_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>0) val = bjets.at(0).Pt();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet pt"); {
        *cutflow << HFTname("bj1_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj1_pt");
            if(bjets.size()>1) val = bjets.at(1).Pt();
            nn_input["bj1_pt"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet pt (no bdef)"); {
        *cutflow << HFTname("bj1_pt_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>1) val = bjets.at(1).Pt();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead bjet pt"); {
        *cutflow << HFTname("bj2_pt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>2) return bjets.at(2).Pt();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead jet eta"); {
        *cutflow << HFTname("j0_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j0_eta");
            if(jets.size()>0) val = jets.at(0)->Eta();
            nn_input["j0_eta"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead jet eta (no bdef)"); {
        *cutflow << HFTname("j0_eta_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(jets.size()>0) val = jets.at(0)->Eta();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet eta"); {
        *cutflow << HFTname("j1_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j1_eta");
            if(jets.size()>1) val = jets.at(1)->Eta();
            nn_input["j1_eta"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet eta (no bdef)"); {
        *cutflow << HFTname("j1_eta_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(jets.size()>1) val = jets.at(1)->Eta();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead jet eta"); {
        *cutflow << HFTname("j2_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>2)  return jets.at(2)->Eta();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead sjet eta"); {
        *cutflow << HFTname("sj0_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->Eta();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead sjet eta"); {
        *cutflow << HFTname("sj1_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>1) return sjets.at(1)->Eta();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead sjet eta"); {
        *cutflow << HFTname("sj2_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>2) return sjets.at(2)->Eta();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead bjet eta"); {
        *cutflow << HFTname("bj0_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj0_eta");
            if(bjets.size()>0) val = bjets.at(0).Eta();
            nn_input["bj0_eta"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead bjet eta (no bdef)"); {
        *cutflow << HFTname("bj0_eta_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>0) val = bjets.at(0).Eta();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet eta"); {
        *cutflow << HFTname("bj1_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj1_eta");
            if(bjets.size()>1) val = bjets.at(1).Eta();
            nn_input["bj1_eta"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet eta (no bdef)"); {
        *cutflow << HFTname("bj1_eta_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>1) val = bjets.at(1).Eta();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead bjet eta"); {
        *cutflow << HFTname("bj2_eta");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>2) return bjets.at(2).Eta();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead jet phi"); {
        *cutflow << HFTname("j0_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j0_phi");
            if(jets.size()>0) val = jets.at(0)->Phi();
            nn_input["j0_phi"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead jet phi (no bdef)"); {
        *cutflow << HFTname("j0_phi_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(jets.size()>0) val = jets.at(0)->Phi();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet phi"); {
        *cutflow << HFTname("j1_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("j1_phi");
            if(jets.size()>1) val = jets.at(1)->Phi();
            nn_input["j1_phi"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead jet phi (no bdef)"); {
        *cutflow << HFTname("j1_phi_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(jets.size()>1) val = jets.at(1)->Phi();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead jet phi"); {
        *cutflow << HFTname("j2_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(jets.size()>2) return jets.at(2)->Phi();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead sjet phi"); {
        *cutflow << HFTname("sj0_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>0) return sjets.at(0)->Phi();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead sjet phi"); {
        *cutflow << HFTname("sj1_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>1)  return sjets.at(1)->Phi();
            else return -10.;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead sjet phi"); {
        *cutflow << HFTname("sj2_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(sjets.size()>2) return sjets.at(2)->Phi();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lead bjet phi"); {
        *cutflow << HFTname("bj0_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj0_phi");
            if(bjets.size()>0) val = bjets.at(0).Phi();
            nn_input["bj0_phi"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lead bjet phi (no bdef)"); {
        *cutflow << HFTname("bj0_phi_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>0) val = bjets.at(0).Phi();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet phi"); {
        *cutflow << HFTname("bj1_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("bj1_phi");
            if(bjets.size()>1) val = bjets.at(1).Phi();
            nn_input["bj1_phi"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sub lead bjet phi (no bdef)"); {
        *cutflow << HFTname("bj1_phi_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>1) val = bjets.at(1).Phi();
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("third lead bjet phi"); {
        *cutflow << HFTname("bj2_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>2) return bjets.at(2).Phi();
            else return -10.;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading jet"); {
        *cutflow << HFTname("dphi_j0_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10.;
            if(jets.size()>0 && leptons.size()>=2) {
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
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10.;
            if(jets.size()>0) {
                out = jets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading sjet"); {
        *cutflow << HFTname("dphi_sj0_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10;
            if(sjets.size()>0 && leptons.size()>=2) {
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
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10;
            if(sjets.size()>0) {
                out = sjets.at(0)->DeltaPhi(*leptons.at(0));
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = var_means.at("dphi_bj0_ll");
            if(bjets.size()>0 && leptons.size()>=2) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = bjets.at(0).DeltaPhi(ll);
            }
            nn_input["dphi_bj0_ll"] = out;
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between dilepton system and leading bjet (no bdef)"); {
        *cutflow << HFTname("dphi_bj0_ll_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10;
            if(bjets.size()>0 && leptons.size()>=2) {
                TLorentzVector l0, l1, ll;
                l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
                l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
                ll = l0 + l1;
                out = bjets.at(0).DeltaPhi(ll);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading bjet"); {
        *cutflow << HFTname("dphi_bj0_l0");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = var_means.at("dphi_bj0_l0");
            if(bjets.size()>0) {
                out = bjets.at(0).DeltaPhi(*leptons.at(0));
            }
            nn_input["dphi_bj0_l0"] = out;
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta phi between leading lepton and leading bjet (no bdef)"); {
        *cutflow << HFTname("dphi_bj0_l0_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double out = -10;
            if(bjets.size()>0) {
                out = bjets.at(0).DeltaPhi(*leptons.at(0));
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
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double val = met.lv().Pt();
            nn_input["met"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("phi coord. of Etmiss"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double metphi = met.lv().Phi();
            nn_input["metPhi"] = metphi;
            return metphi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met TST"); {
        *cutflow << HFTname("metTST");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return met.softTerm_et; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("delta phi between dilepton system and met"); {
        *cutflow << HFTname("dphi_met_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -5;
            TLorentzVector l0, l1, ll;
            l0.SetPtEtaPhiM(leptons.at(0)->Pt(), leptons.at(0)->Eta(), leptons.at(0)->Phi(), leptons.at(0)->M());
            l1.SetPtEtaPhiM(leptons.at(1)->Pt(), leptons.at(1)->Eta(), leptons.at(1)->Phi(), leptons.at(1)->M());
            ll = l0 + l1;
            double dphi = met.lv().DeltaPhi(ll);
            nn_input["dphi_met_ll"] = dphi;
            return dphi;
        };
        *cutflow << SaveVar();
    }

    // MET terms
    *cutflow << NewVar("met_ele_et"); {
        *cutflow << HFTname("met_ele_et");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refEle_et;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_ele_phi"); {
        *cutflow << HFTname("met_ele_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refEle_phi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_ele_sumet"); {
        *cutflow << HFTname("met_ele_sumet");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refEle_sumet;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("met_jet_et"); {
        *cutflow << HFTname("met_jet_et");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refJet_et;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_jet_phi"); {
        *cutflow << HFTname("met_jet_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refJet_phi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_jet_sumet"); {
        *cutflow << HFTname("met_jet_sumet");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refJet_sumet;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_muo_et"); {
        *cutflow << HFTname("met_muo_et");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refMuo_et;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_muo_phi"); {
        *cutflow << HFTname("met_muo_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refMuo_phi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_muo_sumet"); {
        *cutflow << HFTname("met_muo_sumet");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.refMuo_sumet;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_soft_et"); {
        *cutflow << HFTname("met_soft_et");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.softTerm_et;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_soft_phi"); {
        *cutflow << HFTname("met_soft_phi");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.softTerm_phi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("met_soft_sumet"); {
        *cutflow << HFTname("met_soft_sumet");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return met.softTerm_sumet;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("mt2"); {
        *cutflow << HFTname("mt2");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            double mt2 = -10.0;
            if(leptons.size() == 2) {
                mt2 = kin::getMT2(*sl->leptons, *sl->met);
            }
            nn_input["mt2"] = mt2;
            return mt2;
        };
        *cutflow << SaveVar();
    }

    double meff;
    *cutflow << NewVar("meff : scalar sum pt of all jets, leptons, and met"); {
        *cutflow << HFTname("meff");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
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

    ///////////////////////////////////////////////
    // WWBB variables
    ///////////////////////////////////////////////

    // angles
    // dRll
    *cutflow << NewVar("delta R between two leptons"); {
        *cutflow << HFTname("dRll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -1;
            double drll = (leptons.at(0)->DeltaR(*leptons.at(1)));
            nn_input["dRll"] = drll;
            return drll;
        };
        *cutflow << SaveVar();
    }

    // M_bb
    *cutflow << NewVar("invariant mass of di-bjet system"); {
        *cutflow << HFTname("mbb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double mbb = var_means.at("mbb");
            if(bjets.size()>=2) {
                mbb = (bjets.at(0) + bjets.at(1)).M();
            }
            nn_input["mbb"] = mbb;
            return mbb;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("invariant mass of di-bjet system (no bdef)"); {
        *cutflow << HFTname("mbb_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double mbb = -10;
            if(bjets.size()>=2) {
                mbb = (bjets.at(0) + bjets.at(1)).M();
            }
            return mbb;
        };
        *cutflow << SaveVar();
    }

    // dRbb
    *cutflow << NewVar("delta R between two leading b-jets"); {
        *cutflow << HFTname("dRbb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2) {
                return (bjets.at(0).DeltaR(bjets.at(1)));
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }


    // dR_ll_bb
    *cutflow << NewVar("delta R between dilepton system and di-bjet system"); {
        *cutflow << HFTname("dR_ll_bb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2 && leptons.size()>=2) {
                TLorentzVector l0 = (*leptons.at(0));
                TLorentzVector l1 = (*leptons.at(1));
                TLorentzVector b0 = (bjets.at(0));
                TLorentzVector b1 = (bjets.at(1));

                return ( (l0 + l1).DeltaR( (b0 + b1) ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // dphi bb  ll
    *cutflow << NewVar("delta phi between bb and ll systems"); {
        *cutflow << HFTname("dphi_ll_bb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2 && leptons.size()>=2) {
                return ( (bjets.at(0) + bjets.at(1)).DeltaPhi( (*leptons.at(0) + *leptons.at(1)) ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // dphi WW bb
    *cutflow << NewVar("delta phi between WW and bb systems"); {
        *cutflow << HFTname("dphi_WW_bb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2 && leptons.size()>=2) {
                return ( (met.lv() + *leptons.at(0) + *leptons.at(1)).DeltaPhi( (bjets.at(0) + bjets.at(1)) ) );
            }
            return -10.;
        };
        *cutflow << SaveVar();
    }

    // dphi_met_ll
    *cutflow << NewVar("delta phi between MET and dilepton system"); {
        *cutflow << HFTname("dphi_met_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -5;
            return ( (met.lv().DeltaPhi( (*leptons.at(0) + *leptons.at(1)) )) );
        };
        *cutflow << SaveVar();
    }

    // mass_met_ll
    *cutflow << NewVar("mass of met and dilepton system"); {
        *cutflow << HFTname("mass_met_ll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -1;
            return ( (met.lv() + *leptons.at(0) + *leptons.at(1)).M() );
        };
        *cutflow << SaveVar();
    }

    // met_pTll
    *cutflow << NewVar("pT of met + dilepton ssytem"); {
        *cutflow << HFTname("met_pTll");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -1;
            double val = (met.lv() + *leptons.at(0) + *leptons.at(1)).Pt();
            nn_input["met_pTll"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }

    // HT2
    *cutflow << NewVar("HT2"); {
        *cutflow << HFTname("HT2");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float out = var_means.at("HT2");
            if(bjets.size()>=2 && leptons.size()>=2) {
                double HT2 = ( (bjets.at(0) + bjets.at(1)).Pt() +
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );
                out = HT2;
            }
            nn_input["HT2"] = out;
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT2 (no bdef)"); {
        *cutflow << HFTname("HT2_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float out = -10;
            if(bjets.size()>=2 && leptons.size()>=2) {
                double HT2 = ( (bjets.at(0) + bjets.at(1)).Pt() +
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );
                out = HT2;
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    // HT2Ratio
    *cutflow << NewVar("HT2Ratio"); {
        *cutflow << HFTname("HT2Ratio");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float out = var_means.at("HT2Ratio");
            if(bjets.size()>=2 && leptons.size()>=2) {
                double num = ( (bjets.at(0) + bjets.at(1)).Pt() +
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((bjets.at(0)).Pt());
                den += (bjets.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                out = (num/den);
            }
            nn_input["HT2Ratio"] = out;
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("HT2Ratio (no bdef)"); {
        *cutflow << HFTname("HT2Ratio_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float out = -10;
            if(bjets.size()>=2 && leptons.size()>=2) {
                double num = ( (bjets.at(0) + bjets.at(1)).Pt() +
                    (*leptons.at(0) + *leptons.at(1) + met.lv()).Pt() );

                double den = ((bjets.at(0)).Pt());
                den += (bjets.at(1)).Pt();
                den += (*leptons.at(0)).Pt();
                den += (*leptons.at(1)).Pt();
                den += met.lv().Pt();
                out = (num/den);
            }
            return out;
        };
        *cutflow << SaveVar();
    }


    // MT_HWW
    *cutflow << NewVar("MT_HWW"); {
        *cutflow << HFTname("MT_HWW");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(leptons.size()<2) return -1;
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
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2 && leptons.size()>=2) {
                TLorentzVector vis = (*leptons.at(0) + *leptons.at(1)
                            + bjets.at(0) + bjets.at(1));
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
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            if(bjets.size()>=2 && leptons.size()>=2) {
                TLorentzVector bjet_system = (bjets.at(0) + bjets.at(1));
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

    // mt1 (l0, l1), (b0, b1)
    *cutflow << NewVar("mt2_llbb"); {
        *cutflow << HFTname("mt2_llbb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("mt2_llbb");
            if(bjets.size()>=2 && leptons.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( bjets.at(0) + bjets.at(1) );
                val = kin::getMT2(v0,v1,met);
            }
            nn_input["mt2_llbb"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mt2_llbb (no bdef)"); {
        *cutflow << HFTname("mt2_llbb_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>=2 && leptons.size()>=2) {
                const TLorentzVector v0 = ( *leptons.at(0) + *leptons.at(1) );
                const TLorentzVector v1 = ( bjets.at(0) + bjets.at(1) );
                val = kin::getMT2(v0,v1,met);
            }
            return val;
        };
        *cutflow << SaveVar();
    }


    // mt2_bb
    *cutflow << NewVar("mt2_bb"); {
        *cutflow << HFTname("mt2_bb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = var_means.at("mt2_bb");
            if(bjets.size()>=2) {
                const TLorentzVector b0 = (bjets.at(0));
                const TLorentzVector b1 = (bjets.at(1));
                val = kin::getMT2(b0,b1,met);
            }
            nn_input["mt2_bb"] = val;
            return val;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("mt2_bb (no bdef)"); {
        *cutflow << HFTname("mt2_bb_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            float val = -10;
            if(bjets.size()>=2) {
                const TLorentzVector b0 = (bjets.at(0));
                const TLorentzVector b1 = (bjets.at(1));
                val = kin::getMT2(b0,b1,met);
            }
            return val;
        };
        *cutflow << SaveVar();
    }

    // dphi_bb
    *cutflow << NewVar("dphi_bb"); {
        *cutflow << HFTname("dphi_bb");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double dphi_bb = var_means.at("dphi_bb");
            if(bjets.size()>=2) {
                dphi_bb = (bjets.at(0).DeltaPhi(bjets.at(1)));
            }
            nn_input["dphi_bb"] = dphi_bb;
            return dphi_bb;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("dphi_bb (no bdef)"); {
        *cutflow << HFTname("dphi_bb_nobdef");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            double dphi_bb = -10;
            if(bjets.size()>=2) {
                dphi_bb = (bjets.at(0).DeltaPhi(bjets.at(1)));
            }
            return dphi_bb;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mT_W"); {
        *cutflow << HFTname("mT_W");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            TLorentzVector l0 = (*leptons.at(0));
            l0.SetPz(0.0);
            return ( l0 + met.lv() ).Mt();
        };
        *cutflow << SaveVar();
    }

    // NN STUFF
    // NN VARS
    double nn_output_hh = -999;
    double nn_output_tt = -999;
    double nn_output_wt = -999;
    double nn_output_zjets = -999;
    std::map< std::string, std::map< std::string, double >> inputs;

    double nn_output_hh_update = -999;
    double nn_output_top_update = -999;
    double nn_output_zsf_update = -999;
    double nn_output_ztt_update = -999;
    std::map< std::string, std::map< std::string, double >> inputs_update;

    *cutflow << [&](Superlink* /*sl*/, var_void*) {

        inputs["InputLayer"] = nn_input;
        auto output_scores = nn_graph.compute(inputs);
        nn_output_hh = output_scores.at("out_0_hh");
        nn_output_tt = output_scores.at("out_1_tt");
        nn_output_wt = output_scores.at("out_2_wt");
        nn_output_zjets = output_scores.at("out_3_zjets");

        // update R21 network
        inputs_update["InputLayer"] = nn_input;
        auto output_scores_update = nn_graph_update.compute(inputs_update);
        nn_output_hh_update = output_scores_update.at("out_0_hh");
        nn_output_top_update = output_scores_update.at("out_1_top");
        nn_output_zsf_update = output_scores_update.at("out_2_zsf");
        nn_output_ztt_update = output_scores_update.at("out_3_ztt");
    };

    *cutflow << NewVar("NN_p_hh_R20"); {
        *cutflow << HFTname("NN_p_hh_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_hh;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_tt_R20"); {
        *cutflow << HFTname("NN_p_tt_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_tt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_wt_R20"); {
        *cutflow << HFTname("NN_p_wt_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_wt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_zjets_R20"); {
        *cutflow << HFTname("NN_p_zjets_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_zjets;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_hh"); {
        *cutflow << HFTname("NN_p_hh");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_hh_update;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_top"); {
        *cutflow << HFTname("NN_p_top");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_top_update;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_zsf"); {
        *cutflow << HFTname("NN_p_zsf");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_zsf_update;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_p_ztt"); {
        *cutflow << HFTname("NN_p_ztt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return nn_output_ztt_update;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("NN_d_hh_R20"); {
        *cutflow << HFTname("NN_d_hh_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_hh / (nn_output_tt + nn_output_wt + nn_output_zjets) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_tt_R20"); {
        *cutflow << HFTname("NN_d_tt_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_tt / (nn_output_hh + nn_output_wt + nn_output_zjets) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_wt_R20"); {
        *cutflow << HFTname("NN_d_wt_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_wt / (nn_output_hh + nn_output_tt + nn_output_zjets) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_zjets_R20"); {
        *cutflow << HFTname("NN_d_zjets_R20");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_zjets / (nn_output_hh + nn_output_tt + nn_output_wt) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_hh"); {
        *cutflow << HFTname("NN_d_hh");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_hh_update / (nn_output_top_update + nn_output_zsf_update + nn_output_ztt_update) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_top"); {
        *cutflow << HFTname("NN_d_top");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_top_update / (nn_output_hh_update + nn_output_zsf_update + nn_output_ztt_update) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_zsf"); {
        *cutflow << HFTname("NN_d_zsf");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_zsf_update / (nn_output_hh_update + nn_output_top_update + nn_output_ztt_update) );
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("NN_d_ztt"); {
        *cutflow << HFTname("NN_d_ztt");
        *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
            return log( nn_output_ztt_update / (nn_output_hh_update + nn_output_top_update + nn_output_zsf_update) );
        };
        *cutflow << SaveVar();
    }


    // clear the wectors
    *cutflow << [&](Superlink* /* sl */, var_void*) { leptons.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { electrons.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { muons.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { jets.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* /* sl */, var_void*) { met.clear(); };


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

//    *cutflow << NewSystematic("electron ID efficiency"); {
//        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_UP, SupersysWeight::EL_EFF_ID_TOTAL_Uncorr_DN);
//        *cutflow << TreeName("EL_EFF_ID");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("electron ISO efficiency"); {
//        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_UP, SupersysWeight::EL_EFF_Iso_TOTAL_Uncorr_DN);
//        *cutflow << TreeName("EL_EFF_ISO");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("electron RECO efficiency"); {
//        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_UP, SupersysWeight::EL_EFF_Reco_TOTAL_Uncorr_DN);
//        *cutflow << TreeName("EL_EFF_RECO");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("electron Trigger efficiency"); {
//        *cutflow << WeightSystematic(SupersysWeight::EL_EFF_Trigger_TOTAL_UP, SupersysWeight::EL_EFF_Trigger_TOTAL_DN);
//        *cutflow << TreeName("EL_EFF_TRIG");
//        *cutflow << SaveSystematic();
//    }

    // flavor tagging
    *cutflow << NewSystematic("FTAG EFF B"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_B_UP, SupersysWeight::FT_EFF_B_DN);
        *cutflow << TreeName("FT_EFF_B");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("FTAG EFF C"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_C_UP, SupersysWeight::FT_EFF_C_DN);
        *cutflow << TreeName("FT_EFF_C");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("FTAG EFF LIGHT"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_LT_UP, SupersysWeight::FT_EFF_LT_DN);
        *cutflow << TreeName("FT_EFF_LT");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("FTAG EFF EXTRAP"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAP_UP, SupersysWeight::FT_EFF_EXTRAP_DN);
        *cutflow << TreeName("FT_EFF_EXTRAP");
        *cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("FTAG EFF EXTRAPC"); {
        *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAPC_UP, SupersysWeight::FT_EFF_EXTRAPC_DN);
        *cutflow << TreeName("FT_EFF_EXTRAPC");
        *cutflow << SaveSystematic();
    }

    // JVT
 //   *cutflow << NewSystematic("JVT EFF"); {
 //       *cutflow << WeightSystematic(SupersysWeight::JET_JVTEff_UP, SupersysWeight::JET_JVTEff_DN);
 //       *cutflow << TreeName("JVTEff");
 //       *cutflow << SaveSystematic();
 //   }

//    // muon
//    *cutflow << NewSystematic("Muon Bad Muon Stat"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_BADMUON_STAT_UP, SupersysWeight::MUON_EFF_BADMUON_STAT_DN);
//        *cutflow << TreeName("BadMuonStat");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Bad Muon Syst"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_BADMUON_SYS_UP, SupersysWeight::MUON_EFF_BADMUON_SYS_DN);
//        *cutflow << TreeName("BadMuonSys");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff ISO Stat"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_STAT_UP, SupersysWeight::MUON_EFF_ISO_STAT_DN);
//        *cutflow << TreeName("MU_EFF_ISO_STAT");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff ISO Syst"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_SYS_UP, SupersysWeight::MUON_EFF_ISO_SYS_DN);
//        *cutflow << TreeName("MU_EFF_ISO_SYS");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff RECO Stat"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_RECO_STAT_UP, SupersysWeight::MUON_EFF_RECO_STAT_DN);
//        *cutflow << TreeName("MU_EFF_RECO_STAT");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff RECO Syst"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_RECO_SYS_UP, SupersysWeight::MUON_EFF_RECO_SYS_DN);
//        *cutflow << TreeName("MU_EFF_RECO_SYS");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff TTVA STAT"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_TTVA_STAT_UP, SupersysWeight::MUON_EFF_TTVA_STAT_DN);
//        *cutflow << TreeName("MU_EFF_TTVA_STAT");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff TTVA SYS"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_TTVA_SYS_UP, SupersysWeight::MUON_EFF_TTVA_SYS_DN);
//        *cutflow << TreeName("MU_EFF_TTVA_SYS");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff TRIG STAT"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_TrigStat_UP, SupersysWeight::MUON_EFF_TrigStat_DN);
//        *cutflow << TreeName("MU_EFF_TRIG_STAT");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("Muon Eff TRIG SYS"); {
//        *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_TrigSys_UP, SupersysWeight::MUON_EFF_TrigSys_DN);
//        *cutflow << TreeName("MU_EFF_TRIG_SYS");
//        *cutflow << SaveSystematic();
//    }

    // pileup
//    *cutflow << NewSystematic("Pileup"); {
//        *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DN);
//        *cutflow << TreeName("PILEUP");
//        *cutflow << SaveSystematic();
//    }

//    *cutflow << NewSystematic("EG_RESOLUTION_ALL_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
//    	*cutflow << TreeName("EG_RESOLUTION_ALL_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("EG_RESOLUTION_ALL_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
//    	*cutflow << TreeName("EG_RESOLUTION_ALL_UP");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("EG_SCALE_ALL_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
//    	*cutflow << TreeName("EG_SCALE_ALL_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("EG_SCALE_ALL_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
//    	*cutflow << TreeName("EG_SCALE_ALL_UP");
//    	*cutflow << SaveSystematic();
//    }
    *cutflow << NewSystematic("JET_JER_DataVsMC_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_DataVsMC_UP);
    	*cutflow << TreeName("JET_JER_DataVsMC_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_DataVsMC_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_DataVsMC_DN);
    	*cutflow << TreeName("JET_JER_DataVsMC_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_1_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_1_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_1_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_1_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_1_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_1_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_2_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_2_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_2_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_2_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_2_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_2_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_3_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_3_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_3_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_3_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_3_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_3_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_4_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_4_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_4_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_4_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_4_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_4_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_5_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_5_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_5_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_5_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_5_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_5_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_6_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_6_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_6_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_6_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_6_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_6_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_7rest_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_7rest_UP);
    	*cutflow << TreeName("JET_JER_EffectiveNP_7rest_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_JER_EffectiveNP_7rest_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_JER_EffectiveNP_7rest_DN);
    	*cutflow << TreeName("JET_JER_EffectiveNP_7rest_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_1_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
    	*cutflow << TreeName("JET_GroupedNP_1_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_1_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
    	*cutflow << TreeName("JET_GroupedNP_1_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_2_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_2_UP);
    	*cutflow << TreeName("JET_GroupedNP_2_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_2_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_2_DN);
    	*cutflow << TreeName("JET_GroupedNP_2_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_3_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_3_UP);
    	*cutflow << TreeName("JET_GroupedNP_3_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_GroupedNP_3_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_GroupedNP_3_DN);
    	*cutflow << TreeName("JET_GroupedNP_3_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_EtaIntercalibration_UP (UP)"); {
    	*cutflow << EventSystematic(NtSys::JET_EtaIntercalibration_UP);
    	*cutflow << TreeName("JET_EtaIntercalibration_UP");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("JET_EtaIntercalibration_DN (DN)"); {
    	*cutflow << EventSystematic(NtSys::JET_EtaIntercalibration_DN);
    	*cutflow << TreeName("JET_EtaIntercalibration_DN");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET_SoftTrk_ResoPara"); {
    	*cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
    	*cutflow << TreeName("MET_SoftTrk_ResoPara");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET_SoftTrk_ResoPerp"); {
    	*cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
    	*cutflow << TreeName("MET_SoftTrk_ResoPerp");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET_SoftTrk_ScaleDown"); {
    	*cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
    	*cutflow << TreeName("MET_SoftTrk_ScaleDown");
    	*cutflow << SaveSystematic();
    }
    *cutflow << NewSystematic("MET_SoftTrk_ScaleUp"); {
    	*cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
    	*cutflow << TreeName("MET_SoftTrk_ScaleUp");
    	*cutflow << SaveSystematic();
    }
//    *cutflow << NewSystematic("MUON_MS_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_MS_DN);
//    	*cutflow << TreeName("MUON_MS_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_MS_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_MS_UP);
//    	*cutflow << TreeName("MUON_MS_UP");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_ID_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_ID_DN);
//    	*cutflow << TreeName("MUON_ID_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_ID_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_ID_UP);
//    	*cutflow << TreeName("MUON_ID_UP");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SAGITTA_RESBIAS_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SAGITTA_RESBIAS_DN);
//    	*cutflow << TreeName("MUON_SAGITTA_RESBIAS_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SAGITTA_RESBIAS_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SAGITTA_RESBIAS_UP);
//    	*cutflow << TreeName("MUON_SAGITTA_RESBIAS_UP");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SAGITTA_RHO_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SAGITTA_RHO_DN);
//    	*cutflow << TreeName("MUON_SAGITTA_RHO_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SAGITTA_RHO_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SAGITTA_RHO_UP);
//    	*cutflow << TreeName("MUON_SAGITTA_RHO_UP");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SCALE_DN (DN)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SCALE_DN);
//    	*cutflow << TreeName("MUON_SCALE_DN");
//    	*cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("MUON_SCALE_UP (UP)"); {
//    	*cutflow << EventSystematic(NtSys::MUON_SCALE_UP);
//    	*cutflow << TreeName("MUON_SCALE_UP");
//    	*cutflow << SaveSystematic();
//    }

//
//    ////////////////////////////////////
//    // shape systematics
//    ////////////////////////////////////
//
//    // egamma
//    *cutflow << NewSystematic("shift in e-gamma resolution (UP)"); {
//        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
//        *cutflow << TreeName("EG_RESOLUTION_ALL_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("shift in e-gamma resolution (DOWN)"); {
//        *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
//        *cutflow << TreeName("EG_RESOLUTION_ALL_DN");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("shift in e-gamma scale (UP)"); {
//        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
//        *cutflow << TreeName("EG_SCALE_ALL_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("shift in e-gamma scale (DOWN)"); {
//        *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
//        *cutflow << TreeName("EG_SCALE_ALL_DN");
//        *cutflow << SaveSystematic();
//    }
//    // muon
//    *cutflow << NewSystematic("muon ID (UP)"); {
//        *cutflow << EventSystematic(NtSys::MUON_ID_UP);
//        *cutflow << TreeName("MUON_ID_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("muon ID (DOWN)"); {
//        *cutflow << EventSystematic(NtSys::MUON_ID_DN);
//        *cutflow << TreeName("MUON_ID_DN");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("muon MS (UP)"); {
//        *cutflow << EventSystematic(NtSys::MUON_MS_UP);
//        *cutflow << TreeName("MUON_MS_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("muon MS (DOWN)"); {
//        *cutflow << EventSystematic(NtSys::MUON_MS_DN);
//        *cutflow << TreeName("MUON_MS_DN");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("muon scale shift (UP)"); {
//        *cutflow << EventSystematic(NtSys::MUON_SCALE_UP);
//        *cutflow << TreeName("MUON_SCALE_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("muon scale shift (DN)"); {
//        *cutflow << EventSystematic(NtSys::MUON_SCALE_DN);
//        *cutflow << TreeName("MUON_SCALE_DN");
//        *cutflow << SaveSystematic();
//    }
//
//    // jet
//    *cutflow << NewSystematic("JER"); {
//        *cutflow << EventSystematic(NtSys::JER);
//        *cutflow << TreeName("JER");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 1 (up)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
//        *cutflow << TreeName("JET_GroupedNP_1_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 1 (down)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
//        *cutflow << TreeName("JET_GroupedNP_1_DN");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 2 (up)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_UP);
//        *cutflow << TreeName("JET_GroupedNP_2_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 2 (down)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_2_DN);
//        *cutflow << TreeName("JET_GroupedNP_2_DN");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 3 (up)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_UP);
//        *cutflow << TreeName("JET_GroupedNP_3_UP");
//        *cutflow << SaveSystematic();
//    }
//    *cutflow << NewSystematic("JES NP set 3 (down)"); {
//        *cutflow << EventSystematic(NtSys::JET_GroupedNP_3_DN);
//        *cutflow << TreeName("JET_GroupedNP_3_DN");
//        *cutflow << SaveSystematic();
//    }
//    //#warning only setting MET systematics
//    // met
//    //*cutflow << NewSystematic("MET TST Soft-Term resolution (parallel)"); {
//    //    *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
//    //    *cutflow << TreeName("MET_SoftTrk_ResoPara");
//    //    *cutflow << SaveSystematic();
//    //}
//    //*cutflow << NewSystematic("MET TST Soft-Term resolution (perpendicular)"); {
//    //    *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
//    //    *cutflow << TreeName("MET_SoftTrk_ResoPerp");
//    //    *cutflow << SaveSystematic();
//    //}
//    //*cutflow << NewSystematic("MET TST Soft-Term shift in scale (UP)"); {
//    //    *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
//    //    *cutflow << TreeName("MET_SoftTrk_ScaleUp");
//    //    *cutflow << SaveSystematic();
//    //}
//    //*cutflow << NewSystematic("MET TST Soft-Term shift in scale (DOWN)"); {
//    //    *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
//    //    *cutflow << TreeName("MET_SoftTrk_ScaleDown");
//    //    *cutflow << SaveSystematic();
//    //}


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //
    // Superflow methods [END]
    //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

//    delete pu_profile;

    // initialize the cutflow and start the event loop
    chain->Process(cutflow, options.input.c_str(), options.n_events_to_process);
    delete cutflow;
    delete chain;
    cout << "La Fin." << endl;
    exit(0);

}
