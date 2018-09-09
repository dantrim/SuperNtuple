#include "SuperNtuple/NNMaker.h"

// SusyNtuple
#include "SusyNtuple/KinematicTools.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/string_utils.h"
using namespace Susy; // everything in SusyNtuple is in this namespace

//ROOT

// std/stl
#include <iomanip> // setw
#include <iostream>
#include <string>
#include <sstream> // stringstream, ostringstream
#include <fstream>
#include <math.h>
using namespace std;

#define ADD_FEATURE_F( name ) \
    do { \
        chain()->SetBranchAddress(#name, &feature_map_f[#name]); \
    } while (0);

#define ADD_FEATURE_I( name ) \
    do { \
        chain()->SetBranchAddress(#name, &feature_map_i[#name]); \
    } while (0);


//////////////////////////////////////////////////////////////////////////////
NNMaker::NNMaker() :
    m_dbg(0),
    m_suffix(""),
    m_input_name(""),
    m_input_chain(nullptr),
    m_lwt_graph(nullptr),
    m_output_name("")
{
    feature_map_f.clear();
    feature_map_i.clear();
}
//////////////////////////////////////////////////////////////////////////////
bool NNMaker::load_nn(std::string filename)
{
    cout << "NNMaker::load_nn    Loading NN from file: " << filename << endl;
    std::ifstream input_nn_file(filename);
    if(!input_nn_file) {
        cout << "NNMaker::load_nn Provided NN file (=" << filename << ") does not exist" << endl;
        return false;
    }
    string output_layer_name = "OutputLayer";
    auto config = lwt::parse_json_graph(input_nn_file);
    m_lwt_graph = new lwt::LightweightGraph(config, output_layer_name);

    cout << "NNMaker::load_nn    Finished loading NN..." << endl;

    return true;
}
//////////////////////////////////////////////////////////////////////////////
bool NNMaker::set_input_name(std::string input_name)
{
    // expect only ROOT files
    if(!Susy::utils::endswith(input_name, ".root")) {
        cout << "NNMaker::set_input_name    ERROR Input does not end with \".root\", expecing a single ROOT file as input" << endl;
        return false;
    }
    m_input_name = input_name;


    return true;
}
//////////////////////////////////////////////////////////////////////////////
void NNMaker::reset()
{
    m_p_hh = -99;
    m_p_tt = -99;
    m_p_wt = -99;
    m_p_z = -99;

    m_d_hh = -99;
    m_d_tt = -99;
    m_d_wt = -99;
    m_d_z = -99;

    m_pr_hh_tt = -99;
    m_pr_hh_wt = -99;
    m_pr_hh_z = -99;
    m_pr_hh_tt_wt = -99;
    m_pr_hh_tt_wt_z = -99;
}
//////////////////////////////////////////////////////////////////////////////
void NNMaker::setup_nn_branches()
{
    auto splits = Susy::utils::tokenizeString(input_filename(), '/'); 
    //cout << "NNMaker::setup_nn_branches    Found " << splits.size() << " splits: ";
    //for(auto & split : splits) {
    //    cout << " " << split;
    //}
    //cout << endl;

    std::string outname = splits.at(splits.size()-1);
    splits = Susy::utils::tokenizeString(outname, '.');
    outname = splits.at(0);
    if(suffix() != "") outname += "_" + suffix();
    outname += ".root";
    m_output_name = outname;
    cout << "NNMaker::setup_nn_branches    Output filename from input name : " << output_filename() << endl;

    m_output_file = new TFile(output_filename().c_str(), "RECREATE");
//    m_output_tree = new TTree("superNt", "superNt");
    m_output_tree = chain()->CloneTree(0);

    cout << "NNMaker::setup_nn_branches    Copying branches of old (input) TTree" << endl;
    //chain()->SetBranchStatus("*", 1);
    m_output_tree->CopyEntries(0);
    //m_output_tree->CopyEntries(chain());
    cout << "NNMaker::setup_nn_branches    Done copying..." << endl;

    b_p_hh = m_output_tree->Branch("bdef_nn_p_hh", &m_p_hh);
    b_p_tt = m_output_tree->Branch("bdef_nn_p_tt", &m_p_tt);
    b_p_wt = m_output_tree->Branch("bdef_nn_p_wt", &m_p_wt);
    b_p_z = m_output_tree->Branch("bdef_nn_p_z", &m_p_z);

    b_d_hh = m_output_tree->Branch("bdef_nn_d_hh", &m_d_hh);
    b_d_tt = m_output_tree->Branch("bdef_nn_d_tt", &m_d_tt);
    b_d_wt = m_output_tree->Branch("bdef_nn_d_wt", &m_d_wt);
    b_d_z = m_output_tree->Branch("bdef_nn_d_z", &m_d_z);

    b_lr_hh_tt = m_output_tree->Branch("bdef_nn_lr_hh_tt", &m_lr_hh_tt);
    b_lr_hh_wt = m_output_tree->Branch("bdef_nn_lr_hh_wt", &m_lr_hh_wt);
    b_lr_hh_z = m_output_tree->Branch("bdef_nn_lr_hh_z", &m_lr_hh_z);

    b_pr_hh_tt = m_output_tree->Branch("bdef_nn_pr_hh_tt", &m_pr_hh_tt);
    b_pr_hh_wt = m_output_tree->Branch("bdef_nn_pr_hh_wt", &m_pr_hh_wt);
    b_pr_hh_z = m_output_tree->Branch("bdef_nn_pr_hh_z", &m_pr_hh_z);
    b_pr_hh_tt_wt_z = m_output_tree->Branch("bdef_nn_pr_hh_tt_wt_z", &m_pr_hh_tt_wt_z);


    ADD_FEATURE_F(met)
    ADD_FEATURE_F(metPhi)
    ADD_FEATURE_F(mll)
    ADD_FEATURE_F(dRll)
    ADD_FEATURE_F(pTll)
    ADD_FEATURE_F(dphi_ll)
    ADD_FEATURE_F(dphi_bb)
    ADD_FEATURE_F(dphi_met_ll)
    ADD_FEATURE_F(met_pTll)
    ADD_FEATURE_I(nJets)
    ADD_FEATURE_I(nSJets)
    ADD_FEATURE_I(nBJets)
    ADD_FEATURE_I(isEE)
    ADD_FEATURE_I(isMM)
    ADD_FEATURE_I(isSF)
    ADD_FEATURE_I(isDF)
    ADD_FEATURE_F(HT2)
    ADD_FEATURE_F(HT2Ratio)
    ADD_FEATURE_F(l0_pt)
    ADD_FEATURE_F(l1_pt)
    ADD_FEATURE_F(l0_phi)
    ADD_FEATURE_F(l1_phi)
    ADD_FEATURE_F(j0_pt)
    ADD_FEATURE_F(j1_pt)
    ADD_FEATURE_F(j0_phi)
    ADD_FEATURE_F(j1_phi)
    ADD_FEATURE_F(bj0_pt)
    ADD_FEATURE_F(bj1_pt)
    ADD_FEATURE_F(bj0_phi)
    ADD_FEATURE_F(bj1_phi)

}
//////////////////////////////////////////////////////////////////////////////
float NNMaker::feature(string fname)
{
    if(feature_map_f.count(fname)) {
        return feature_map_f.at(fname);
    }
    else if(feature_map_i.count(fname)) {
        return static_cast<float>(feature_map_i.at(fname));
    }
    cout << "NNMaker::feature    WARNING Did not find requested feature " << fname << " in either feature map" << endl;
    return -999;
}
//////////////////////////////////////////////////////////////////////////////
std::map<std::string, double> NNMaker::lwt_map()
{
    for(auto m : feature_map_f) {
        m_lwt_map[m.first] = static_cast<double>(m.second);
    }
    for(auto m : feature_map_i) {
        m_lwt_map[m.first] = static_cast<double>(m.second);
    }
    return m_lwt_map;
}
//////////////////////////////////////////////////////////////////////////////
bool NNMaker::process(long int n_to_process)
{
    // before the event loop starts we should setup the connections and outputs
    setup_nn_branches();
    std::map< std::string, std::map< std::string, double >> lwt_inputs;

    // start looping
    long int  n_total = chain()->GetEntries();
    if(n_to_process < 0) n_to_process = n_total;
    for(long int entry = 0; entry < n_to_process; entry++) {

        // clear out the observables we are writing new for each new event
        reset();

        chain()->GetEntry(entry);
        if(dbg() || entry%1000==0) {
            cout << "NNMaker::process     **** Processing entry " << setw(6) << entry << " ****" << endl;
        }

        lwt_inputs["InputLayer"] = lwt_map(); 
        auto scores = graph()->compute(lwt_inputs);

        m_p_hh = scores.at("out_0_hh");
        m_p_tt = scores.at("out_1_tt");
        m_p_wt = scores.at("out_2_wt");
        m_p_z = scores.at("out_3_zjets");

        m_d_hh = log( m_p_hh / (m_p_tt + m_p_wt + m_p_z) );
        m_d_tt = log( m_p_tt / (m_p_hh + m_p_wt + m_p_z) );
        m_d_wt = log( m_p_wt / (m_p_hh + m_p_tt + m_p_z) );
        m_d_z = log( m_p_z / (m_p_hh + m_p_tt + m_p_wt) );

        m_lr_hh_tt = log( m_p_hh / m_p_tt );
        m_lr_hh_wt = log( m_p_hh / m_p_wt );
        m_lr_hh_z = log( m_p_hh / m_p_z );

        m_pr_hh_tt = m_p_hh / m_p_tt;
        if(std::isnan(m_pr_hh_tt)) m_pr_hh_tt = -99;
        m_pr_hh_wt = m_p_hh / m_p_wt;
        if(std::isnan(m_pr_hh_wt)) m_pr_hh_wt = -99;
        m_pr_hh_z = m_p_hh / m_p_z;
        if(std::isnan(m_pr_hh_z)) m_pr_hh_z = -99;
        m_pr_hh_tt_wt = m_p_hh / (m_p_tt + m_p_wt);
        if(std::isnan(m_pr_hh_tt_wt)) m_pr_hh_tt_wt = -99;
        m_pr_hh_tt_wt_z = m_p_hh / (m_p_tt + m_p_wt + m_p_z);
        if(std::isnan(m_pr_hh_tt_wt_z)) m_pr_hh_tt_wt_z = -99;

        // at the end of each event, push the output TTree
        m_output_tree->Fill();
    } // entry
    
    // all done
    m_output_tree->Write();

    return true;
}
//////////////////////////////////////////////////////////////////////////////
//Bool_t NNMaker::Process(Long64_t entry)
//{
//    static Long64_t chain_entry = -1;
//    chain_entry++;
//    if(chain_entry == 0) {
//        setup_nn_branches();
//    }
//    chain()->GetEntry(chain_entry);
//
//    if(dbg() || chain_entry % 1000 == 0) {
//        cout << "NNMaker::Process    **** Processing entry " << setw(6) << chain_entry << " **** " << endl;
//    }
//
//    m_output_tree->Fill();
//    return kTRUE;
//}
//////////////////////////////////////////////////////////////////////////////
//void NNMaker::Terminate()
//{
//    if(m_lwt_graph) delete m_lwt_graph;
//
//    if(m_output_file) {
//     //   m_output_file->cd();
//        m_output_tree->Write();
//        m_output_file->Write();
//        //m_output_file->Close();
//    }
//
//    return;
//}
//////////////////////////////////////////////////////////////////////////////
