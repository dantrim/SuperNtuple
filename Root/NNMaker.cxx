#include "SuperNtuple/NNMaker.h"

// SusyNtuple
#include "SusyNtuple/KinematicTools.h"
#include "SusyNtuple/SusyDefs.h"
using namespace Susy; // everything in SusyNtuple is in this namespace

//ROOT

// std/stl
#include <iomanip> // setw
#include <iostream>
#include <string>
#include <sstream> // stringstream, ostringstream
#include <fstream>
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
    m_input_chain(nullptr),
    m_lwt_graph(nullptr)
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
void NNMaker::setup_nn_branches()
{
    string output_name = "my_test";
    if(suffix()!="") output_name += "_" + suffix();
    output_name += ".root";

    m_output_file = new TFile(output_name.c_str(), "RECREATE");
    m_output_tree = new TTree("superNt", "superNt");

    cout << "NNMaker::setup_nn_branches    Copying branches of old (input) TTree" << endl;
    //chain()->SetBranchStatus("*", 1);
    m_output_tree->CopyEntries(0);
    m_output_tree->CopyEntries(chain());
    cout << "NNMaker::setup_nn_branches    Done copying..." << endl;

    b_p_hh = m_output_tree->Branch("nn_p_hh", &m_p_hh);
    b_p_tt = m_output_tree->Branch("nn_p_tt", &m_p_tt);
    b_p_wt = m_output_tree->Branch("nn_p_wt", &m_p_wt);
    b_p_z = m_output_tree->Branch("nn_p_z", &m_p_z);


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
void NNMaker::Init(TTree* tree)
{
    cout << "NNMaker::Init" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void NNMaker::Begin(TTree* /*tree*/)
{
    // call base class' Begin method
    if(dbg()) cout << "NNMaker::Begin" << endl;

    return;
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
Bool_t NNMaker::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    if(chain_entry == 0) {
        setup_nn_branches();
    }
    chain()->GetEntry(chain_entry);

    if(dbg() || chain_entry % 1000 == 0) {
        cout << "NNMaker::Process    **** Processing entry " << setw(6) << chain_entry << " **** " << endl;
    }

    m_output_tree->Fill();
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void NNMaker::Terminate()
{
    if(m_lwt_graph) delete m_lwt_graph;

    if(m_output_file) {
     //   m_output_file->cd();
        m_output_tree->Write();
        m_output_file->Write();
        //m_output_file->Close();
    }

    return;
}
//////////////////////////////////////////////////////////////////////////////
