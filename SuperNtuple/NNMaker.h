#ifndef SusyNtuple_NNMaker_h
#define SusyNtuple_NNMaker_h

//ROOT
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"

//SusyNtuple
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"

//std/stl
#include <fstream>
#include <map>

//lwtnn
#include "SuperNtuple/LightweightGraph.hh"
#include "SuperNtuple/parse_json.hh"

/////////////////////////////////////////////////////////////
//
// NNMaker
// Class auto-generated with SusyNtuple/make_susy_skeleton on 2018-09-07 07:44
//
/////////////////////////////////////////////////////////////

// for TSelector analysis loopers processing susyNt you MUST inherit from SusyNtAna
// in order to pick up the susyNt class objects
class NNMaker
{

    public :
        NNMaker();
        virtual ~NNMaker() {};

        void set_debug(int dbg) { m_dbg = dbg; }
        int dbg() { return m_dbg; }

        void set_chain(TChain* chain) { m_input_chain = chain; }
        TChain* chain() { return m_input_chain; }

        void set_suffix(std::string val) { m_suffix = val; }
        std::string suffix() { return m_suffix; }

        bool set_input_name(std::string input_filename);
        std::string input_filename() { return m_input_name; }
        std::string output_filename() { return m_output_name; }

        void reset();

        // load the provided neural network (for now just do one)
        bool load_nn(std::string filename);
        void setup_nn_branches();
        float feature(std::string name);
        std::map<std::string, double> lwt_map();
        lwt::LightweightGraph* graph() { return m_lwt_graph; }

        ////////////////////////////////////////////
        // analysis methods
        ////////////////////////////////////////////
        bool process(long int n_to_process = -1); // main process loop

        ////////////////////////////////////////////
        // TSelector methods override
        ////////////////////////////////////////////
        //virtual void Init(TTree* tree);
        //virtual void Begin(TTree* tree); // Begin is called before looping on entries
        //virtual Bool_t Process(Long64_t entry); // Main event loop function called on each event
        //virtual void Terminate(); // Terminate is called after looping has finished
        //virtual Int_t Version() const { return 2; }


    private :
        int m_dbg;
        std::string m_suffix;
        std::string m_input_name;
        TChain* m_input_chain; // the TChain object we are processing
        std::string m_lwt_nn_jsonfile;
        lwt::LightweightGraph* m_lwt_graph;

        std::string m_output_name;
        TFile* m_output_file;
        TTree* m_output_tree;

        TBranch* b_p_hh;
        TBranch* b_p_tt;
        TBranch* b_p_wt;
        TBranch* b_p_z;

        TBranch* b_d_hh;
        TBranch* b_d_tt;
        TBranch* b_d_wt;
        TBranch* b_d_z;

        TBranch* b_lr_hh_tt;
        TBranch* b_lr_hh_wt;
        TBranch* b_lr_hh_z;

        TBranch* b_pr_hh_tt;
        TBranch* b_pr_hh_wt;
        TBranch* b_pr_hh_z;
        TBranch* b_pr_hh_tt_wt;
        TBranch* b_pr_hh_tt_wt_z;

        double m_p_hh;
        double m_p_tt;
        double m_p_wt;
        double m_p_z;

        double m_d_hh;
        double m_d_tt;
        double m_d_wt;
        double m_d_z;

        double m_lr_hh_tt;
        double m_lr_hh_wt;
        double m_lr_hh_z;

        double m_pr_hh_tt;
        double m_pr_hh_wt;
        double m_pr_hh_z;
        double m_pr_hh_tt_wt;
        double m_pr_hh_tt_wt_z;

        std::map<std::string, float> feature_map_f;
        std::map<std::string, int> feature_map_i;
        float f_met;
        float f_metPhi;
        float f_mll;
        float f_dRll;
        float f_pTll;
        float f_dphi_ll;
        float f_dphi_bb;
        float f_dphi_met_ll;
        float f_met_pTll;
        float f_nJets;
        float f_nSJets;
        float f_nBJets;
        float f_isEE;
        float f_isMM;
        float f_isSF;
        float f_isDF;
        float f_HT2;
        float f_HT2Ratio;
        float f_l0_pt;
        float f_l1_pt;
        float f_l0_phi;
        float f_l1_phi;
        float f_j0_pt;
        float f_j1_pt;
        float f_j0_phi;
        float f_j1_phi;
        float f_bj0_pt;
        float f_bj1_pt;
        float f_bj0_phi;
        float f_bj1_phi;

        std::map<std::string, double> m_lwt_map;


}; //class


#endif
