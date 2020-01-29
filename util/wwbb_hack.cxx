#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/parse_json.hh"

int main(int argc, char* argv[])
{
    if(argc<2)
    {
        cout << "error: must provide one argument for the LWTNN input" << endl;
        exit(1);
    }
    string nn_filename = argv[1];
    std::ifstream input_nn_file(nn_filename);
    string input_layer_name = "InputLayer";
    string output_layer_name = "OutputLayer";
    auto config = lwt::parse_json_graph(input_nn_file);
    lwt::LightweightGraph nn_graph(config, output_layer_name);
    std::map<std::string, double> nn_inputs;

    // imagine that this is a single event where you have filled these variables
    nn_inputs["met"] = 4.0;
    nn_inputs["metPhi"] = 4.0;
    nn_inputs["mll"] = 4.0;
    nn_inputs["dRll"] = 4.0;
    nn_inputs["pTll"] = 4.0;
    nn_inputs["dphi_ll"] = 4.0;
    nn_inputs["dphi_bb"] = 4.0;
    nn_inputs["dphi_met_ll"] = 4.0;
    nn_inputs["met_pTll"] = 4.0;
    nn_inputs["nJets"] = 4.0;
    nn_inputs["nSJets"] = 4.0;
    nn_inputs["nBJets"] = 4.0;
    nn_inputs["isEE"] = 4.0;
    nn_inputs["isMM"] = 4.0;
    nn_inputs["isSF"] = 4.0;
    nn_inputs["isDF"] = 4.0;
    nn_inputs["HT2"] = 4.0;
    nn_inputs["HT2Ratio"] = 4.0;
    nn_inputs["l0_pt"] = 4.0;
    nn_inputs["l1_pt"] = 4.0;
    nn_inputs["l0_phi"] = 4.0;
    nn_inputs["l1_phi"] = 4.0;
    nn_inputs["j0_pt"] = 4.0;
    nn_inputs["j1_pt"] = 4.0;
    nn_inputs["j0_phi"] = 4.0;
    nn_inputs["j1_phi"] = 4.0;
    nn_inputs["bj0_pt"] = 4.0;
    nn_inputs["bj1_pt"] = 4.0;
    nn_inputs["bj0_phi"] = 4.0;
    nn_inputs["bj1_phi"] = 4.0;
    nn_inputs["l0_eta"] = 4.0;
    nn_inputs["l1_eta"] = 4.0;
    nn_inputs["bj0_eta"] = 4.0;
    nn_inputs["bj1_eta"] = 4.0;
    nn_inputs["mt2_bb"] = 4.0;

    std::map< std::string, std::map<std::string, double> > layer_inputs;
    layer_inputs[input_layer_name] = nn_inputs;
    auto output_scores = nn_graph.compute(layer_inputs);
    double p_hh = output_scores.at("out_0_hh");
    double p_tt = output_scores.at("out_1_tt");
    double p_wt = output_scores.at("out_2_wt");
    double p_z  = output_scores.at("out_3_zll");

    cout << "=============================" << endl;
    cout << "NN scores for the event:" << endl;
    cout << "  p_hh = " << p_hh << endl;
    cout << "  p_tt = " << p_tt << endl;
    cout << "  p_wt = " << p_wt << endl;
    cout << "  p_z  = " << p_z << endl;
    cout << "=============================" << endl;

}
