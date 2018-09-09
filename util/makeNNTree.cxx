//SusyNtuple
#include "SuperNtuple/NNMaker.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

//std/stl
#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

//ROOT
#include "TChain.h"

//////////////////////////////////////////////////////
//
// makeNNTree 
// Executable auto-generated with SusyNtuple/make_susy_skeleton on 2018-09-07 07:44
//
//
//////////////////////////////////////////////////////


void help()
{
    cout << "----------------------------------------------------------" << endl;
    cout << " makeNNTree" << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "   -n          number of events to process (default: all)" << endl;
    cout << "   -s|--suffix suffix to attach to end of output (default: NN)" << endl;
    cout << "   -d          debug level (integer) (default: 0)" << endl;
    cout << "   -i          input file (ROOT file, *.txt file, or directory)" << endl;
    cout << "   --nn        input NN json descriptor file for LWTNN" << endl;
    cout << "   --tree      name of TTree object in file (default: superNt)" << endl;
    cout << "   -h          print this help message" << endl;
    cout << endl;
    cout << "  Example Usage:" << endl;
    cout << "   makeNNTree -i susyNt.root -n 500 --nn <file>" << endl;
    cout << "----------------------------------------------------------" << endl;
}

int main(int argc, char** argv)
{

    /////////////////////////
    // cmd line options
    /////////////////////////

    int n_events = -1;
    int dbg = 0;
    string input = "";
    string nn_file = "";
    string suffix = "NN";

    for(int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-n") == 0) n_events = atoi(argv[++i]);
        else if (strcmp(argv[i], "-d") == 0) dbg = atoi(argv[++i]);
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--suffix") == 0) suffix = argv[++i];
        else if (strcmp(argv[i], "-i") == 0) input = argv[++i];
        else if (strcmp(argv[i], "--nn") == 0) nn_file = argv[++i];
        else if (strcmp(argv[i], "-h") == 0) { help(); return 0; }
        else {
            cout << "makeNNTree    Unknown command line argument '" << argv[i] << "', exiting" << endl;
            help();
            return 1;
        }
    } // i

    if(input.empty()) {
        cout << "makeNNTree    You must specify an input" << endl;
        return 1;
    }

    if(nn_file.empty()) {
        cout << "makeNNTree    You must specify an input NN JSON file" << endl;
        return 1;
    }


    /////////////////////////////////////////////////////////
    // Build the TChain object
    // For SusyNtuple analysis, the chain name is susyNt
    /////////////////////////////////////////////////////////
    TChain* chain = new TChain("superNt");

    // use ChainHelper to infer the input type (ROOT file, *.txt, or dir/)
    // and build the full chain of the input files
    // (c.f. SusyNtuple/ChainHelper.h)
    ChainHelper::addInput(chain, input, dbg>0);
    Long64_t n_entries_in_chain = chain->GetEntries();
    // let's see what it looks like
    chain->ls();

    /////////////////////////////////////////////////////////
    // Build the TSelector object
    // SusyNt analyses inheriting from SusyNtAna must
    // build their own TSelector looper
    /////////////////////////////////////////////////////////
    NNMaker* analysis = new NNMaker();

    analysis->set_debug(dbg);
    analysis->set_chain(chain); // propagate the TChain to the analysis
    analysis->set_suffix(suffix);
    analysis->load_nn(nn_file);
    analysis->set_input_name(input);
    if(n_events < 0) n_events = n_entries_in_chain;

    cout << "---------------------------------------------------------" << endl;
    cout << " Total entries in input chain          : " << n_entries_in_chain << endl;
    cout << " Total entries to process for analysis : " << n_events << endl;
    cout << "---------------------------------------------------------" << endl;
    
    // call TChain Process to star the TSelector looper over the input TChain
    //if(n_events > 0) chain->Process(analysis, input.c_str(), n_events);
    if(n_events > 0) analysis->process(n_events);

    cout << endl;
    cout << "makeNNTree    Analysis loop done" << endl;

    delete chain;
    return 0;
}
