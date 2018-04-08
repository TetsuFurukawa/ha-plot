#include "getExon.h"
#include <stdarg.h>
/*
FUSの解析を行うためのmainファイル
FUS_Anolis
FUS_Human
FUS_Ciona
FUS_Chimp
FUS_Danio
FUS_Gallus
FUS_Mus_Strain1
FUS_Mus_Strain2
FUS_Xenopus

*/

int main(int argc, char *argv[]){
  /*雛形その1(構造抽出)*/
  if(argc < 2){
    cout << "This program needs *.gp file." << endl;
    cout << "EXAMPLE:" << endl;
    cout << "   ./GE example.gp" << endl;
    cout << "This program is aborted." << endl;
    exit(1);
  }


  /*vector<string> spacies = { "Armadillo", "Chimp", "Ciona", "Danio",
                            "Gallus", "Gorilla", "Human", "Macaca", "Mus",
                            "Opossum", "Platypes", "Rattus", "Xenopus"};
  string proteinName = "TDP43";*/
  string filename = string(argv[1]);
  ReadFile r;
  r.filename = filename;
  vector<string> file = r.readFile();
  Extract e;
  e.content = file;
  e.outfile = filename + ".gnst";
  e.extractStructure();
  /*for(auto s : spacies){
    cout << s << endl;
    ReadFile r;
    r.filename = s+"/"+proteinName+"_"+s+".gb";
    vector<string> file = r.readFile();
    Extract e(proteinName);
    e.makeSequencePair(s);
    e.content = file;
    e.outfile = s+"/"+proteinName+"_"+s+"_gnst";
    e.extractStructure();
    e.addShell(s);
  }
*/
  //
  /*雛形その2(構造抽出+配列抽出)
  ReadFile r;
  r.filename = "FUS_Human.gb";
  vector<string> file = r.readFile();
  Extract e;
  e.content = file;
  e.outfile = "test";
  e.extractStructure();
  vector<int> area{11, 12, 13, 14, 15};
  e.which_area = area;
  e.isoform = 0;
  e.outfile_2 = "FUS_extract.fasta";
  e.extractSequence();
*/
  return 0;
}
