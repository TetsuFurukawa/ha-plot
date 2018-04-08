#ifndef __getExon_H_INCLUDED__
#define __getExon_H_INCLUDED__
/*
GenBank形式(.gb)のファイルから、コーディング領域を抽出して、
HPLT6によるドットプロット解析に用いる-ginf1のデータを作成するプログラム

NC_000070.6 (148612382..148626996, complement)
ToDo

配列抽出のクラスを作成する
なるべく他のものの依存性が少なくなるようにしておく。
*/

#include <tuple>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <thread>
#include <list>
#include <iomanip>
#include <stdlib.h>

using namespace std;


class ReadFile{
public:
  string filename;
  vector<string> readFile();
};

vector<string> ReadFile::readFile(){
  ifstream ifst;
  ifst.open(filename);
  string buff;
  vector<string> lines;
  while(getline(ifst, buff)){
    lines.push_back(buff);
  }
  return lines;
}

class Extract{
public:
  string outfile;
  string outfile_2;
  vector<string> content;
  string proteinName;
  int structureCount = 0;
  Extract(){
  }
  Extract(string input){
    proteinName = input;
    structureCount = 0;
  }
  void extractStructure();
  void extractSequence();
  void addShell(string spacie);
  void makeSequencePair(string spacie);
  int isoform = 0;
  bool mRNAorNot = 1;//これがtrueのときはexonではなくmRNAの情報を利用する
  vector<string> annotation;
  vector<vector<int>> wholepins;//Intron, Exon, nc Exon全ての開始終了点を記述したベクトル
  vector<vector<vector<int>>> allPins;
  vector<vector<string>> allAnnotation;
  vector<int> which_area; //配列を抽出する場所
};

void Extract::extractStructure(){
  vector<int> exons;
  int gene_initial;
  int gene_end;
  vector<int> e_initials;
  vector<int> e_ends;
  vector<int> CDSs;
  vector<vector<int>> CDSS;
  vector<int> mRNAs;
  vector<vector<int>> mRNAS;
  int line = 0;
  bool gene_hit = 1;
  for(auto con:content){
    cout << con << endl;
    if(((int)con.size() > 10)&&(con.substr(5,4) == "gene")
    &&(con.substr(21, 1) != "c")){//遺伝子の領域を取り出す。
      string buff = con.substr(21, (int)con.size()-21);

      if(buff.find(">") != string::npos)
      buff.erase(buff.begin()+(int)buff.rfind(">"), buff.begin()+(int)buff.rfind(">")+1);
      if(buff.find("<") != string::npos)
      buff.erase(buff.begin()+(int)buff.rfind("<"), buff.begin()+(int)buff.rfind("<")+1);

      gene_initial = stoi(buff.substr(0, buff.find("..")));
      gene_end = stoi(buff.substr(buff.find("..") + 2,
       (int)buff.size() - (buff.find("..") + 2)));
      gene_hit = 0;
    }
    if(((int)con.size() > 10)&&(con.substr(5,3) == "CDS")
    &&(con.substr(21, 1) != "c")){//遺伝子の領域を取り出す。
      string CDS = "";
      CDSs.clear();
      //cout << con << endl;
      int line2 = line;
      string first_buff;
      if(content[line2].substr(21, 1) == "j"){
        first_buff = content[line2].substr(26, (int)con.size()-26);
      }else{
        first_buff = content[line2].substr(21, (int)con.size()-20);
      }
      //first_buff.pop_back();
      CDS += first_buff;
      line2++;
      while(content[line2].substr(21, 1) != "/"){
        string buff = content[line2].substr(21, (int)content[line2].size()-21);
        CDS += buff;
        line2++;
      }
      CDS.pop_back();
      if(CDS.find(">") != string::npos)
      CDS.erase(CDS.begin()+(int)CDS.rfind(">"), CDS.begin()+(int)CDS.rfind(">")+1);
      if(CDS.find("<") != string::npos)
      CDS.erase(CDS.begin()+(int)CDS.rfind("<"), CDS.begin()+(int)CDS.rfind("<")+1);
      int init_char = 0;
      while(CDS.find(',') < (int)CDS.size()){
        CDSs.push_back(stoi(CDS.substr(init_char, CDS.find(".."))));
        CDS.erase(CDS.begin() + init_char, CDS.begin() + CDS.find("..")+2);
        CDSs.push_back(stoi(CDS.substr(init_char, CDS.find(','))));
        CDS.erase(CDS.begin() + init_char, CDS.begin() + CDS.find(',') + 1);
      }
      CDSs.push_back(stoi(CDS.substr(init_char, CDS.find(".."))));
      CDS.erase(CDS.begin() + init_char, CDS.begin() + CDS.find("..")+2);
      CDSs.push_back(stoi(CDS));
      CDSS.push_back(CDSs);
      CDSs.clear();
    }
    if(((int)con.size() > 10)&&(con.substr(5,4) == "exon")){
      string buff = con.substr(21, (int)con.size()-21);
      if(buff.find(">") != string::npos)
      buff.erase(buff.begin()+(int)buff.rfind(">"), buff.begin()+(int)buff.rfind(">")+1);
      if(buff.find("<") != string::npos)
      buff.erase(buff.begin()+(int)buff.rfind("<"), buff.begin()+(int)buff.rfind("<")+1);
      int initial = stoi(buff.substr(0, buff.find("..")));
      int end = stoi(buff.substr(buff.find("..") + 2,
       (int)buff.size() - (buff.find("..") + 2)));
      exons.push_back(initial);
      exons.push_back(end);
    }
    if(((int)con.size() > 10)&&(con.substr(5,4) == "mRNA")
    &&(con.substr(21, 1) != "c")){//遺伝子の領域を取り出す。f
      string mRNA = "";
      mRNAs.clear();
      //cout << con << endl;
      int line2 = line;
      string first_buff;
      if(content[line2].substr(21, 1) == "j"){
        first_buff = content[line2].substr(26, (int)con.size()-26);
      }else{
        first_buff = content[line2].substr(21, (int)con.size()-20);
      }
      //first_buff.pop_back();
      mRNA += first_buff;
      line2++;
      while(content[line2].substr(21, 1) != "/"){
        string buff = content[line2].substr(21, (int)content[line2].size()-21);
        if(buff.find(">") != string::npos)
        buff.erase(buff.begin()+(int)buff.rfind(">"), buff.begin()+(int)buff.rfind(">")+1);
        if(buff.find("<") != string::npos)
        buff.erase(buff.begin()+(int)buff.rfind("<"), buff.begin()+(int)buff.rfind("<")+1);
        mRNA += buff;
        line2++;
      }
      mRNA.pop_back();
      int init_char = 0;
      while(mRNA.find(',') < (int)mRNA.size()){
        mRNAs.push_back(stoi(mRNA.substr(init_char, mRNA.find(".."))));
        mRNA.erase(mRNA.begin() + init_char, mRNA.begin() + mRNA.find("..")+2);
        mRNAs.push_back(stoi(mRNA.substr(init_char, mRNA.find(','))));
        mRNA.erase(mRNA.begin() + init_char, mRNA.begin() + mRNA.find(',') + 1);
      }
      mRNAs.push_back(stoi(mRNA.substr(init_char, mRNA.find(".."))));
      mRNA.erase(mRNA.begin() + init_char, mRNA.begin() + mRNA.find("..")+2);
      mRNAs.push_back(stoi(mRNA));
      mRNAS.push_back(mRNAs);
      mRNA.clear();
      mRNAs.clear();
    }
    line++;

  }

  cout << endl << "CDS:" << endl;
  for(auto c:CDSs){
    cout << c << endl;
  }
  cout << endl << "Exons:" << endl;
  for(auto e:exons){
    cout << e << endl;
  }
  cout << endl << "mRNAs:" << endl;
  for(auto m:mRNAs){
    cout << m << endl;
  }
  cout << endl << "CDSS: size = " << (int)mRNAS.size() << endl;

  int l = 0;
  for(auto C:CDSS){
    cout << l << ":" << endl;
    for(auto c:C){
      cout << c << endl;
    }
    cout << endl;
    l++;
    break;
  }
  l = 0;
  cout << endl << "mRNAS: size = " << (int)mRNAS.size() << endl;
  for(auto M:mRNAS){
    cout << l << ":" << endl;
    for(auto m:M){
      cout << m << endl;
    }
    cout << endl;
    l++;
    break;
  }

  
  //対応しない領域(gap)には0を代入してアラインメントする。
  // exon:CDSでアラインメントする
  if((int)mRNAS.size() > 1){
    mRNAs = mRNAS[isoform];
    CDSs = CDSS[isoform];
  }else{
    mRNAs = mRNAS[0]; CDSs = CDSS[0];
  }
  if(mRNAorNot){
    exons = mRNAs;
  }

  for(int p = 0; p < (int)mRNAS.size(); p++){
    //cout << "p = " << p << endl;
    exons = mRNAS[p];
    CDSs = CDSS[p];
  int index_of_CDSs = 0;
  int index_of_exons = 0;
  vector<vector<int>> nc_exon;

  int j = 0;
  while((index_of_CDSs < (int)CDSs.size())&&(index_of_exons < (int)exons.size())){
    /*cout << "CDS = "<< CDSs[index_of_CDSs] << endl;
    cout << "ioc = " << index_of_CDSs << endl;
    cout << "exon = " << exons[index_of_exons] << endl;
    cout << "ioe = " << index_of_exons << endl;*/
    vector<int> row;
    vector<int> row2;
    if(exons[index_of_exons] < CDSs[index_of_CDSs]){
      if(index_of_exons % 2 == 0){
        annotation.push_back("n");
        //cout << j << ":nce_1" << endl;
        row2.push_back(exons[index_of_exons]);
        if(exons[index_of_exons + 1] < CDSs[index_of_exons]){
          row2.push_back(exons[index_of_exons + 1]);

        }else{
          row2.push_back(CDSs[index_of_exons] - 1);
        }
      }else if(index_of_exons % 2 == 1){
        annotation.push_back("i");
        //cout << j << ":nc_1" << endl;
        row2.push_back(exons[index_of_exons]+1);
        row2.push_back(exons[index_of_exons+1]-1);
      }
      row.push_back(exons[index_of_exons]);
      row.push_back(0);
      index_of_exons++;
    }else if(exons[index_of_exons] == CDSs[index_of_CDSs]){
      if(index_of_CDSs % 2 == 0){
        annotation.push_back("c");
        //cout << j << ":c_2" << endl;
        row2.push_back(CDSs[index_of_CDSs]);
        row2.push_back(CDSs[index_of_CDSs+1]);
      }else{
        annotation.push_back("i");
        //cout << j << ":nc_2" << endl;
        row2.push_back(CDSs[index_of_CDSs]+1);
        row2.push_back(CDSs[index_of_CDSs+1]-1);
      }
      row.push_back(exons[index_of_exons]);
      row.push_back(CDSs[index_of_CDSs]);
      index_of_CDSs++;    index_of_exons++;
    }else if(exons[index_of_exons] > CDSs[index_of_CDSs]){
      if(index_of_CDSs % 2 == 0){
        annotation.push_back("c");
        //cout << j << ":c_3" << endl;
        row2.push_back(CDSs[index_of_CDSs]);
        row2.push_back(CDSs[index_of_CDSs+1]);
      }else{
        if(index_of_CDSs == (int)CDSs.size()-1){
          //cout << j << ":nce_4" << endl;
          annotation.push_back("n");
          row2.push_back(CDSs[index_of_CDSs]+1);
          row2.push_back(exons[index_of_exons]);
        }else{
          annotation.push_back("i");

          //cout << j << ":nc_3" << endl;
          row2.push_back(CDSs[index_of_CDSs]+1);
          row2.push_back(CDSs[index_of_CDSs+1]-1);
        }
      }
      row.push_back(0);
      row.push_back(CDSs[index_of_CDSs]);

      index_of_CDSs++;
    }
    //cout << row[0] << " : " << row[1] << endl;
    nc_exon.push_back(row);
    wholepins.push_back(row2);
    j++;
  }
  if((int)exons.size() - index_of_exons > 1 ){
    while(index_of_exons < (int)exons.size()){
      vector<int> row2;
      if(index_of_exons % 2 == 0){
        annotation.push_back("n");
        //cout << j << ":nce_5" << endl;
        row2.push_back(exons[index_of_exons]);
        row2.push_back(exons[index_of_exons + 1]);
      }else if(index_of_exons % 2 == 1){
        annotation.push_back("i");
        //cout << j << ":nc_5" << endl;
        row2.push_back(exons[index_of_exons]+1);
        row2.push_back(exons[index_of_exons+1]-1);
      }
      vector<int> row;
      row.push_back(exons[index_of_exons]);
      row.push_back(0);
      nc_exon.push_back(row);
      index_of_exons++;
      wholepins.push_back(row2);
    }
  }
  //分けた区画が、イントロンかエキソンかncエキソンかを区分する。
  int line_count=0;
  vector<int> exon_intron; //端っこのエキソンに挟まれているイントロンを挿入する
  for(int i = 0; i < (int)exons.size()-1; i++){
    exon_intron.push_back(exons[i]);
    if(i % 2 == 1){
      exon_intron.push_back(exons[i] + 1);
      exon_intron.push_back(exons[i+1] - 1);
    }
  }
  //cout << gene_initial << endl;
  //cout << gene_end << endl << endl;

  exon_intron.push_back(exons[(int)exons.size()-1]);
  vector<string> out;

  int k = 0;

  for(auto w:wholepins){
    //cout << annotation[k] << "\t" << w[0] << " : " << w[1] << endl;
    cout << setw(1) << left << annotation[k] << setw(6) << right << w[0]
    << setw(6) << right << w[1]  << endl;
    k++;
  }
  cout << endl;



  //ファイル出力
  ofstream ofst;
  //ofst.open(outfile + "_" + to_string(p) + ".txt", ios::trunc);
  ofst.open(outfile, ios::trunc);
  ofst << "\t" << (int)annotation.size() << "\t" << outfile << endl;
  k = 0;
  for(auto w:wholepins){
    //cout << annotation[k] << "\t" << w[0] << " : " << w[1] << endl;
    ofst << setw(1) << annotation[k] << setw(6) << right << to_string(w[0])
    << setw(6) << right << to_string(w[1]) << endl;
    k++;
  }
  allPins.push_back(wholepins);
  allAnnotation.push_back(annotation);
  wholepins.clear();
  annotation.clear();
  structureCount = p;
}
}

//読み込んだファイルから必要な配列情報のみを取得する。

void Extract::extractSequence(){
  cout << "::::extractSequence::::" << endl;
  string sequence;
  int line = 0;
  for(auto c:content){
    //cout << c << endl;
    if(c.substr(0,6) == "ORIGIN"){
      line++;
      while((content[line].substr(0,2) != "//")){
        sequence += content[line].substr(10, (int)content[line].size()-10);
        line++;
      }
      while(sequence.find(" ") != string::npos){
        sequence.erase(sequence.begin()+(int)sequence.rfind(" "),
        sequence.begin()+(int)sequence.rfind(" ")+1);
      }  break;
    }  line++;
  }
  string extract;
  string part;
  for(auto wh:which_area){
    cout << allPins[isoform][wh][0] << " : " << allPins[isoform][wh][1] << endl;
    extract += sequence.substr(allPins[isoform][wh][0] - 1,
       allPins[isoform][wh][1] - allPins[isoform][wh][0] + 1);
    part += to_string(wh) + "," ;
  }

  string header = ">" + content[0] + "\t" + ":extracted: isoform_index: "
  + to_string(isoform) + ", part:" + part;

  ofstream ofst;
  ofst.open(outfile_2, ios::trunc);
  for(int i = (int)extract.size()/70; i > 0; i--){
    extract.insert(i*70, "\n");
  }
  ofst << header << endl;
  ofst << extract << endl;

}

//HPLT6実行のためのシェルスクリプトファイルに書き込みを行うプログラム
//引数は生物種名
void Extract::addShell(string spacie){
  ofstream ofs;
  ofs.open("executeHPLT6_1.sh", ios::app);
  ofstream ofs1;
  ofs1.open("loadGNU.txt", ios::app);
  for(int i = 0; i < structureCount+1; i++){
    string inFileName = spacie + "/"+ proteinName + "_" + spacie + "_auto";
    string outFileName = "gp_files/" + proteinName + "_" + spacie + "_auto_"
                          + to_string(i);
    string script = "./HPLT6 -infile "+inFileName+".txt -outfile "+outFileName
                    +".gp -pltdata "+outFileName+".plt -ginf1 "+outfile+"_"
                    +to_string(i)+".txt -ginf2 "+outfile+"_"+to_string(i)
                    +".txt -shift 250 -window 10 -thr 0.9 -optaa 0";
    string script_with_human = "./HPLT6 -infile Human_"+spacie+"/"+proteinName
    +"_Human_"+spacie+".txt -outfile gp_files/"+proteinName+"_Human_"+spacie
    +"_0_"+to_string(i)+".gp -pltdata gp_files/"+proteinName+"_Human_"+spacie+"_0_"
    +to_string(i)+".plt -ginf1 Human/"+proteinName+"_Human_gnst_0.txt -ginf2 "
    +spacie+"/"+proteinName+"_"+spacie+"_gnst_"+to_string(i)
    +".txt -shift 250 -window 10 -thr 0.9 -optaa 0";

    ofs << script << endl;
    ofs << script_with_human << endl;

    string scriptGNU = "load \'" + outFileName + ".gp\'";
    ofs1 << scriptGNU << endl;
    string scriptGNU_with_human = "load \'gp_files/" + proteinName + "_human_" + spacie + "_0_"
                                  + to_string(i) + ".gp\'";
    ofs1 << scriptGNU_with_human << endl;

  }
}


class WriteFile{
public:
  vector<string> output;
  string filename;
  void writeFile();

};

void Extract::makeSequencePair(string spacie){
  ofstream ofst;
  ofst.open("catinateSequence.sh", ios::app);
  string inputFile = spacie + "/" + proteinName + "_" + spacie + ".fasta";
  string outputFile = spacie + "/" + proteinName + "_" + spacie + "_auto.txt";
  string script_1 = "cat " + inputFile + " " + inputFile + " > " + outputFile;
  ofst << script_1 << endl;
  //ここからはヒトとの保存度を比較を行うためのシェルスクリプトを記述するためのプログラム
  string dirName = "Human_" + spacie;
  string systemScript = "mkdir " + dirName;
  system(systemScript.c_str());
  string script_2 = "cat Human/" + proteinName + "_Human.fasta " + inputFile + " > "
                    + dirName + "/" + proteinName + "_Human_" + spacie + ".txt";
  ofst << script_2 << endl;
  ofst.close();
}

void WriteFile::writeFile(){
  ofstream ofst;
  ofst.open(filename, ios::trunc);
}

#endif
