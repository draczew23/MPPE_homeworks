#include <stdio.h>

#include <TROOT.h>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>

#include <TCanvas.h>
#include <TStyle.h>

#include <TRandom3.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

const long N_GEN_EVENTS = 10000;

const double MEAN_p = 2.0; // [GeV]

float p1,p2;
float the1,the2;
float phi1,phi2;

long N_EVENTS = 0;

vector<Float_t> P1;
vector<Float_t> P2;

vector<Float_t> THE1;
vector<Float_t> THE2;

vector<Int_t> TRG1;
vector<Int_t> TRG2;

void read_MC(string file_name, long& N) {

  fstream in_file;

  float p1,p2;
  float the1,the2;

  int trg1,trg2;
  
  in_file.open(file_name, ios::in);

  N = 0;
  while(in_file >> p1 >> the1 >> trg1 >> p2 >> the2 >> trg2) {

    //cout << p1 << " " << the1 << " " << trg1 << endl;
    //cout << p2 << " " << the2 << " " << trg2 << endl;

    P1.push_back(p1);
    P2.push_back(p2);

    THE1.push_back(the1);
    THE2.push_back(the2);
    
    TRG1.push_back(trg1);
    TRG2.push_back(trg2);
    
    N++;

    if((N % 1000) == 0) {cout << "read_MC: N = " << N << endl;}
    
  }
  
  in_file.close();

}

void printx2(Float_t a)
{
    cout << a  << endl;
}


void first() {

  cout << endl;
  cout << "Hello world" << endl;
  cout << endl;  
}