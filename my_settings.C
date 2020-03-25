#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

using namespace std;

// For observable of interest
const Int_t filter_observable =8;// 8; // 0: TPC-only + refit, 1: itstpc2011
// filter 8: trackcuts2015 with open DCA
// Event selection
const Int_t trigger_sel = 1; //    kINT7=1, kCINT5=2
const Char_t *estimator = "CombTPCITS08"; // "CombTPCITS08", "V0M"
const Int_t multOjjects = 0; // 0 combinedTPCITStracklets, 1 TPC-only tracks 
const Bool_t isbinpercentile = kFALSE;// this is for multiplicity, it should be kTRUE when estimator = "V0M"
const Int_t  evShape = 0;//0 So, 1 St
const Bool_t isHybrid = kFALSE;  //Hybrid tracks in the So calculation
const Bool_t isGolden = kFALSE;  //Golden tracks in the So calculation, else TPConly+refit
const Bool_t IsMC = kTRUE;
const Bool_t IsESPerc = kTRUE;
const Bool_t IsSo4perc =kTRUE;

Bool_t removePileUp = kTRUE;
Bool_t removeBkg    = kTRUE;//cut on the correlation SDPclusters vs. tracklets?
Bool_t removePF     = kFALSE;//past and future protection
Float_t multiplicity = -10;
Float_t multiplicityMidEta = -10;

//For TPC-only
Bool_t fCutDCAToVertex2D = kTRUE;
Double_t fCutMaxDCAToVertexXY = 2.4;
Double_t fCutMaxDCAToVertexZ = 3.2;

// Spherocity settings

const Int_t nSobins = 10;
const Int_t nSobins2 = 8;

//Double_t Sobins[nSobins+1] = { 0.0, 0.3, 0.5, 0.8, 1.0 };
Double_t Sobins[nSobins+1] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
//Double_t Sobins[nSobins+1] = { 0.0, 0.2, 0.4, 0.8, 1.0 };
//Double_t Sobins[nSobins+1] = { 0.0, 0.1, 0.5, 0.9, 1.0 };

//Double_t Sobinst[12+1][nSobins2+1]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//Double_t Sobinsr[12+1][nSobins2+1]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

// Para Pythia Monash  como data



Int_t etacriteria=0;    // 0 for abs, 1 for positve eta, 2 for negative

Double_t size_step = 0.1;   //this is for the computation of spherocity
Double_t ptminSo = 0.15;
Double_t ptmaxSo = 1e8;
Int_t minMultSo = 3;
Int_t type_track = 0; // 0: no sel based on MC true, 1: only primary, 2: primary+weak decays


// Binning

const Int_t nMultbins2 = 13;
Double_t Multbins2[nMultbins2+1]={ 0.0, 1.0, 4.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 140.0 };


const Int_t nMultbins =100;
Double_t Multbins[nMultbins+1]={0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                               10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,
			       20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,
			       30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,
			       40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,
			       50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,
			       60.0,61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,
			       70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,
			       80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,
			       90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,
			      100.0 };

const Int_t nMultbinsTmp = 29;
Double_t MultbinsTmp[nMultbinsTmp+1]={ 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,35,40,45,50,55,60,65,70,75,80,100,120,140 };

const Int_t nV0Mbins = 11;
Double_t V0Mbins[nV0Mbins+1]={ 0.00, 0.10, 1.00, 5.00, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100 };

/*
const Int_t nPtBins      = 59;
Double_t xBins[nPtBins+1] = {

        0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
        1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
        4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20

};
*/

const Int_t nPtBins      = 42;
Double_t xBins[nPtBins+1] = {

        0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
        1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.5,4,5,10,20

};


const Char_t * inPercent = "./HistPercPP13000.root";
const Char_t * inPercentData = "./HistPercPP13000Data.root";
const Char_t * inSo4perc = "./SoIntervals_Per0asdat_100Nchbins.root";
TFile * finPercent;
TFile * finSo4perc;
TH1D  * hSOMPerc[nMultbins];
TH1D  * hSOTPerc[nMultbins];
TH1D  * hSTMPerc[nMultbins];
TH1D  * hSTTPerc[nMultbins];

TH1D * hHelperPerc;
TH1D * hVtxSpd;

// true mid-eta vs rec
TH2D * hRMmultCombAll0;
TH2D * hRMmultCombAllUn;
TH2D * hRMmultCombAll; // combined estimator
TH2D * hRMmultTPCoAll0; 
TH2D * hRMmultTPCoAll; // track counting TPC-only
TH2D * hTSoVsmultComb;
TH2D * hMSoVsmultComb;

TH2D * hSTSO[nMultbins];
TH2D * hSoRM[nMultbins];
TH2D * hSoRMNB[nMultbins];
TH2D * hSoRMNB2[nMultbins];
TH2D * hSompercVsSom[nMultbins];

TH2D * hRMmultComb[nSobins]; // combined estimator
TH2D * hRMmultTPCo[nSobins]; // track counting TPC-only
TH2D * hRMmultCombAllUnSom[nSobins]; 
TH2D * hRMmultCombAllUnSot[nSobins];


// Histos MakeSo
//
TH2D * hMidEtaMultVsMultInclusive;
TH1D * hcontadoreventos;

TH1D * hSoM[nMultbins];
TH1D * hSoT[nMultbins];

TH1D * hSOM[nMultbins];
TH1D * hSTM[nMultbins];
TH1D * hSOT[nMultbins];
TH1D * hSTT[nMultbins];


TH2D * hSotVsMult[nSobins];
TH2D * hMidEtaMultVsMultSot[nSobins];
TH2D * hSomVsMult[nSobins];
TH2D * hMidEtaMultVsMultSom[nSobins];
TH2D * hSotVsMultFake;
TH2D * hMidEtaMultVsMultSotFake;
TH2D * hSomVsMultFake;
TH2D * hMidEtaMultVsMultSomFake;

TH1D * hNtrk;
TH1D * hNch;
//------------changeH
TH1D * hMcOutAll;
TH1D * hMcOutAllnp; 
TH1D * hMcSecAll;

TH1D * hdcaXY;
TH1D * hdcaXYMC;
TH1D * hdcaXYwdec;
TH1D * hdcaXYsecmat;
TH1D * hdcaXYn[nMultbins2+1];
TH1D * hdcaXYMCn[nMultbins2+1];
TH1D * hdcaXYwdecn[nMultbins2+1];
TH1D * hdcaXYsecmatn[nMultbins2+1];
TH1D * hdcaXYns[nMultbins2+1][nSobins+1];
TH1D * hdcaXYMCns[nMultbins2+1][nSobins+1];
TH1D * hdcaXYwdecns[nMultbins2+1][nSobins+1];
TH1D * hdcaXYsecmatns[nMultbins2+1][nSobins+1];

TH1D * hMcOutMult[nMultbins];
TH1D * hMcOutMultnp[nMultbins];
TH2D * hPhiOutMult[nMultbins];
TH1D * hSoMInNtrk[nMultbins][nSobins+1];
TH1D * hSoTInNtrk[nMultbins][nSobins+1];
TH1D * hNtrkInNtrkSoM[nMultbins][nSobins+1]; 
TH1D * hMcOut[nMultbins][nSobins+1];
TH1D * hMcOutnp[nMultbins][nSobins+1];
TH2D * hPhiOut[nMultbins][nSobins+1];
TH1D * hMcSec[nMultbins][nSobins+1];

TH1D * hMcInAll;
TH1D * hMcInMult[nMultbins];
TH2D * hPhiInMult[nMultbins];
TH1D * hMcInMultTrue[nMultbins];
TH1D * hSoTInNch[nMultbins][nSobins+1];
TH1D * hNchInNchSoT[nMultbins][nSobins+1];
TH1D * hMcIn[nMultbins][nSobins+1];
TH2D * hPhiIn[nMultbins][nSobins+1];
TH1D * hMcTrue[nMultbins][nSobins+1];
