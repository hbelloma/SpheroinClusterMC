#include <TFile.h>
#include <TRandom.h>
#include <TList.h>
#include <TTree.h>
#include <TProfile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include "AliLog.h"
#include <AliXRDPROOFtoolkit.h>
#include "DebugClassesMultESA2016.C"
#include "my_tools.C"
#include "my_settings.C"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;



Int_t GetMeasMultiplicityBin(Double_t multIn);
Int_t GetMeasMultiplicityBin2(Double_t multIn);
Int_t GetMultiplicityParticles(TClonesArray* TrackArray, Double_t etaCutlower, Double_t etaCutupper, Bool_t isMC);
Float_t GetMultiplicity(ESAEvent* event);
Float_t CutDCAToVertex2D(ESATrack* track,Double_t fCutMaxDCAToVertexXY,Double_t fCutMaxDCAToVertexZ, Bool_t fCutDCAToVertex2D);
Double_t GetESPerc( Double_t valES, Int_t multbin, Bool_t isSo, Bool_t isMC );
TString CompleteDirName( TString nameOut3 );
void DeclareHists();
void DeclareHistsData();
void DeclareHistsEqual();
Bool_t isGoodTrack( ESATrack* trackIn );
Double_t GetSphericity(TClonesArray* TrackArray, Double_t ptCutlower, Double_t ptCutupper, Double_t etaCutlower, Double_t etaCutupper, Int_t minMult, Bool_t isMC, Bool_t ispPb);
Double_t GetSpherocity(TClonesArray* TrackArray, Double_t ptCutlower, Double_t ptCutupper, Double_t etaCutlower, Double_t etaCutupper, Int_t minMult, Bool_t isMC, Bool_t ispPb);

//___________________________________________________
void AnalyzeData(const Char_t* dataFileName, const Char_t* outFileName, Int_t maxEvents,  Int_t nmaxfiles, Int_t startidx)
{

	if(isHybrid==kTRUE && isGolden==kTRUE){
		cout<<"                                          "<<endl;
		cout<<"!!!!!!!!!!!!!!!!   bye, bye, please select only one option; isHybrid and isGolden can not be true simoulateously"<<endl;
		return;
	}


	if(IsESPerc){
		finPercent = 0;
		if(IsMC)
			finPercent = TFile::Open(inPercent);
		else
			finPercent = TFile::Open(inPercentData);
		for( Int_t i_mult = 0; i_mult < nMultbins; ++i_mult ){
			hSOMPerc[i_mult] = 0;
			hSOMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOMPerc%d",i_mult));
			hSTMPerc[i_mult] = 0;
			hSTMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSTMPerc%d",i_mult));
		}
	}


	TString nameOut2 = "";
	if(IsESPerc)
		nameOut2 += "Perc";
	nameOut2 = CompleteDirName( nameOut2 );
	nameOut2 += "_maxJet";
	nameOut2 += Form("%1.1f",Sobins[1]);
	nameOut2 += "minIso";
	nameOut2 += Form("%1.1f",Sobins[3]);



	//create a directory
	CreateDir(nameOut2.Data());
	TFile* dataOutFile = new TFile(Form("%s/%s", nameOut2.Data(),outFileName), "RECREATE");


	//declare histograms
	DeclareHistsData();

	//Open files, create the tree
	TTree* Tree = 0;
	if(strstr(dataFileName, ".list"))
	{
		AliXRDPROOFtoolkit tool;
		TChain* chain = tool.MakeChain(dataFileName,"tree", "", nmaxfiles, startidx);
		//TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 8000);
		chain->Lookup();
		Tree = chain;
	} 
	else
	{
		TFile* dataFile = FindFileFresh(dataFileName);
		if(!dataFile) return;
		Tree = (TTree*)dataFile->Get("tree");
	}

	ESAEvent* event = 0;
	Tree->SetBranchAddress("event", &event);
	TClonesArray* trackArray = 0;
	if(isHybrid)
		Tree->SetBranchAddress("trackHybridPar"  , &trackArray);
	else
		Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);



	Int_t nEvents = Tree->GetEntries();
	cout << "Number of events: " << nEvents << endl;

	if(maxEvents>0 && maxEvents < nEvents)
	{
		nEvents = maxEvents;
		cout << "N events was reduced to: " << maxEvents << endl;
	}

	Int_t currentRun = 0;
	Int_t nBad = 0;


	//Start analysis

	for(Int_t n = 0; n < nEvents; n++)
	{
		Tree->GetEntry(n);

		if((n+1)%1000000==0)
			cout << "Event: " << n+1 << "/" << nEvents << endl;

		if(event->run == -1) { nBad++; continue; }
		if(event->run != currentRun) { cout << "New run: " << event->run << endl; currentRun = event->run; }
		hcontadoreventos->Fill(1);

		if(!(event->trig&trigger_sel)) 
			continue;
		hcontadoreventos->Fill(2);

		// Data Adquisition (DAQ)
		if( event->isincompletedaq ) //==0
			continue;
		hcontadoreventos->Fill(3);


		//if you want to remove pile-up using SPD
		if(removePileUp)
			if(event->pileup)
				continue;
		hcontadoreventos->Fill(4);


		//cut on SPDclusters vs. tracklets
		if(removeBkg)
			if( event->bgreject ) //==1
				continue;
		hcontadoreventos->Fill(5);


		//only good vertex
		if(event->vtxstatus2015pp<1) // only fill tracks for events with vtx inside cuts
			continue;		


		hcontadoreventos->Fill(6);

		if( TMath::Abs(event->zvtxSpd) > 10.0 )
			continue;

		hcontadoreventos->Fill(7);

		hVtxSpd->Fill(event->zvtxTrk);

		//mid-rapidity multiplicity
		multiplicityMidEta = -10;
		multiplicityMidEta = event->trackmult08;
		//multiplicty selection
		multiplicity = -10;
		multiplicity = GetMultiplicity(event);

		// rec multiplicity using TPC-only track cuts
		Int_t multR = 0;
		multR = GetMultiplicityParticles( trackArray, 0, 0.8, kFALSE );

		if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==1) )
			multiplicity = multR;



		// spherocity
		Double_t Som = -1; // measured
		Double_t SomPerc = -1; // measured

		Int_t binNtrk = GetMeasMultiplicityBin(multiplicity);

		if(evShape==0){ 
			Som = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );

			if(IsESPerc){
				SomPerc = GetESPerc(Som,binNtrk,kTRUE,kFALSE);
			}


		}
		if(evShape==1){
			Som = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
			if(IsESPerc){
				SomPerc = GetESPerc(Som,binNtrk,kFALSE,kFALSE);
			}

		}

		hNtrk->Fill(multiplicity);

		if(binNtrk>=0){

			Double_t Som2 = -1;

			if(evShape==0){
				Som2 = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
				hSTSO[binNtrk]->Fill(Som2,Som);


			}
			if(evShape==1){
				Som2 = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
				hSTSO[binNtrk]->Fill(Som,Som2);
			}


		}


		if( multiplicityMidEta > 0 ){
			hcontadoreventos->Fill(8);
			if(Som>=0)
				hMSoVsmultComb->Fill(multiplicityMidEta,Som);
		}


		Int_t BinSom = -1;

		for( Int_t i_so = 0; i_so < nSobins; ++i_so ){

			if(IsESPerc){
				if( SomPerc >= Sobins[i_so] && SomPerc < Sobins[i_so+1] )
					BinSom = i_so;
			}else{
				if( Som >= Sobins[i_so] && Som < Sobins[i_so+1] )
					BinSom = i_so;
			}
		}


		if(binNtrk>=0){
			hSoM[binNtrk]->Fill(Som);
			if( BinSom >= 0 && BinSom < nSobins ){
				hSoMInNtrk[binNtrk][BinSom]->Fill(Som);
				hNtrkInNtrkSoM[binNtrk][BinSom]->Fill(multiplicityMidEta);
			}

		}
		hMidEtaMultVsMultInclusive->Fill(multiplicity,multiplicityMidEta);

		//
		//                             LOOP OVER CHARGED PARTICLES
		//                                   

		const Int_t nTracks = trackArray->GetEntries();

		for(Int_t i = 0; i < nTracks; i++) {

			ESATrack* track = (ESATrack*)trackArray->At(i);

			if(filter_observable==0){
				if(!isGoodTrack(track))
					continue;
			}
			else{
				if(!(track->filter&filter_observable))
					continue;
			}


                        //changed by hector

                        Double_t eta_cms = track->eta;

                        if(etacriteria==0){
                        if(TMath::Abs(eta_cms)>0.8)
                                continue;

                        }

                        if(etacriteria==1){
                        if(eta_cms<0 || eta_cms>0.8)
                                continue;

                        }
                        if(etacriteria==2){
                        if(eta_cms>0 || eta_cms<-0.8)
                                continue;

                        }
                        //


			hMcOutAll->Fill(track->pt);

			if( binNtrk >= 0 ){
				hMcOutMult[binNtrk]->Fill(track->pt);
				if( BinSom >= 0 && BinSom < nSobins ){
					hMcOut[binNtrk][BinSom]->Fill(track->pt);
				}

			}

		}




	}//end loop over events

	dataOutFile->cd();
	dataOutFile->Write();
	dataOutFile->Close();
	delete dataOutFile;
	dataOutFile = 0;

	cout << "Nbad (runno == -1) : " << nBad << endl;

}

//_____________________________________________
void GetHistosForEqual(const Char_t* dataFileName, const Char_t* outFileName, Int_t maxEvents, Int_t nmaxfiles, Int_t startidx)
{

	if(isHybrid==kTRUE && isGolden==kTRUE){
		cout<<"                                          "<<endl;
		cout<<"!!!!!!!!!!!!!!!!   bye, bye, please select only one option; isHybrid and isGolden can not be true simoulateously"<<endl;
		return;
	}

	TString nameOut2 = "Equal_";
	nameOut2 = CompleteDirName( nameOut2 );

	//create a directory
	CreateDir(nameOut2.Data());
	TFile* dataOutFile = new TFile(Form("%s/%s", nameOut2.Data(),outFileName), "RECREATE");

	//declare histograms
	DeclareHistsEqual();

	//Open files, create the tree
	TTree* Tree = 0;
	if(strstr(dataFileName, ".list"))
	{
		AliXRDPROOFtoolkit tool;
		TChain* chain = tool.MakeChain(dataFileName,"tree", "", nmaxfiles, startidx);
		chain->Lookup();
		Tree = chain;
	} 
	else
	{
		TFile* dataFile = FindFileFresh(dataFileName);
		if(!dataFile) return;
		Tree = (TTree*)dataFile->Get("tree");
	}

	ESAEvent* event = 0;
	Tree->SetBranchAddress("event", &event);
	TClonesArray* trackArray = 0;
	if(isHybrid)
		Tree->SetBranchAddress("trackHybridPar"  , &trackArray);
	else
		Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);
	TClonesArray* trackArrayMC = 0;
	if(IsMC)
		Tree->SetBranchAddress("trackMC"  , &trackArrayMC);



	Int_t nEvents = Tree->GetEntries();
	cout << "Number of events: " << nEvents << endl;

	if(maxEvents>0 && maxEvents < nEvents)
	{
		nEvents = maxEvents;
		cout << "N events was reduced to: " << maxEvents << endl;
	}

	Int_t currentRun = 0;
	Int_t nBad = 0;


	//Start analysis

	for(Int_t n = 0; n < nEvents; n++)
	{
		Tree->GetEntry(n);

		if((n+1)%1000000==0)
			cout << "Event: " << n+1 << "/" << nEvents << endl;

		if(event->run == -1) { nBad++; continue; }
		if(event->run != currentRun) { cout << "New run: " << event->run << endl; currentRun = event->run; }


		hcontadoreventos->Fill(1);


		if(!(event->trig&trigger_sel)) 
			continue;
		hcontadoreventos->Fill(2);

		// Data Adquisition (DAQ)
		if( event->isincompletedaq ) //==0
			continue;
		hcontadoreventos->Fill(3);


		//if you want to remove pile-up using SPD
		if(removePileUp)
			if(event->pileup)
				continue;
		hcontadoreventos->Fill(4);


		//cut on SPDclusters vs. tracklets
		if(removeBkg)
			if( event->bgreject ) //==1
				continue;
		hcontadoreventos->Fill(5);


		//only good vertex
		if(event->vtxstatus2015pp<1) // only fill tracks for events with vtx inside cuts
			continue;		


		hcontadoreventos->Fill(6);

		if( TMath::Abs(event->zvtxSpd) > 10.0 )
			continue;

                if(IsMC)
                        if( TMath::Abs(event->zvtxMC) > 10.0 )
                                continue;



		hcontadoreventos->Fill(7);

		hVtxSpd->Fill(event->zvtxTrk);

		multiplicityMidEta = -10;
		multiplicityMidEta = event->trackmult08;
		multiplicity = -10;
		multiplicity = GetMultiplicity(event);

		// rec multiplicity using TPC-only track cuts
               Int_t multT = 0;
                if(IsMC)
                        multT = GetMultiplicityParticles( trackArrayMC, 0, 0.8, kTRUE );

		Int_t multR = 0;
		multR = GetMultiplicityParticles( trackArray, 0, 0.8, kFALSE );

		if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==1) )
			multiplicity = multR;


		// spherOcity
		Double_t SOm = -1.0;
		Double_t SOt = -1.0;
		// spherIcity
		Double_t STm = -1.0;
                Double_t STt = -1.0; 
		SOm = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
		STm = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
		if(IsMC){
                        SOt = GetSpherocity( trackArrayMC, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kTRUE, kFALSE );
                        STt = GetSphericity( trackArrayMC, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kTRUE, kFALSE );
                }

		Int_t binNtrk = -1; 
		Int_t binNch  = -1;
		binNtrk = GetMeasMultiplicityBin(multiplicity);
		hNtrk->Fill(multiplicity);
		if(binNtrk>=0){
			hSOM[binNtrk]->Fill(SOm);
			hSTM[binNtrk]->Fill(STm);
		}
                if(IsMC){

                        binNch = GetMeasMultiplicityBin(multT);
                        hNch->Fill(multT);
                        if( (multiplicityMidEta > 0) && (multT > 0) ){
                                hRMmultCombAll->Fill(multT,multiplicityMidEta);
                        }

                        if(binNch>=0){

                                hSOT[binNch]->Fill(SOt);
                                hSTT[binNch]->Fill(STt);

                        }

                }


	}//end loop over events

	dataOutFile->cd();
	dataOutFile->Write();
	dataOutFile->Close();
	delete dataOutFile;
	dataOutFile = 0;

	cout << "Nbad (runno == -1) : " << nBad << endl;

}
//------------------------
void GetSpherocitySpectrumData(const Char_t* dataFileName, const Char_t* outFileName, Int_t maxEvents,  Int_t nmaxfiles, Int_t startidx)
{
        if(isHybrid==kTRUE && isGolden==kTRUE){
                cout<<"                                          "<<endl;
                cout<<"!!!!!!!!!!!!!!!!   bye, bye, please select only one option; isHybrid and isGolden can not be true simoulateously"<<endl;
                return;
        }
        if(IsESPerc){
                finPercent = 0;
                finPercent = TFile::Open(inPercent);
                for( Int_t i_mult = 0; i_mult < nMultbins; ++i_mult ){
                        hSOMPerc[i_mult] = 0;
                        hSOMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOMPerc%d",i_mult));
                        hSOTPerc[i_mult] = 0;
                        hSOTPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOTPerc%d",i_mult));
                        hSTMPerc[i_mult] = 0;
                        hSTMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSTMPerc%d",i_mult));
                        hSTTPerc[i_mult] = 0;
                        hSTTPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSTTPerc%d",i_mult));
                }
        }

       if(IsSo4perc){
                cout<<"IsSo4perc"<<endl;
                finSo4perc = 0;
                finSo4perc = TFile::Open(inSo4perc);
                cout<<"Openfile"<<endl;
                for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
                   hSoRMNB[ibin] = 0;
                   hSoRMNB[ibin] = (TH2D *)finSo4perc->Get(Form("h2newbins%d",ibin));
		}
	}
	TString nameOut2 = "";
        if(IsESPerc)
                nameOut2 += "Perc";
        nameOut2 = CompleteDirName( nameOut2 );
        nameOut2 += "_maxJet";
        nameOut2 += Form("%1.1f",Sobins[1]);
        nameOut2 += "minIso";
        nameOut2 += Form("%1.1f",Sobins[3]);

        CreateDir(nameOut2.Data());
        TFile* dataOutFile = new TFile(Form("%s/%s", nameOut2.Data(),outFileName), "RECREATE");
        DeclareHists();
        TTree* Tree = 0;
        if(strstr(dataFileName, ".list"))
        {
                AliXRDPROOFtoolkit tool;
                TChain* chain = tool.MakeChain(dataFileName,"tree", "", nmaxfiles, startidx);
                chain->Lookup();
                Tree = chain;
        }
        else
        {
                TFile* dataFile = FindFileFresh(dataFileName);
                if(!dataFile) return;
                Tree = (TTree*)dataFile->Get("tree");
        }
       ESAEvent* event = 0;
        Tree->SetBranchAddress("event", &event);
        TClonesArray* trackArray = 0;
        if(isHybrid)
                Tree->SetBranchAddress("trackHybridPar"  , &trackArray);
        else
                Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);

        Int_t nEvents = Tree->GetEntries();
        cout << "Number of events: " << nEvents << endl;
        if(maxEvents>0 && maxEvents < nEvents)
        {
                nEvents = maxEvents;
                cout << "N events was reduced to: " << maxEvents << endl;
        }

        Int_t currentRun = 0;
        Int_t nBad = 0;

        for(Int_t n = 0; n < nEvents; n++)
        {
                Tree->GetEntry(n);

                if((n+1)%1000000==0)
                        cout << "Event: " << n+1 << "/" << nEvents << endl;

                if(event->run == -1) { nBad++; continue; }
                if(event->run != currentRun) { cout << "New run: " << event->run << endl; currentRun = event->run; }
                hcontadoreventos->Fill(1);
                if(!(event->trig&trigger_sel))
                        continue;
                hcontadoreventos->Fill(2);
               if( event->isincompletedaq ) //==0
                        continue;
                hcontadoreventos->Fill(3);
               if(removePileUp)
                        if(event->pileup)
                                continue;
                hcontadoreventos->Fill(4);
                if(removeBkg)
                        if( event->bgreject ) //==1
                                continue;
                hcontadoreventos->Fill(5);
                
                if(event->vtxstatus2015pp<1) // only fill tracks for events with vtx inside cuts
                        continue;
                hcontadoreventos->Fill(6);
                
                if( TMath::Abs(event->zvtxSpd) > 10.0 )
                        continue;

                hcontadoreventos->Fill(7);
                hVtxSpd->Fill(event->zvtxTrk);
               //  cout<<"calculating Nch"<<endl;
                multiplicityMidEta = -10;
                multiplicityMidEta = event->trackmult08;
                multiplicity = -10;
                multiplicity = GetMultiplicity(event);

                Int_t multR = 0;
                multR = GetMultiplicityParticles( trackArray, 0, 0.8, kFALSE );

                if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==1) )
                        multiplicity = multR;
                Double_t Som = -1; // measured
                Double_t SomPerc = -1; // measured

                Int_t binNtrk = GetMeasMultiplicityBin(multiplicity);
		Int_t binNtrk2 = GetMeasMultiplicityBin2(multiplicity);
               // cout<<"binNtrk="<<binNtrk<<endl;

                if(evShape==0){
                        Som = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
		        if(IsESPerc){
                                SomPerc = GetESPerc(Som,binNtrk,kTRUE,kFALSE);
			}
                }
                if(evShape==1){
                        Som = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
                        if(IsESPerc){
                                SomPerc = GetESPerc(Som,binNtrk,kFALSE,kFALSE);
                        }
                }
                hNtrk->Fill(multiplicity);
                if(binNtrk>=0 && binNtrk2>=0){
                        Double_t Som2 = -1;
                        if(evShape==0){
                                Som2 = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
                                hSTSO[binNtrk]->Fill(Som2,Som);
                        }
                        if(evShape==1){
                                Som2 = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
                                hSTSO[binNtrk]->Fill(Som,Som2);
                        }

                                    if(IsESPerc && Som>=0 && SomPerc>=0 ){
                                   hSompercVsSom[binNtrk]->Fill(Som,SomPerc);
                        }
                }

                Int_t BinSom = -1;
                Int_t BinSom1cut = -1;

                for( Int_t i_so = 0; i_so < nSobins; ++i_so ){
                        if(IsESPerc){
                               if( ( SomPerc >= Sobins[i_so] && SomPerc < Sobins[i_so+1] ) )  //changeH
                                 {
                                     BinSom = i_so;
                                     BinSom1cut = i_so;
                                 }
                        }else{
                             if( ( Som >= Sobins[i_so] && Som < Sobins[i_so+1] ) )  //changeH
                                 {       BinSom1cut = i_so;
                                        BinSom = i_so;
                                 }
                        }
                }


                if(binNtrk>=0 && binNtrk2>=0){
                        hSoM[binNtrk]->Fill(Som);
                        if( BinSom >= 0 && BinSom < nSobins ){ //nSobins changeH
                                hSoMInNtrk[binNtrk][BinSom1cut]->Fill(Som);
                        }
                }


                const Int_t nTracks = trackArray->GetEntries();
           //     cout<<"nTracks"<<nTracks<<endl;
                for(Int_t i = 0; i < nTracks; i++) {
                        ESATrack* track = (ESATrack*)trackArray->At(i);
                        if(filter_observable==0){
                                if(!isGoodTrack(track))
                                        continue;
                        }
                        else{
                                if(!(track->filter&filter_observable))
                                        continue;
				//	cout<<"no pass"<<endl;
                                
                        }
           //             cout<<"pass the filter"<<endl;
                        Double_t eta_cms = track->eta;
                        if(etacriteria==0){
                                if(TMath::Abs(eta_cms)>0.8)
                                        continue;
                        }
                        if(etacriteria==1){
                                if(eta_cms<0 || eta_cms>0.8)
                                        continue;
                        }
                        if(etacriteria==2){
                                if(eta_cms>0 || eta_cms<-0.8)
                                        continue;
                        }
                        hMcOutAllnp->Fill(track->pt);
                        hdcaXY->Fill(track->dcaxy);
                        if( multiplicityMidEta > 0  ){
                          hdcaXYn[binNtrk2]->Fill(track->dcaxy);
                          hMcOutMultnp[binNtrk]->Fill(track->pt);
                          if( BinSom >= 0 && BinSom < nSobins ){
                            hMcOutnp[binNtrk][BinSom]->Fill(track->pt);
                            hdcaXYns[binNtrk2][BinSom]->Fill(track->dcaxy);
                          }
                        }

                        //if(track->primary==1)
                                hMcOutAll->Fill(track->pt);


                        if( multiplicityMidEta > 0 ){
                                        hMcOutMult[binNtrk]->Fill(track->pt);
                                        hPhiOutMult[binNtrk]->Fill(track->pt,track->phi);
                                if( BinSom >= 0 && BinSom < nSobins ){ //nSobins changeH
                                                        hMcOut[binNtrk][BinSom1cut]->Fill(track->pt);
                                                        hPhiOut[binNtrk][BinSom]->Fill(track->pt,track->phi);

                                }

                        }

                }
                cout<<"end loop of tracks"<<endl; 

        }//end loop over events
        cout<<"end loop of events"<<endl;
        dataOutFile->cd();
        dataOutFile->Write();
        //for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
        //           hSoRMNB[ibin]->Write();
        //}
        dataOutFile->Close();
        delete dataOutFile;
        dataOutFile = 0;

        cout << "Nbad (runno == -1) : " << nBad << endl;

}
//_____________________________________________
void GetSpherocitySpectrum(const Char_t* dataFileName, const Char_t* outFileName, Int_t maxEvents,  Int_t nmaxfiles, Int_t startidx)
{

        if(isHybrid==kTRUE && isGolden==kTRUE){
		cout<<"                                          "<<endl;
		cout<<"!!!!!!!!!!!!!!!!   bye, bye, please select only one option; isHybrid and isGolden can not be true simoulateously"<<endl;
		return;
	}

	if(IsESPerc){
		finPercent = 0;
		finPercent = TFile::Open(inPercent);
		for( Int_t i_mult = 0; i_mult < nMultbins; ++i_mult ){
			hSOMPerc[i_mult] = 0;
			hSOMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOMPerc%d",i_mult));
			hSOTPerc[i_mult] = 0;
			hSOTPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOTPerc%d",i_mult));
			hSTMPerc[i_mult] = 0;
			hSTMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSTMPerc%d",i_mult));
			hSTTPerc[i_mult] = 0;
			hSTTPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSTTPerc%d",i_mult));
		}
	}
        
        if(IsSo4perc){
                cout<<"IsSo4perc"<<endl;
                finSo4perc = 0;
                finSo4perc = TFile::Open(inSo4perc);
                cout<<"Openfile"<<endl;
                for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
                   hSoRMNB[ibin] = 0;
                   hSoRMNB[ibin] = (TH2D *)finSo4perc->Get(Form("h2newbins%d",ibin));
                    //   cout<<"get ibin="  <<ibin<<endl;                
                   //hSoRMNB[ibin] = (TH2D *)hSoRMNB2->Clone(Form("hSoRMNB_%d",ibin));
                }
        }

	TString nameOut2 = "";
	if(IsESPerc)
		nameOut2 += "Perc";
	nameOut2 = CompleteDirName( nameOut2 );
	nameOut2 += "_maxJet";
	nameOut2 += Form("%1.1f",Sobins[1]);
	nameOut2 += "minIso";
	nameOut2 += Form("%1.1f",Sobins[3]);


	//create a directory
	CreateDir(nameOut2.Data());
	TFile* dataOutFile = new TFile(Form("%s/%s", nameOut2.Data(),outFileName), "RECREATE");

	//declare histograms
	DeclareHists();

	//Open files, create the tree
	TTree* Tree = 0;
	if(strstr(dataFileName, ".list"))
	{
		AliXRDPROOFtoolkit tool;
		TChain* chain = tool.MakeChain(dataFileName,"tree", "", nmaxfiles, startidx);
		//TChain* chain = tool.MakeChain(dataFileName,"tree", 0, 8000);
		chain->Lookup();
		Tree = chain;
	} 
	else
	{
		TFile* dataFile = FindFileFresh(dataFileName);
		if(!dataFile) return;
		Tree = (TTree*)dataFile->Get("tree");
	}

	ESAEvent* event = 0;
	Tree->SetBranchAddress("event", &event);
	TClonesArray* trackArray = 0;
	if(isHybrid)
		Tree->SetBranchAddress("trackHybridPar"  , &trackArray);
	else
		Tree->SetBranchAddress("trackGlobalPar"  , &trackArray);
	TClonesArray* trackArrayMC = 0;
	Tree->SetBranchAddress("trackMC"  , &trackArrayMC);

	Int_t nEvents = Tree->GetEntries();
	cout << "Number of events: " << nEvents << endl;

	if(maxEvents>0 && maxEvents < nEvents)
	{
		nEvents = maxEvents;
		cout << "N events was reduced to: " << maxEvents << endl;
	}

	Int_t currentRun = 0;
	Int_t nBad = 0;


	//Start analysis

	for(Int_t n = 0; n < nEvents; n++)
	{
		Tree->GetEntry(n);

		if((n+1)%1000000==0)
			cout << "Event: " << n+1 << "/" << nEvents << endl;

		if(event->run == -1) { nBad++; continue; }
		if(event->run != currentRun) { cout << "New run: " << event->run << endl; currentRun = event->run; }
		hcontadoreventos->Fill(1);
		if(!(event->trig&trigger_sel)) 
			continue;
		hcontadoreventos->Fill(2);

		// Data Adquisition (DAQ)
		if( event->isincompletedaq ) //==0
			continue;
		hcontadoreventos->Fill(3);


		//if you want to remove pile-up using SPD
		if(removePileUp)
			if(event->pileup)
				continue;
		hcontadoreventos->Fill(4);


		//cut on SPDclusters vs. tracklets
		if(removeBkg)
			if( event->bgreject ) //==1
				continue;
		hcontadoreventos->Fill(5);

		if(event->vtxstatus2015pp<1) // only fill tracks for events with vtx inside cuts
			continue;		
		hcontadoreventos->Fill(6);
		if( TMath::Abs(event->zvtxSpd) > 10.0 )
			continue;
		if( TMath::Abs(event->zvtxMC) > 10.0 )
			continue;

		hcontadoreventos->Fill(7);
		hVtxSpd->Fill(event->zvtxTrk);

		multiplicityMidEta = -10;
		multiplicityMidEta = event->trackmult08;
		multiplicity = -10;
		multiplicity = GetMultiplicity(event);

		Int_t multT = 0;
                multT = GetMultiplicityParticles( trackArrayMC, 0, 0.8, kTRUE );
		Int_t multR = 0;
		multR = GetMultiplicityParticles( trackArray, 0, 0.8, kFALSE );

		if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==1) )
			multiplicity = multR;

		// spherocity
		Double_t Som = -1; // measured
		Double_t Sot = -1;
		Double_t SomPerc = -1; // measured
		Double_t SotPerc = -1;

		Int_t binNtrk = GetMeasMultiplicityBin(multiplicity);
                Int_t binNtrk2 = GetMeasMultiplicityBin2(multiplicity);
		Int_t binNch  = GetMeasMultiplicityBin(multT);
	        Int_t binNch2  = GetMeasMultiplicityBin2(multT);
        //        cout<<"binNtrk"<<binNtrk<<endl;
        //        cout<<"binNch"<<binNch<<endl;
        //        cout<<"starting ESA sel"<<endl;
		if(evShape==0){ 
			Som = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
			Sot = GetSpherocity( trackArrayMC, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kTRUE, kFALSE );
          //              cout<<"Som="<<Som<<endl;
	//		cout<<"Sot="<<Sot<<endl;
			if(IsESPerc){
				SomPerc = GetESPerc(Som,binNtrk,kTRUE,kFALSE);
	                	SotPerc = GetESPerc(Sot,binNch,kTRUE,kTRUE);
			}

            //       cout<<"end So"<<endl;
		}
		if(evShape==1){
			Som = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
			Sot = GetSphericity( trackArrayMC, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kTRUE, kFALSE );
			if(IsESPerc){
				SomPerc = GetESPerc(Som,binNtrk,kFALSE,kFALSE);
				SotPerc = GetESPerc(Sot,binNch,kFALSE,kTRUE);
			}

		}

		hNtrk->Fill(multiplicity);
		hNch->Fill(multT);
              //   cout<<"pased eventshape"<<endl;
		if(binNtrk>=0 && binNtrk2>=0){
			Double_t Som2 = -1;
			if(evShape==0){
				Som2 = GetSphericity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
				hSTSO[binNtrk]->Fill(Som2,Som);
			}
			if(evShape==1){
				Som2 = GetSpherocity( trackArray, ptminSo, ptmaxSo, 0, 0.8, minMultSo, kFALSE, kFALSE );
				hSTSO[binNtrk]->Fill(Som,Som2);
			}

			hSoRM[binNtrk]->Fill(Sot,Som);
                        if(IsSo4perc){
		           hSoRMNB[binNtrk]->Fill(Sot,Som);
                        }
                                    if(IsESPerc && Som>=0 && SomPerc>=0 ){
                          	   hSompercVsSom[binNtrk]->Fill(Som,SomPerc);
			} 
		}

		   if( (multiplicityMidEta > 0) && (multT > 0) ){
			hcontadoreventos->Fill(8);
			hRMmultCombAll0->Fill(multT,multiplicityMidEta);
                        hRMmultCombAllUn->Fill(multT,multR);
			hRMmultCombAll->Fill(multT,multiplicityMidEta);
			hRMmultTPCoAll0->Fill(multT,multR);
			hRMmultTPCoAll->Fill(multT,multR);
			if(Sot>=0)
				hTSoVsmultComb->Fill(multiplicityMidEta,Sot);
			if(Som>=0)
				hMSoVsmultComb->Fill(multiplicityMidEta,Som);
                         
		   }
                
 
		Int_t BinSom = -1;
		Int_t BinSot = -1;
                Int_t BinSom1cut = -1;
                Int_t BinSot1cut = -1;

		for( Int_t i_so = 0; i_so < nSobins; ++i_so ){
			if(IsESPerc){
			//	if( ( SomPerc >= Sobins[i_so] && SomPerc < Sobins[i_so+1] ) && ( SotPerc >= Sobins[i_so] && SotPerc < Sobins[i_so+1] ) )  //changeH
			//		BinSom = i_so;
			//	if( ( SotPerc >= Sobins[i_so] && SotPerc < Sobins[i_so+1]) && ( SomPerc >= Sobins[i_so] && SomPerc < Sobins[i_so+1] ) ) //---- added m cut ------------changeH
			//		BinSot = i_so;
                                if( ( SomPerc >= Sobins[i_so] && SomPerc < Sobins[i_so+1] ) )  //changeH
                                 {   
				     BinSom = i_so;
                                     BinSom1cut = i_so;
                                 }
                                 if( ( SotPerc >= Sobins[i_so] && SotPerc < Sobins[i_so+1] ) )  //changeH
                                 {
				     BinSot1cut = i_so;
				     BinSot = i_so;
			         }
			}else{
			//	if( ( Som >= Sobins[i_so] && Som < Sobins[i_so+1] )&&( Sot >= Sobins[i_so] && Sot < Sobins[i_so+1]) ) //changeH
			//		BinSom = i_so;
			//	if( (Sot >= Sobins[i_so] && Sot < Sobins[i_so+1]) && ( Som >= Sobins[i_so] && Som < Sobins[i_so+1] ) ) //----added m cut ---------changeH
			//		BinSot = i_so;
                                if( ( Som >= Sobins[i_so] && Som < Sobins[i_so+1] ) )  //changeH
                                 {       BinSom1cut = i_so;
 					BinSom = i_so;
				 }
			
                                 if( ( Sot >= Sobins[i_so] && Sot < Sobins[i_so+1] ) )  //changeH
                                 {       BinSot1cut = i_so;
					BinSot = i_so;
			         }
			
			}
		}


		if(binNtrk>=0 && binNtrk2>=0){
			hSoM[binNtrk]->Fill(Som);
			if( BinSom >= 0 && BinSom < nSobins ){ //nSobins changeH
                                hSoMInNtrk[binNtrk][BinSom1cut]->Fill(Som);
                                hSoTInNtrk[binNtrk][BinSom1cut]->Fill(Sot);
				if( BinSot >= 0 && BinSot <= nSobins ){  //Added NEw changeH
					hNtrkInNtrkSoM[binNtrk][BinSom]->Fill(multiplicity);
				}
			}
		}
		 if(binNch>=0){
			hSoT[binNch]->Fill(Sot);
			if( BinSot >= 0 && BinSot < nSobins ){ //nSobins changeH
				hSoTInNch[binNch][BinSot1cut]->Fill(Sot);
				hNchInNchSoT[binNch][BinSot]->Fill(multT);
			}
		 } 
	

		if( BinSom >= 0 && BinSom < nSobins ){

			hSotVsMult[BinSom]->Fill(multiplicity,Sot); // Sot: -1 to 1
			hMidEtaMultVsMultSot[BinSom]->Fill(multiplicity,multiplicityMidEta);

			hRMmultComb[BinSom]->Fill(multT,multiplicityMidEta);
			hRMmultTPCo[BinSom]->Fill(multT,multR);
                        hRMmultCombAllUnSom[BinSom1cut]->Fill(multT,multR);


		}
		else{

			hSotVsMultFake->Fill(multiplicity,Sot); // Sot: -1 to 1
			hMidEtaMultVsMultSotFake->Fill(multiplicity,multiplicityMidEta);

		}
                
		hMidEtaMultVsMultInclusive->Fill(multiplicity,multiplicityMidEta);

		if( BinSot >= 0 && BinSot < nSobins ){

			hSomVsMult[BinSot]->Fill(multiplicity,Som); // Som: -1 to 1
			hMidEtaMultVsMultSom[BinSot]->Fill(multiplicity,multiplicityMidEta);
                        hRMmultCombAllUnSot[BinSot1cut]->Fill(multT,multR);
		}
		else{
			hSomVsMultFake->Fill(multiplicity,Som); // Som: -1 to 1
			hMidEtaMultVsMultSomFake->Fill(multiplicity,multiplicityMidEta);
		}
		//
		//                             LOOP OVER CHARGED PARTICLES
		//
		//
//		cout<< "statring loop of MC tracks" <<endl
		const Int_t nMcTracks = trackArrayMC->GetEntries();
	 //       Printf("Enter loop of tracks");
         	for(Int_t i = 0; i < nMcTracks; i++) {

			ESATrackMC* trackMC = (ESATrackMC*)trackArrayMC->At(i);

			if(trackMC->qMC == 0)
				continue;

			// Lines added by Hector
			Double_t eta_cms = trackMC->etaMC;

			if(etacriteria==0){
				if(TMath::Abs(eta_cms)>0.8)
					continue;

			}

			if(etacriteria==1){
				if(eta_cms<0 || eta_cms>0.8)
					continue; 

			}
			if(etacriteria==2){
				if(eta_cms>0 || eta_cms<-0.8)
					continue;

			}

			
			//
                        // ?if(trackMC->primary==1)
			hMcInAll->Fill(trackMC->ptMC);
			// all charged particles, true info
			if(binNch>=0){
				hMcInMultTrue[binNch]->Fill(trackMC->ptMC);
				if( BinSot >= 0 && BinSot < nSobins ){  //nSobins changeH
					hMcTrue[binNch][BinSot]->Fill(trackMC->ptMC);
				}
			}
			if( (multiplicityMidEta > 0) && (multT > 0) ){
				hMcInMult[binNtrk]->Fill(trackMC->ptMC);
				hPhiInMult[binNtrk]->Fill(trackMC->ptMC,trackMC->phiMC);
				if( BinSom >= 0 && BinSom < nSobins ){ //nSobins changeH
					if( BinSot >= 0 && BinSot < nSobins ){
                                                //if(trackMC->primary==1){
						hMcIn[binNtrk][BinSom1cut]->Fill(trackMC->ptMC);
						hPhiIn[binNtrk][BinSom]->Fill(trackMC->ptMC,trackMC->phiMC);
					}
				}

			}

		}

		const Int_t nTracks = trackArray->GetEntries();
		for(Int_t i = 0; i < nTracks; i++) {
			ESATrack* track = (ESATrack*)trackArray->At(i);
      			if(filter_observable==0){
				if(!isGoodTrack(track))
					continue;
			}
			else{
				if(!(track->filter&filter_observable))
					continue;
			} 
			Double_t eta_cms = track->eta;
			if(etacriteria==0){
				if(TMath::Abs(eta_cms)>0.8)
					continue;
			}
			if(etacriteria==1){
				if(eta_cms<0 || eta_cms>0.8)
					continue;
			}
			if(etacriteria==2){
				if(eta_cms>0 || eta_cms<-0.8)
					continue;
			}
 			hMcOutAllnp->Fill(track->pt);
                        if( (multiplicityMidEta > 0) && (multT > 0) ){
  			  hMcOutMultnp[binNtrk]->Fill(track->pt);
 			  if( BinSom >= 0 && BinSom < nSobins ){
			    if( BinSot >= 0 && BinSot < nSobins ){
			    hMcOutnp[binNtrk][BinSom]->Fill(track->pt);
			    }
			  }
                        }

			if(track->primary==1){
				hMcOutAll->Fill(track->pt);
                                hdcaXYMC->Fill(track->dcaxy);
                        }
			else{
			    hMcSecAll->Fill(track->pt);
                            if(track->primary==2)
			      hdcaXYwdec->Fill(track->dcaxy);
			    if(track->primary==3)
                              hdcaXYsecmat->Fill(track->dcaxy);
                        }
                       
			if( (multiplicityMidEta > 0) && (multT > 0) ){
				if(track->primary==1){           
					hMcOutMult[binNtrk]->Fill(track->pt);
					hPhiOutMult[binNtrk]->Fill(track->pt,track->phi);
                                        hdcaXYMCn[binNtrk2]->Fill(track->dcaxy);
				}
                                else{
					if(track->primary==2)
                            		   hdcaXYwdecn[binNtrk2]->Fill(track->dcaxy);
                            		if(track->primary==3)
                            		   hdcaXYsecmatn[binNtrk2]->Fill(track->dcaxy);
				}
				if( BinSom >= 0 && BinSom < nSobins ){ //nSobins changeH
					if( BinSot >= 0 && BinSot < nSobins ){
						if(track->primary==1){
							hMcOut[binNtrk][BinSom]->Fill(track->pt);
							hPhiOut[binNtrk][BinSom]->Fill(track->pt,track->phi);
							hdcaXYMCns[binNtrk2][BinSom]->Fill(track->dcaxy);
						}
						else{
							hMcSec[binNtrk][BinSom]->Fill(track->pt);
						    if(track->primary==2)
                                                        hdcaXYwdecns[binNtrk2][BinSom]->Fill(track->dcaxy);
                                                    if(track->primary==3)
                                                         hdcaXYsecmatns[binNtrk2][BinSom]->Fill(track->dcaxy);						

						}

					}
				}

			}

		}// end loop overp
          //Printf("end loop of part");
	}//end loop over events
        //Printf("end loop of events");
	dataOutFile->cd();
	dataOutFile->Write();
        for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
                   hSoRMNB[ibin]->Write();;
        }
	dataOutFile->Close();
	delete dataOutFile;
	dataOutFile = 0;

	cout << "Nbad (runno == -1) : " << nBad << endl;
}
//______________________________________
Double_t GetESPerc( Double_t valES, Int_t multbin, Bool_t isSo, Bool_t isMC ){

	if(multbin<0)
		return -1.0;

	Double_t porcentaje = -1.0;
	if(isSo){
		if(isMC){
			for(Int_t bin = 1; bin <= hSOTPerc[multbin]->GetNbinsX(); ++bin){

				Double_t minEs = 0.0;
				Double_t maxEs = 0.0;
				minEs = hSOTPerc[multbin]->GetBinLowEdge(bin);
				maxEs = minEs + hSOTPerc[multbin]->GetBinWidth(bin);
				if( (valES >= minEs) && (valES < maxEs) )
					porcentaje = hSOTPerc[multbin]->GetBinContent(bin);

			}
		}else{
			for(Int_t bin = 1; bin <= hSOMPerc[multbin]->GetNbinsX(); ++bin){

				Double_t minEs = 0.0;
				Double_t maxEs = 0.0;
				minEs = hSOMPerc[multbin]->GetBinLowEdge(bin);
				maxEs = minEs + hSOMPerc[multbin]->GetBinWidth(bin);
				if( (valES >= minEs) && (valES < maxEs) )
					porcentaje = hSOMPerc[multbin]->GetBinContent(bin);

			}
		}
	}
	else{
		if(isMC){
			for(Int_t bin = 1; bin <= hSTTPerc[multbin]->GetNbinsX(); ++bin){

				Double_t minEs = 0.0;
				Double_t maxEs = 0.0;
				minEs = hSTTPerc[multbin]->GetBinLowEdge(bin);
				maxEs = minEs + hSTTPerc[multbin]->GetBinWidth(bin);
				if( (valES >= minEs) && (valES < maxEs) )
					porcentaje = hSTTPerc[multbin]->GetBinContent(bin);

			}
		}else{
			for(Int_t bin = 1; bin <= hSTMPerc[multbin]->GetNbinsX(); ++bin){

				Double_t minEs = 0.0;
				Double_t maxEs = 0.0;
				minEs = hSTMPerc[multbin]->GetBinLowEdge(bin);
				maxEs = minEs + hSTMPerc[multbin]->GetBinWidth(bin);
				if( (valES >= minEs) && (valES < maxEs) )
					porcentaje = hSTMPerc[multbin]->GetBinContent(bin);

			}
		}


	}


	return porcentaje;

}

//____________________________________________
void DeclareHistsEqual(){

	hcontadoreventos = 0;
	hcontadoreventos = new TH1D("hcontadoreventos","",10,0,10);

	hVtxSpd = 0;
	hVtxSpd = new TH1D("hVtxSpd","",4000,-20,20);

	hNtrk = 0;
	hNtrk = new TH1D("hNtrk","",300,0,300);

	for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
		hSOM[ibin] = 0;
		hSOM[ibin] = new TH1D(Form("hSOM%d",ibin),"measured spherocity", 1000,0,1.0);
		hSTM[ibin] = 0;
		hSTM[ibin] = new TH1D(Form("hSTM%d",ibin),"measured sphericity", 1000,0,1.0);
	}

	if(IsMC){

		hNch = 0;
		hNch = new TH1D("hNch","",300,0,300);

		hRMmultCombAll = 0;
		hRMmultCombAll = new TH2D( "hRMmultCombAll", "Combined estimator;d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
				300,0,300,300,0,300 );

		for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
			hSOT[ibin] = 0;
			hSOT[ibin] = new TH1D(Form("hSOT%d",ibin),"true spherocity", 1000,0,1.0);
			hSTT[ibin] = 0;
			hSTT[ibin] = new TH1D(Form("hSTT%d",ibin),"true sphericity", 1000,0,1.0);
		}

	}


}
//______________________________________
void DeclareHistsData(){

	hHelperPerc = 0;
	hHelperPerc = new TH1D("hHelperPerc","",nSobins,Sobins);

	hVtxSpd = 0;
	hVtxSpd = new TH1D("hVtxSpd","",4000,-20,20);

	hNtrk = 0;
	hNtrk = new TH1D("hNtrk","",nMultbins, Multbins);

	hMSoVsmultComb = 0;
	hMSoVsmultComb = new TH2D( "hMSoVsmultComb", "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbins, Multbins, nSobins, Sobins );

	hMidEtaMultVsMultInclusive = 0;
	hMidEtaMultVsMultInclusive = new TH2D( "hMidEtaMultVsMultInclusive", "",
			nV0Mbins, V0Mbins, nMultbins, Multbins );

	hcontadoreventos = 0;
	hcontadoreventos = new TH1D("hcontadoreventos","",10,0,10);

	hMcOutAll = 0;
	hMcOutAll = new TH1D("hMcOutAll","",nPtBins,xBins);

        hMcOutAllnp = 0;
        hMcOutAllnp = new TH1D("hMcOutAllnp","",nPtBins,xBins);
	// arrays vs multiplicity
	for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){

		hSoM[ibin] = 0;
		hSoM[ibin] = new TH1D(Form("hSoM%d",ibin),"", 1000, 0.0, 1.0);


		hSTSO[ibin] = 0;
		hSTSO[ibin] = new TH2D( Form("hSTvsSOMult%d",ibin), Form("Mult bin %d;#it{S}_{T}; #it{S}_{O}", ibin ),
				40, 0.0, 1.0, 40, 0.0, 1.0 );


		hMcOutMult[ibin] = 0;
		hMcOutMult[ibin] = new TH1D(Form("hMcOutMult%d",ibin),"",nPtBins,xBins);

                hMcOutMultnp[ibin] = 0;
                hMcOutMultnp[ibin] = new TH1D(Form("hMcOutMultnp%d",ibin),"",nPtBins,xBins);

		for(Int_t i_So = 0; i_So < nSobins; ++i_So ){ //changeH


			hSoMInNtrk[ibin][i_So] = 0;
			hSoMInNtrk[ibin][i_So]  = new TH1D(Form("hESMInNtrk%dES%d",ibin,i_So),"",1000,0.0,1.0);

			hNtrkInNtrkSoM[ibin][i_So] = 0;
			hNtrkInNtrkSoM[ibin][i_So] = new TH1D(Form("hNtrkInNtrkESM%dES%d",ibin,i_So),"",300,0,300);


			hMcOut[ibin][i_So] = 0;
			hMcOut[ibin][i_So] = new TH1D(Form("hMcOutMult%dES%d",ibin,i_So),"",nPtBins,xBins);
 
       			hMcOutnp[ibin][i_So] = 0;
                        hMcOutnp[ibin][i_So] = new TH1D(Form("hMcOutMultnp%dES%d",ibin,i_So),"",nPtBins,xBins);

		}


	}

}
//________________________________________
void DeclareHists(){


	hHelperPerc = 0;
	hHelperPerc = new TH1D("hHelperPerc","",nSobins,Sobins);

	hVtxSpd = 0;
	hVtxSpd = new TH1D("hVtxSpd","",4000,-20,20);

	hNtrk = 0;
	hNtrk = new TH1D("hNtrk","",nMultbins, Multbins);
	hNch = 0;
	hNch = new TH1D("hNch","",nMultbins, Multbins);

	hRMmultCombAll = 0;
	hRMmultCombAll = new TH2D( "hRMmultCombAll", "Combined estimator;d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbins, Multbins, nMultbins, Multbins );

        hRMmultCombAllUn = 0;
        hRMmultCombAllUn = new TH2D( "hRMmultCombAllUn", "Combined estimator;d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
                        200, 0, 200, 200, 0, 200 );

	hRMmultCombAll0 = 0;
	hRMmultCombAll0 = new TH2D( "hRMmultCombAll0", "Combined estimator;d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbinsTmp, MultbinsTmp, 200, 0, 200 );


	hRMmultTPCoAll = 0;
	hRMmultTPCoAll = new TH2D( "hRMmultTPCoAll", "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbins, Multbins, nMultbins, Multbins );


	hRMmultTPCoAll0 = 0;
	hRMmultTPCoAll0 = new TH2D( "hRMmultTPCoAll0", "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbinsTmp, MultbinsTmp, 200, 0, 200 );


	hTSoVsmultComb = 0;
	hTSoVsmultComb = new TH2D( "hTSoVsmultComb", "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbins, Multbins, nSobins, Sobins );

	hMSoVsmultComb = 0;
	hMSoVsmultComb = new TH2D( "hMSoVsmultComb", "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
			nMultbins, Multbins, nSobins, Sobins );


	hMidEtaMultVsMultInclusive = 0;
	hMidEtaMultVsMultInclusive = new TH2D( "hMidEtaMultVsMultInclusive", "",
			nV0Mbins, V0Mbins, nMultbins, Multbins );

	hcontadoreventos = 0;
	hcontadoreventos = new TH1D("hcontadoreventos","",10,0,10);

	//hPhi = 0;
	//hPhi = new TH1D("hPhi","",64,0,2*TMath::Pi());

	//hEta = 0;
	//hEta = new TH1D("hEta","",20,-1,1);

	hMcOutAll = 0;
	hMcOutAll = new TH1D("hMcOutAll","",nPtBins,xBins);
     
        hMcOutAllnp = 0;
        hMcOutAllnp = new TH1D("hMcOutAllnp","",nPtBins,xBins);

	hMcSecAll = 0;
	hMcSecAll = new TH1D("hMcSecAll","",nPtBins,xBins);

	hMcInAll = 0;
	hMcInAll = new TH1D("hInOutAll","",nPtBins,xBins);

        hdcaXY= new TH1D("hdcaXY", "dcaxy distribution ; dca xy  [cm]; Counts", 800, -4, 4);
        hdcaXYMC= new TH1D("hdcaXYMC", "dcaxyMC distribution ; dca xy  [cm]; Counts", 800, -4, 4);
        hdcaXYsecmat= new TH1D("hdcaXYsecmat", "dcaxysec distribution ; dca xy  [cm]; Counts", 800, -4, 4);
        hdcaXYwdec= new TH1D("hdcaXYwdec", "dcaxywdec distribution ; dca xy  [cm]; Counts", 800, -4, 4);
        for( Int_t ibin = 0; ibin < nMultbins2+1; ++ibin ){
          hdcaXYn[ibin]= new TH1D(Form("hdcaXYn%d",ibin), "dcaxy distribution ; dca xy  [cm]; Counts", 800, -4, 4);
          hdcaXYMCn[ibin]= new TH1D(Form("hdcaXYMCn%d",ibin), "dcaxyMC distribution ; dca xy  [cm]; Counts", 800, -4, 4);
          hdcaXYsecmatn[ibin]= new TH1D(Form("hdcaXYsecmatn%d",ibin), "dcaxysec distribution ; dca xy  [cm]; Counts", 800, -4, 4);
          hdcaXYwdecn[ibin]= new TH1D(Form("hdcaXYwdecn%d",ibin), "dcaxywdec distribution ; dca xy  [cm]; Counts", 800, -4, 4);
	  for(Int_t i_So = 0; i_So < nSobins; ++i_So ){
            hdcaXYns[ibin][i_So]= new TH1D(Form("hdcaXYns%d_%d",ibin,i_So), "dcaxy distribution ; dca xy  [cm]; Counts", 800, -4, 4);
            hdcaXYMCns[ibin][i_So]= new TH1D(Form("hdcaXYMCns%d_%d",ibin,i_So), "dcaxyMC distribution ; dca xy  [cm]; Counts", 800, -4, 4);
	    hdcaXYsecmatns[ibin][i_So]= new TH1D(Form("hdcaXYsecmatns%d_%d",ibin,i_So), "dcaxysec distribution ; dcaMC xy  [cm]; Counts", 800, -4, 4);
            hdcaXYwdecns[ibin][i_So]= new TH1D(Form("hdcaXYwdecns%d_%d",ibin,i_So), "dcaxywdec distribution ; dcaMC xy  [cm]; Counts", 800, -4, 4);
 	  }
	}       
        
	for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){

		hSoM[ibin] = 0;
		hSoM[ibin] = new TH1D(Form("hSoM%d",ibin),"", 1000, 0.0, 1.0);

		hSoT[ibin] = 0;
		hSoT[ibin] = new TH1D(Form("hSoT%d",ibin),"", 1000, 0.0, 1.0);

		hSTSO[ibin] = 0;
		hSTSO[ibin] = new TH2D( Form("hSTvsSOMult%d",ibin), Form("Mult bin %d;#it{S}_{T}; #it{S}_{O}", ibin ),
				40, 0.0, 1.0, 40, 0.0, 1.0 );

		hSoRM[ibin] = 0;
		hSoRM[ibin] = new TH2D( Form("hesRMMult%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),                      
				40, 0.0, 1.0, 40, 0.0, 1.0 );
               //new 
/*    
		hSoRMNB[ibin] = 0;
                if(ibin==0)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc0, nSobins, SobinsMRsopc0 );
                if(ibin==1)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc1, nSobins, SobinsMRsopc1 );                
                if(ibin==2)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc2, nSobins, SobinsMRsopc2 );
                if(ibin==3)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc3, nSobins, SobinsMRsopc3 );
                if(ibin==4)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc4, nSobins, SobinsMRsopc4 );
                if(ibin==5)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc5, nSobins, SobinsMRsopc5 );
                if(ibin==6)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc6, nSobins, SobinsMRsopc6 );
                if(ibin==7)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc7, nSobins, SobinsMRsopc7 );
                if(ibin==8)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc8, nSobins, SobinsMRsopc8 );
                if(ibin==9)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc9, nSobins, SobinsMRsopc9 );
                if(ibin==10)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc10, nSobins, SobinsMRsopc10 );
                if(ibin==11)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc11, nSobins, SobinsMRsopc11 );
                if(ibin==12)
                hSoRMNB[ibin] = new TH2D( Form("hesRMMultNB%d",ibin), Form("Mult bin %d; true event shape; measured event shape", ibin ),
                                nSobins, SobinsMRsopc12, nSobins, SobinsMRsopc12 );
              
  */   
                hSompercVsSom[ibin] = new TH2D( Form("hSompercVsSom%d",ibin), "", 
                                 1000,0,1.0, nSobins, Sobins); 


		hMcOutMult[ibin] = 0;
		hMcOutMult[ibin] = new TH1D(Form("hMcOutMult%d",ibin),"",nPtBins,xBins);

   		hMcOutMultnp[ibin] = 0;
                hMcOutMultnp[ibin] = new TH1D(Form("hMcOutMultnp%d",ibin),"",nPtBins,xBins);

		hPhiOutMult[ibin] = 0;
		hPhiOutMult[ibin] = new TH2D(Form("hMcPhiOutMult%d",ibin),"",nPtBins,xBins,64,0,2.0*TMath::Pi());

		hMcInMult[ibin]  = 0;
		hMcInMult[ibin]  = new TH1D(Form("hMcInMult%d",ibin),"",nPtBins,xBins);

		hPhiInMult[ibin]  = 0;
		hPhiInMult[ibin]  = new TH2D(Form("hMcPhiInMult%d",ibin),"",nPtBins,xBins,64,0,2.0*TMath::Pi());

		hMcInMultTrue[ibin]  = 0;
		hMcInMultTrue[ibin]  = new TH1D(Form("hMcInMultTrue%d",ibin),"",nPtBins,xBins);

		for(Int_t i_So = 0; i_So < nSobins; ++i_So ){ //changeH

			hSoTInNch[ibin][i_So] = 0;
			hSoTInNch[ibin][i_So]  = new TH1D(Form("hESTInNch%dES%d",ibin,i_So),"",1000,0.0,1.0);

			hSoMInNtrk[ibin][i_So] = 0;
			hSoMInNtrk[ibin][i_So]  = new TH1D(Form("hESMInNtrk%dES%d",ibin,i_So),"",1000,0.0,1.0);

			hSoTInNtrk[ibin][i_So] = 0;
			hSoTInNtrk[ibin][i_So]  = new TH1D(Form("hESTInNtrk%dES%d",ibin,i_So),"",1000,0.0,1.0);

			hNtrkInNtrkSoM[ibin][i_So] = 0;
			hNtrkInNtrkSoM[ibin][i_So] = new TH1D(Form("hNtrkInNtrkESM%dES%d",ibin,i_So),"",300,0,300);

			hNchInNchSoT[ibin][i_So] = 0;
			hNchInNchSoT[ibin][i_So] = new TH1D(Form("hNchInNchEST%dES%d",ibin,i_So),"",300,0,300);

			hMcOut[ibin][i_So] = 0;
			hMcOut[ibin][i_So] = new TH1D(Form("hMcOutMult%dES%d",ibin,i_So),"",nPtBins,xBins);

 			hMcOutnp[ibin][i_So] = 0;
                        hMcOutnp[ibin][i_So] = new TH1D(Form("hMcOutMultnp%dES%d",ibin,i_So),"",nPtBins,xBins);

			hPhiOut[ibin][i_So] = 0;
			hPhiOut[ibin][i_So] = new TH2D(Form("hMcPhiOutMult%dES%d",ibin,i_So),"",nPtBins,xBins,64,0,2.0*TMath::Pi());

			hMcSec[ibin][i_So] = 0;
			hMcSec[ibin][i_So] = new TH1D(Form("hMcSecMult%dES%d",ibin,i_So),"",nPtBins,xBins);

			hMcTrue[ibin][i_So] = 0;
			hMcTrue[ibin][i_So] = new TH1D(Form("hMcTrueMult%dES%d",ibin,i_So),"",nPtBins,xBins);

			hMcIn[ibin][i_So]  = 0;
			hMcIn[ibin][i_So]  = new TH1D(Form("hMcInMult%dES%d",ibin,i_So),"",nPtBins,xBins);

			hPhiIn[ibin][i_So] = 0;
			hPhiIn[ibin][i_So] = new TH2D(Form("hMcPhiInMult%dES%d",ibin,i_So),"",nPtBins,xBins,64,0,2.0*TMath::Pi());

		}


	}


	for(Int_t i_So = 0; i_So < nSobins; ++i_So ){

                hRMmultCombAllUnSom[i_So] = 0;
                hRMmultCombAllUnSom[i_So] = new TH2D(Form("hRMmultCombAllUnSom%d",i_So), "Combined estimator;#it{N}_{ch};#it{N}_{trk}",200, 0, 200, 200, 0, 200 );

                hRMmultCombAllUnSot[i_So] = 0;
                hRMmultCombAllUnSot[i_So] = new TH2D(Form("hRMmultCombAllUnSot%d",i_So), "Combined estimator;#it{N}_{ch};#it{N}_{trk}",200, 0, 200, 200, 0, 200 );

		hRMmultComb[i_So] = 0;
		hRMmultComb[i_So] = new TH2D( Form("hRMmultCombSom%d",i_So), "Combined estimator;d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
				nMultbins, Multbins, nMultbins, Multbins );

		hRMmultTPCo[i_So] = 0;
		hRMmultTPCo[i_So] = new TH2D( Form("hRMmultTPCoSom%d",i_So), "Track counting (TPC-only track cuts + TPC refit);d#it{N}_{ch}/d#eta;d#it{N}_{trk}/d#eta",
				nMultbins, Multbins, nMultbins, Multbins );


		hSotVsMult[i_So] = 0; // Sot: -1 to 1
		hMidEtaMultVsMultSot[i_So] = 0;

		hSomVsMult[i_So] = 0; // Sot: -1 to 1
		hMidEtaMultVsMultSom[i_So] = 0;

		if(isbinpercentile){

			hSotVsMult[i_So] = new TH2D( Form( "hSotVsMult%d", i_So ), "",
					nV0Mbins, V0Mbins, nSobins, Sobins );
			hMidEtaMultVsMultSot[i_So] = new TH2D( Form( "hMidEtaMultVsMultSot%d", i_So ), "", 
					nV0Mbins, V0Mbins, nMultbins, Multbins );

			hSomVsMult[i_So] = new TH2D( Form( "hSomVsMult%d", i_So ), "",
					nV0Mbins, V0Mbins, nSobins, Sobins );
			hMidEtaMultVsMultSom[i_So] = new TH2D( Form( "hMidEtaMultVsMultSom%d", i_So ), "",
					nV0Mbins, V0Mbins, nMultbins, Multbins );


		}else{

			hSotVsMult[i_So] = new TH2D( Form( "hSotVsMult%d", i_So), "",
					nMultbins, Multbins, nSobins, Sobins );  
			hMidEtaMultVsMultSot[i_So] = new TH2D( Form( "hMidEtaMultVsMultSot%d", i_So), "", 
					nMultbins, Multbins, nMultbins, Multbins );

			hSomVsMult[i_So] = new TH2D( Form( "hSomVsMult%d", i_So), "",
					nMultbins, Multbins, nSobins, Sobins );  
			hMidEtaMultVsMultSom[i_So] = new TH2D( Form( "hMidEtaMultVsMultSom%d", i_So), "",
					nMultbins, Multbins, nMultbins, Multbins );

		}

	}

	//Fake spherocity
	if(isbinpercentile){

		hSotVsMultFake = new TH2D( "hSotVsMultFake", "",
				nV0Mbins, V0Mbins, nSobins, Sobins );
		hMidEtaMultVsMultSotFake = new TH2D( "hMidEtaMultVsMultSotFake", "",
				nV0Mbins, V0Mbins, nMultbins, Multbins );

		hSomVsMultFake = new TH2D( "hSomVsMultFake", "",
				nV0Mbins, V0Mbins, nSobins, Sobins );
		hMidEtaMultVsMultSomFake = new TH2D( "hMidEtaMultVsMultSomFake", "",
				nV0Mbins, V0Mbins, nMultbins, Multbins );


	}else{

		hSotVsMultFake = new TH2D( "hSotVsMultFake", "",
				nMultbins, Multbins, nSobins, Sobins );
		hMidEtaMultVsMultSotFake = new TH2D( "hMidEtaMultVsMultSotFake", "",
				nMultbins, Multbins, nMultbins, Multbins );

		hSomVsMultFake = new TH2D( "hSomVsMultFake", "",
				nMultbins, Multbins, nSobins, Sobins );
		hMidEtaMultVsMultSomFake = new TH2D( "hMidEtaMultVsMultSomFake", "",
				nMultbins, Multbins, nMultbins, Multbins );

	}



}
//______________________________
TString CompleteDirName( TString nameOut3 ){

	TString nameOut = nameOut3;


	if(IsMC)
		nameOut += "Mc_";
	else
		nameOut += "Data_";
	if(evShape==0)
		nameOut += "So_";
	if(evShape==1)
		nameOut += "St_";

	if(etacriteria==0)
		nameOut += "AbsEta_";
	else if(etacriteria==1)
		nameOut += "PosEta_";
	else if(etacriteria==2)
		nameOut += "NegEta_";

	nameOut += estimator;
	if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==0) )
		nameOut += "Multcombined";
	if( (strstr(estimator, "CombTPCITS08")!= 0) && (multOjjects==1) )
		nameOut += "Multtpconly";
	nameOut += "_";
	nameOut += minMultSo;
	nameOut += "_ptMin";
	nameOut += ptminSo*1000;
	nameOut += "MeV";
	if(ptmaxSo<1.0e7){
		nameOut += "_ptMax";
		nameOut += ptmaxSo;
		nameOut += "GeV";
	}
	if(minMultSo>3){
		nameOut += "_minNtrks";
		nameOut += minMultSo;
	}
	if((type_track>0)&&(IsMC==kTRUE)){
		if(type_track==1)
			nameOut += "_prim";
		else if(type_track==2)
			nameOut += "_prim_wdecays";
	}
	if(isHybrid){
		nameOut += "_Hybrid";
	}
	if(isGolden){
		nameOut += "_Golden";
	}
	if(filter_observable!=1){
		nameOut += "_ObsFilter";
		nameOut += filter_observable;
	}

	return nameOut;
}
//_______________________________
Bool_t isGoodTrack( ESATrack* trackIn ){

	Bool_t isOk = kFALSE;

	Double_t phiI  = trackIn->phi;
	Double_t dcaXYI = trackIn->dcaxy;
	Double_t dcaZI = trackIn->dcaz;
	Int_t nclustersI = trackIn->nCl;
	Float_t chi2tpcI = trackIn->chi2perclustertpc;
	Bool_t kinkdaughtersI = trackIn->iskinkdaughter;
	Bool_t tpcrefitI = trackIn->istpcrefit;
	Bool_t itsrefitI = trackIn->isitsrefit;
	Bool_t itshitsI = trackIn->hasspdpoint;

	/* TPC only

	   Min. TPC Clusters         |     70
	   Max. chi^2 per TPC Cls.   |     4
	   Reject kink daughters     |    yes
	   TPC and ITS refit         |     no
	   Require hits in ITS       |     no
	   Max. DCA to vertex xy     |    2.4 cm.
	   Max. DCA to vertex z      |    3.2 cm.
	   */

	//same kin cuts as in the paper and analysis note

	if( nclustersI < 70 )
		return isOk;
	if( chi2tpcI > 4 )
		return isOk;
	if( kinkdaughtersI )
		return isOk;
	if( !(tpcrefitI) )
		return isOk;

	Float_t dcaToVertexI = CutDCAToVertex2D(trackIn,fCutMaxDCAToVertexXY,fCutMaxDCAToVertexZ,fCutDCAToVertex2D);
	if( fCutDCAToVertex2D && dcaToVertexI > 1 )
		return isOk;
	if( TMath::Abs(dcaXYI) > 2.4 )
		return isOk;
	if( TMath::Abs(dcaZI) > 3.2 )
		return isOk;

	return kTRUE;

}
//____________________________________________________________________
Float_t GetMultiplicity(ESAEvent* event){

	Float_t multiplicity = -10.0;

	if(strstr(estimator, "V0A") != 0) { multiplicity=event->v0Aperc; }
	else if(strstr(estimator, "V0M") != 0) { multiplicity=event->v0Mperc; }
	else if(strstr(estimator, "ADM") != 0) { multiplicity=event->adMperc; }
	else if(strstr(estimator, "CombTPCITS08") != 0) { multiplicity=event->trackmult08; }
	else if(strstr(estimator, "CombTPCITS05") != 0) { multiplicity=event->trackmult05; }

	return multiplicity;

}
//____________________________________________________________________
Float_t CutDCAToVertex2D(ESATrack* track,Double_t fCutMaxDCAToVertexXY,Double_t fCutMaxDCAToVertexZ, Bool_t fCutDCAToVertex2D){

	Float_t dcaToVertexXY = track->dcaxy;
	Float_t dcaToVertexZ = track->dcaz;

	Float_t dcaToVertex = -1;

	if (fCutDCAToVertex2D)
	{
		dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY/fCutMaxDCAToVertexXY/fCutMaxDCAToVertexXY + dcaToVertexZ*dcaToVertexZ/fCutMaxDCAToVertexZ/fCutMaxDCAToVertexZ);
	}
	else
		dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY + dcaToVertexZ*dcaToVertexZ);

	return dcaToVertex;


}
//______________________________________________
Int_t GetMultiplicityParticles(TClonesArray* TrackArray, Double_t etaCutlower, Double_t etaCutupper, Bool_t isMC){

	Int_t mult=-1;
	const Int_t nTracks = TrackArray->GetEntries();

	if(!isMC){  
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrack* track = (ESATrack*)TrackArray->At(i_a);

			if(!isGoodTrack(track))
				continue;

			Float_t eta_a = track->eta;

			//select leading particle in +/- etaCut

			//Added by Hector
			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			}
			//

			//if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
			//	continue;

			Float_t pt  = track->pt;

			if( pt < 0.15 )
				continue;

			mult++;


		}//end loop over particles
	}
	else{
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrackMC* track = (ESATrackMC*)TrackArray->At(i_a);

			//Float_t eta_a = track->etaMC;
			Float_t eta_a = track->etaMC;
			//select leading particle in +/- etaCut

			// added by Hector
			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			} 
			//

			//if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
			//	continue;

			Float_t pt  = track->ptMC;

			if( pt <= 0.0 )
				continue;

			//only charged
			if( TMath::Abs(track->qMC) < 0.1 )
				continue;      




			mult++;


		}//end loop over particles 
	}

	return mult;

}
//_______________________________________________________
Double_t GetSpherocity(TClonesArray* TrackArray, Double_t ptCutlower, Double_t ptCutupper, Double_t etaCutlower, Double_t etaCutupper, Int_t minMult, Bool_t isMC, Bool_t ispPb){

	Double_t spherocity=-10;
	Int_t mult=0;
	//  cout<<"empezando so"<<endl;
	const Int_t nTracks = TrackArray->GetEntries();

	if(!isMC){  
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrack* track = (ESATrack*)TrackArray->At(i_a);

			if(isHybrid){
				if(track->filter==0) 
					continue;
			}
			else if(isGolden){
				if(!(track->filter&1)) 
					continue;
			}
			else{
				if(!isGoodTrack(track))
					continue;
			}

			if( type_track > 0 && type_track < 3 ){

				if(type_track==1){
					if(track->primary!=type_track)
						continue;
				}
				if(type_track==2){
					if(track->primary==3)
						continue;
				}

			}

			Float_t eta_a = track->eta;
			if(ispPb){
				eta_a += 0.465;
			}

			//select leading particle in +/- etaCut

			//Added by Hector

			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			}			

			//


			Float_t pt  = track->pt;
			//if( pt < 0.15 )
			//	continue;

			//hPhi->Fill(track->phi);
			//hEta->Fill(eta_a);

			if( pt < ptCutlower || pt > ptCutupper )
				continue;

			mult++;


		}//end loop over particles
	}
	else{
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrackMC* track = (ESATrackMC*)TrackArray->At(i_a);

			//Float_t eta_a = track->etaMC;
			Float_t eta_a = track->etaMC;
			if(ispPb){
				eta_a += 0.465;
			}   
			//select leading particle in +/- etaCut

			//Added by Hector

			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			}
			//

			Float_t pt  = track->ptMC;

			if( pt < ptCutlower || pt > ptCutupper )
				continue;

			//only charged
			if( TMath::Abs(track->qMC) < 0.1 )
				continue;      

			mult++;


		}//end loop over particles 
	}


	if(mult<minMult)
		return -0.5;
	else{

		Double_t *pxA=new Double_t[mult];
		Double_t *pyA=new Double_t[mult];
		Double_t sumapt=0;
		Int_t counter=0;

		if(!isMC){
			for (Int_t i_a = 0; i_a < nTracks; ++i_a){


				ESATrack* track = (ESATrack*)TrackArray->At(i_a);

				if(isHybrid){
					if(track->filter==0)     
						continue;
				}
				else if(isGolden){
					if(!(track->filter&1))          
						continue;
				}
				else{   
					if(!isGoodTrack(track))
						continue;
				} 


				if( type_track > 0 && type_track < 3 ){

					if(type_track==1){
						if(track->primary!=type_track)
							continue;
					}
					if(type_track==2){
						if(track->primary==3)
							continue;
					}

				}

				//Float_t eta_a = track->eta;
				Float_t eta_a = track->eta;
				if(ispPb){
					eta_a += 0.465;
				}   
				//select leading particle in +/- etaCut
				// added by Hector

				if(etacriteria==0){
					if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
						continue;

				}

				if(etacriteria==1){
					if(eta_a<etaCutlower || eta_a>etaCutupper)
						continue;

				}
				if(etacriteria==2){
					if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
						continue;

				}
				//

				Float_t pt  = track->pt;
				Float_t phi = track->phi;    

				Double_t px=pt*TMath::Cos(phi);
				Double_t py=pt*TMath::Sin(phi);

				if( pt < ptCutlower || pt > ptCutupper )
					continue;

				pxA[counter]=px;
				pyA[counter]=py;
				sumapt+=pt;
				counter++;



			}//end loop over particles
		}
		else{
			for (Int_t i_a = 0; i_a < nTracks; ++i_a){


				ESATrackMC* track = (ESATrackMC*)TrackArray->At(i_a);

				//Float_t eta_a = track->etaMC;
				Float_t eta_a = track->etaMC;
				if(ispPb){
					eta_a += 0.465;
				}   
				//select leading particle in +/- etaCut
				//Added by Hector

				if(etacriteria==0){
					if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
						continue;

				}

				if(etacriteria==1){
					if(eta_a<etaCutlower || eta_a>etaCutupper)
						continue;

				}
				if(etacriteria==2){
					if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
						continue;

				}

				//

				Float_t pt  = track->ptMC;
				Float_t phi = track->phiMC;    

				//only charged
				if( TMath::Abs(track->qMC) < 0.1 )
					continue;      

				Double_t px=pt*TMath::Cos(phi);
				Double_t py=pt*TMath::Sin(phi);

				if( pt < ptCutlower || pt > ptCutupper )
					continue;

				pxA[counter]=px;
				pyA[counter]=py;
				sumapt+=pt;
				counter++;



			}//end loop over particles
		}

		Double_t pFull = 0;
		Double_t Spherocity = 2;
		//Getting spherocity
		for(Int_t i = 0; i < 360/(size_step); ++i){
			Double_t numerador = 0;
			Double_t phiparam  = 0;
			Double_t nx = 0;
			Double_t ny = 0;
			phiparam=((TMath::Pi()) * i * size_step ) / 180; // parametrization of the angle
			nx = TMath::Cos(phiparam);            // x component of an unitary vector n
			ny = TMath::Sin(phiparam);            // y component of an unitary vector n
			for(Int_t i1 = 0; i1 < mult; ++i1){
				numerador += TMath::Abs(ny * pxA[i1] - nx * pyA[i1]);//product between momentum proyection in XY plane and the unitari vector.
			}
			pFull=TMath::Power( (numerador / sumapt),2 );
			if(pFull < Spherocity)//maximization of pFull
			{
				Spherocity = pFull;
			}
		}

		spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;

		if(pxA){// clean up array memory used for TMath::Sort
			delete[] pxA; 
			pxA=0;
		}
		if(pyA){// clean up array memory used for TMath::Sort
			delete[] pyA; 
			pyA=0;
		}

	}
	//  cout<<"fin so"<<endl;

	return spherocity;
}
//______________________________________________________________________
Double_t GetSphericity(TClonesArray* TrackArray, Double_t ptCutlower, Double_t ptCutupper, Double_t etaCutlower, Double_t etaCutupper, Int_t minMult, Bool_t isMC, Bool_t ispPb){

	Double_t sphericity=-10;
	Double_t s00=0;
	Double_t s01=0;
	Double_t s11=0;
	Double_t totalpt=0;
	Int_t mult=0;
	const Int_t nTracks = TrackArray->GetEntries();

	if(!isMC){
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrack* track = (ESATrack*)TrackArray->At(i_a);

			if(isHybrid){
				if(track->filter==0)
					continue;
			}
			else if(isGolden){
				if(!(track->filter&1))
					continue;
			}
			else{
				if(!isGoodTrack(track))
					continue;
			}
			if( type_track > 0 && type_track < 3 ){

				if(type_track==1){
					if(track->primary!=type_track)
						continue;
				}
				if(type_track==2){
					if(track->primary==3)
						continue;
				}

			}

			Float_t eta_a = track->eta;
			if(ispPb)
				eta_a += 0.465;
			// Added by Hector

			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			}
			//


			Float_t pt  = track->pt;
			Float_t phi = track->phi;

			Double_t px=pt*TMath::Cos(phi);
			Double_t py=pt*TMath::Sin(phi);

			if( pt < ptCutlower || pt > ptCutupper )
				continue;

			s00 += (px * px)/pt;
			s01 += (py * px)/pt;
			s11 += (py * py)/pt;
			totalpt += pt;
			mult++;

		}
	}
	else{
		for (Int_t i_a = 0; i_a < nTracks; ++i_a){

			ESATrackMC* track = (ESATrackMC*)TrackArray->At(i_a);
			Float_t eta_a = track->etaMC;
			if(ispPb)
				eta_a += 0.465;

			// Added by Hector

			if(etacriteria==0){
				if( TMath::Abs(eta_a) < etaCutlower || TMath::Abs(eta_a) > etaCutupper )
					continue;

			}

			if(etacriteria==1){
				if(eta_a<etaCutlower || eta_a>etaCutupper)
					continue;

			}
			if(etacriteria==2){
				if(eta_a>etaCutlower || eta_a<(-1*etaCutupper))
					continue;

			}
			//

			Float_t pt  = track->ptMC;
			Float_t phi = track->phiMC;

			Double_t px=pt*TMath::Cos(phi);
			Double_t py=pt*TMath::Sin(phi);

			if( pt < ptCutlower || pt > ptCutupper )
				continue;

			if( TMath::Abs(track->qMC) < 0.1 )
				continue;

			s00 += (px * px)/pt;
			s01 += (py * px)/pt;
			s11 += (py * py)/pt;
			totalpt += pt;
			mult++;

		}
	}
	if(mult<minMult){
		sphericity = -0.5;
	}
	else{

		Double_t S00=s00/totalpt;
		Double_t S01=s01/totalpt;
		Double_t S11=s11/totalpt;

		Float_t lambda1=((S00+S11)+TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
		Float_t lambda2=((S00+S11)-TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
		if((lambda2==0)&&(lambda1==0))
			sphericity=0;
		if(lambda1+lambda2!=0)
			sphericity=2*TMath::Min( lambda1,lambda2 )/( lambda1+lambda2 );
	}

	return sphericity;
}
//____________________________________________________________
Int_t GetMeasMultiplicityBin(Double_t multIn){

	Int_t binOut = -1;

	if( strstr(estimator, "CombTPCITS08")!= 0 ){
		for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){

			if( (multIn >= Multbins[ibin]) && (multIn < Multbins[ibin+1]) )
				binOut = ibin;

		}
	}
	else if( strstr(estimator, "V0M")!= 0 ){
		for( Int_t ibin = 0; ibin < nV0Mbins; ++ibin ){

			if( (multIn >= V0Mbins[ibin]) && (multIn < V0Mbins[ibin+1]) )
				binOut = ibin;

		}
	}

	return binOut;
}

//____________________________________________________________
Int_t GetMeasMultiplicityBin2(Double_t multIn){

        Int_t binOut = -1;

        if( strstr(estimator, "CombTPCITS08")!= 0 ){
                for( Int_t ibin = 0; ibin < nMultbins2; ++ibin ){

                        if( (multIn >= Multbins2[ibin]) && (multIn < Multbins2[ibin+1]) )
                                binOut = ibin;

                }
        }
        else if( strstr(estimator, "V0M")!= 0 ){
                for( Int_t ibin = 0; ibin < nV0Mbins; ++ibin ){

                        if( (multIn >= V0Mbins[ibin]) && (multIn < V0Mbins[ibin+1]) )
                                binOut = ibin;

                }
        }

        return binOut;
}

