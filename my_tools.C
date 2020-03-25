#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

using namespace std;



//_____________________________________________________________________________


TFile* FindFileFresh(const Char_t* fileName);
TFile* FindFile(const Char_t* fileName);
void CutHistogram(TH1* hist, Double_t xMin, Double_t xMax);
void SetHistError(TH1* hist, Double_t error);
void CreateDir(const Char_t* dirName);

//______________________________________________________________________
TFile* FindFileFresh(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

//______________________________________________________________________
TFile* FindFile(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    return file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

//______________________________________________________________________
void CutHistogram(TH1* hist, Double_t xMin, Double_t xMax)
{
  const Int_t n = hist->GetNbinsX();
  
  for(Int_t bin = 1; bin <= n; bin++) {
    
    Float_t x = hist->GetXaxis()->GetBinCenter(bin);
    if(x < xMin) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    } else if(x > xMax) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    }

  }
}

//______________________________________________________________________
void SetHistError(TH1* hist, Double_t error)
{
  const Int_t n = hist->GetNbinsX();
  
  for(Int_t bin = 1; bin <= n; bin++) {
    
    //    Float_t x = hist->GetXaxis()->GetBinCenter(bin);
    hist->SetBinError(bin, error);
  }
}

//______________________________________________________________________
void CreateDir(const Char_t* dirName)
{
  TString pwd(gSystem->pwd());
  gSystem->cd(pwd.Data());
  
  if(gSystem->cd(dirName)) {
    gSystem->cd(pwd.Data());
  } else {
    gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
  }
}
