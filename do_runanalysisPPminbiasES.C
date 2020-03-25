void do_runanalysisPPminbiasES(Char_t* rootname = "", Int_t nmaxfiles, Int_t startidx)
{

	gSystem->Exec("rm -f *.so");
	gSystem->Exec("rm -f *.d");

	gSystem->AddIncludePath("-I$ALICE_ROOT/include");
	gSystem->AddIncludePath("-I$ALICE/aliroot/git/TPC/Base"); 
        gROOT->ProcessLine(".L my_tools.C+");
	gROOT->ProcessLine(".L DebugClassesMultESA2016.C+");
	gROOT->ProcessLine(".L analysisPPminbiasSo.C+");

	TString cmd1 = Form("GetSpherocitySpectrum(\"%s\",\"SOpp13TeV_spec.root\",0,%d,%d)",rootname,nmaxfiles,startidx);
	//TString cmd1 = Form("GetHistosForEqual(\"%s\",\"SOpp13TeV.root\",0,%d,%d)",rootname,nmaxfiles,startidx);

	//TString cmd1 = Form("GetSpherocitySpectrum(\"test.list\",\"TestSo.root\",0,20,0)");
	//TString cmd1 = Form("GetHistosForEqual(\"test.list\",\"TestSo.root\",0,20,0)");

	gROOT->ProcessLine(cmd1);

}

