/*
 * This macro find bad runs and draw comparison of distributions 
 * until and after RanQA. One can use root file that was 
 * obtain after work FemtoDstQA example.
 *
 * authors: Alexey Povarov <povarovas@gmail.com>
 * date September 22, 2019.
 */

#include <iostream>
#include <math.h>
#include <vector>
#include <TProfile.h>
#include <TGraph.h>

Float_t RoundValue( Double_t value );

std::vector<int> GetBadRuns( TFile *f, const Char_t *tprofile, std::vector<Int_t> BadRuns,
									   const Char_t *inFileRunQA, const Char_t *energy, 
									   const Char_t *pathPics, const Char_t *format);


void FindBadRuns(const Char_t *inFileNoRunQA = "QAtest14gev.root",
			 	 const Char_t *inFileRunQA = "",
			     const Char_t *energy = "14gev",
			     const Char_t *pathPics = "",
			     const Char_t *format = "") {

	TFile *f1 = new TFile(inFileNoRunQA,"READ");
	std::vector<Int_t> BadRunList;
	
	for(Int_t i = 0; i < 10; i++) {
		BadRunList = GetBadRuns(f1, Form("hEventProfile_%i",i), BadRunList, inFileRunQA,energy,pathPics,format);
		if( i < 6 ) {
			BadRunList = GetBadRuns(f1, Form("hTrackProfile_%i",i), BadRunList, inFileRunQA,energy,pathPics,format);
		}
		if( i < 3) {
			BadRunList = GetBadRuns(f1, Form("hSinPhi%i",i+1), BadRunList, inFileRunQA,energy,pathPics,format);
			BadRunList = GetBadRuns(f1, Form("hCosPhi%i",i+1), BadRunList, inFileRunQA,energy,pathPics,format);
		}
	}

	if( strncmp(inFileRunQA,"",1) == 0 ) {

		std::sort(BadRunList.begin(), BadRunList.end() );
		BadRunList.erase( unique( BadRunList.begin(), BadRunList.end() ), BadRunList.end() );

		for(Int_t i = 0; i < BadRunList.size(); i++) {
			if( i != 0 && i%5 == 0 ) cout << "\n";
			cout << BadRunList[i];
			if( i < BadRunList.size() - 1) cout << ",";
		}
		cout << endl;
	}
	

}

Float_t RoundValue( Double_t value ) {
	
	if( abs(value) < 1.0 ) {

		Int_t n = 1;
		while( abs(value) < 10.0 ) {
			value*=10.0;
			n*=10;
		};

		value = ceil(value)/(Float_t)n;
	}
	if( abs(value) >= 1.0 ) value = ceil(value*100.0)/100.0;
	

	return value;
}

std::vector<int> GetBadRuns(TFile *f, const Char_t *tprofile, std::vector<Int_t> BadRuns, 
									  const Char_t *inFileRunQA, const Char_t *energy, 
									  const Char_t *pathPics, const Char_t *format) {
	cout << "Working on " << tprofile << endl;
	Int_t runIdBins;
  	Int_t runIdRange[2];

  	if( strncmp(energy, "39gev",5) == 0) {
	    runIdRange[0] = 11095000;
	    runIdRange[1] = 11115000;
	    runIdBins = runIdRange[1] - runIdRange[0];
  	}// if( strncmp(energy, "39GeV",5) == 0 )

  	if( strncmp(energy, "27gev",5) == 0 ){
	    runIdRange[0] = 12171000;
	    runIdRange[1] = 12180000;
	    runIdBins = runIdRange[1] - runIdRange[0];
	}// if( strncmp( energy, "27GeV",5) == 0 )

	if( strncmp(energy, "19GeV",5) == 0) {
	    runIdRange[0] = 12110000;
	    runIdRange[1] = 12123000;
	    runIdBins = runIdRange[1] - runIdRange[0];
  	}// if( strncmp(energy, "19GeV",5) == 0 )

	if( strncmp(energy, "14gev",5) == 0 ) {
	    runIdRange[0] = 15045000;
	    runIdRange[1] = 15075000;
	    runIdBins = runIdRange[1] - runIdRange[0];
  	}// if( strncmp(energy, "14GeV",5) == 0)

  	if( strncmp(energy, "11GeV",5) == 0 ) {
	    runIdRange[0] = 11145000;
	    runIdRange[1] = 11165000;
	    runIdBins = runIdRange[1] - runIdRange[0];
    }// if( strncmp(energy, "11GeV",5) == 0 )

    if(strncmp(energy, "7GeV",4)==0){
	    runIdRange[0] = 11110000;
	    runIdRange[1] = 11150000;
	    runIdBins = runIdRange[1] - runIdRange[0];
  	}// if(strncmp(energy, "7GeV",4)==0)

  TProfile *tp = (TProfile*)f -> Get( tprofile );
	std::vector<Float_t> VbinContent , VbinError;
	Float_t binContent = 0, binError = 0;
	Float_t Mcontent = 0, Merror = 0, Scontent = 0, Serror = 0;
	Float_t Sc = 0, Serr = 0, sigmaC, sigmaErr; 
	Int_t n = 0, k = 0;

	for( Int_t bin = 1; bin <= runIdBins; bin++ ) { //get content and error from TProfile and calc mean
		binContent  = tp -> GetBinContent(bin);
		binError = tp -> GetBinError(bin);
		if( binContent  == 0 || binError == 0 ) continue;
		VbinContent.push_back(binContent);
		VbinError.push_back(binError);
		Scontent  += binContent ;
		Serror += binError; 
		n++;
	}
	Mcontent  = Scontent/n;
	Merror = Serror/n;

	for( Int_t i = 0; i < VbinContent.size(); i++ ) { // root mean square
		Sc += VbinContent[i]*VbinContent[i] - Mcontent*Mcontent;
		Serr += VbinError[i]*VbinError[i] - Merror*Merror;
	} 

	Sc = Sc/(n-1);
	Serr = Serr/(n-1);
	sigmaC = sqrt(Sc);     //sigma content 
	sigmaErr = sqrt(Serr); //sigma error

	k = 0;
	Int_t b = 0;
	for( Int_t bin = 1; bin <= runIdBins; bin++ ) { //bad nambers Run on desktop 
		binContent  = tp -> GetBinContent(bin);
		binError = tp -> GetBinError(bin);
		k++;
		if( binContent  == 0 || binError == 0 ) continue;
		if ( TMath::Abs( binContent - Mcontent ) > 3*sigmaC || TMath::Abs( binError - Merror ) > 3*sigmaErr ) {
			BadRuns.push_back( bin + runIdRange[0] - 1 );
			b++;
		}
	}
	
	Float_t badRun[b];             // arrays for graphics
	Float_t badContent[b];
	Float_t badCErr[b], badRErr[b];

	k =0;
	b = 0;
	for( Int_t bin = 1; bin <= runIdBins; bin++ ) {       //record Bad Run in array for graph 
		binContent  = tp -> GetBinContent(bin);
		binError = tp -> GetBinError(bin);
		k++;
		if( binContent  == 0 || binError == 0 ) continue;
		if ( TMath::Abs( binContent - Mcontent ) > 3*sigmaC || TMath::Abs( binError - Merror ) > 3*sigmaErr ) {
			badRun[b] = bin + runIdRange[0] - 1;
			badContent[b] = binContent;
			badCErr[b] = binError;
			badRErr[b] = 0;
			b++;
		}
	}

	if( strncmp(inFileRunQA,"",1) != 0 ) {

		TCanvas *c1 = new TCanvas("c1", "c1", 1366, 768);
  		c1->Divide(2);
  		c1 -> cd(1);

		tp -> SetStats(0);
		tp -> SetMarkerStyle(20); 
		tp -> SetMarkerColor(kBlack); 
		tp -> SetAxisRange(Mcontent - 15*sigmaC,Mcontent + 15*sigmaC,"Y");
		tp -> Draw();

		TLegend *legC2_1 = new TLegend(0.54,0.7,0.89,.89);
	    legC2_1->AddEntry(tp,Form("Mean = %f",RoundValue(Mcontent) ),"");
	    legC2_1->AddEntry(tp,Form("Mean error = %f",RoundValue(Merror) ),"");    
	    legC2_1->AddEntry(tp,Form("Sigma = %f", RoundValue(sigmaC) ),"");
	    legC2_1->AddEntry(tp,Form("Sigma error = %f", RoundValue(sigmaErr) ),"");
	    legC2_1->SetFillColor(kWhite);
	    legC2_1->SetBorderSize(0); 
	    legC2_1->Draw();
		
		Float_t x[2], mean[2]; 
		x[0] = runIdRange[0] - 1000, x[1] = runIdRange[1] + 1000;
		mean[0] = Mcontent, mean[1] = Mcontent;
		TGraph *gmean = new TGraph(2,x,mean);  //paint calculated mean
		gmean -> SetLineColor(kRed);
		gmean -> SetLineWidth(2);
		gmean -> Draw("CPsame");
		
		Float_t sigmapos[2];
		sigmapos[0] = Mcontent + 3*sigmaC, sigmapos[1] = Mcontent + 3*sigmaC;
		TGraph *gsigmapos = new TGraph(2,x,sigmapos);  //paint top sigma 
		gsigmapos -> SetLineColor(kRed);
		gsigmapos -> SetLineStyle(9);
		gsigmapos -> SetLineWidth(2);
		gsigmapos -> Draw("CPsame");
		
		Float_t sigmaneg[2];
		sigmaneg[0] = Mcontent - 3*sigmaC, sigmaneg[1] = Mcontent - 3*sigmaC;
		TGraph *gsigmaneg = new TGraph(2,x,sigmaneg);  //paint down sigma 
		gsigmaneg -> SetLineColor(kRed);
		gsigmaneg -> SetLineStyle(9);
		gsigmaneg-> SetLineWidth(2);
		gsigmaneg -> Draw("CPsame");
		
		TGraphErrors* gbaddot = new TGraphErrors(b,badRun,badContent,badRErr,badCErr); //paint red dot for visualisation
		gbaddot -> SetMarkerStyle(20); 
		gbaddot -> SetLineWidth(2);
		gbaddot -> SetLineColor(kRed);
		gbaddot -> SetMarkerColor(kRed); 
		gbaddot -> Draw("Psame");

		TFile *f2 = new TFile(inFileRunQA,"READ");
		c1 -> cd(2);
		TProfile *tp1 = (TProfile*) f2 -> Get( tprofile );
		tp1 -> SetStats(0);
		tp1 -> SetMarkerStyle(20); 
		tp1 -> SetMarkerColor(kBlack); 
		tp1 -> SetAxisRange(Mcontent - 15*sigmaC,Mcontent + 15*sigmaC,"Y");
		tp1 -> Draw();
		gmean -> Draw("CPsame");
		gsigmaneg -> Draw("CPsame");
		gsigmapos -> Draw("CPsame");
		c1->SaveAs( Form("%s%s.%s",pathPics,tprofile,format) );
	}
	
	cout << "Done!" << endl; 
	return BadRuns;
}

