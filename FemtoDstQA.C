/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading STAR FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \authors: Grigory Nigmatkulov, Povarov Alexey, Demanov Alexandr
 * \date May 29, 2018
 */

// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream> 

// ROOT headers
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"

// FemtoDst headers

#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoDstReader.h"
#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoDst.h"
#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoEvent.h"
#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoTrack.h"
#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoV0.h"
#include "/mnt/pool/1/aspovarov/basov/StFemtoEvent/StFemtoXi.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(/mnt/pool/1/aspovarov/basov/StFemtoEvent/libStFemtoDst)
#endif

// Forward declarations
// Check event and track
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZcut, Double_t VtxRcut, Double_t VtxYdelta );
Bool_t isGoodTrack(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCAcut);

const Float_t electron_mass = 0.0005485799;
const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Float_t electron_mass_sqr = 0.000000301;
const Float_t pion_mass_sqr = 0.019479955;
const Float_t kaon_mass_sqr = 0.24371698;
const Float_t proton_mass_sqr = 0.880354499;

Double_t mZDCSMDCenterex = 4.72466;
Double_t mZDCSMDCenterey = 5.53629;
Double_t mZDCSMDCenterwx = 4.39604;
Double_t mZDCSMDCenterwy = 5.19968;



// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.FemtoDst.root files
//_________________
void FemtoDstQA(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                const Char_t *outFileName = "oTest.root",
                const Char_t *energy = "14GeV") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  gSystem->Load("/mnt/pool/1/aspovarov/basov/StFemtoEvent/libStFemtoDst.so");
  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  gSystem->Load("/mnt/pool/1/aspovarov/basov/StFemtoEvent/libStFemtoDst.so");
  #endif

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  //Long64_t events2read = femtoReader->chain()->GetEntries();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInt_tree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInt_tree: "  << eventsInt_tree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  Int_t runIdBins;
  Int_t runIdRange[2];

  Double_t VtxZcut;
  Double_t VtxRcut;
  Double_t VtxYdelta;
  Double_t DCAcut;

  std::vector<Int_t> badRuns;
  Bool_t useRunQA;
  Bool_t CUTS;
  //false - RunQA OFF
  //true - RunQA ON

  Int_t c = 9;    // cent9()
  Int_t cent, RunID;
  Double_t Phi, pt, q, SqM;

  if( strncmp(energy, "39GeV",5) == 0) {

    runIdRange[0] = 11095000;
    runIdRange[1] = 11115000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut = 40.0;
    VtxRcut = 2.0;
    VtxYdelta = 0.0;
    DCAcut = 2.0;
  
    CUTS = false;
    useRunQA = false;
    badRuns ={11099102, 11099103, 11099104, 11099106, 11099107, 
    		      11099125, 11100004, 11100005, 11100008, 11100010, 
    		      11100011, 11100016, 11100020, 11100071, 11101014, 
    		      11101104, 11102098, 11103009, 11103047, 11103065, 
    		      11105011, 11105029, 11106026, 11106027, 11106028, 
    		      11106029, 11106030, 11106040, 11106041, 11107008, 
    		      11107046, 11107083, 11108040, 11108053, 11108065, 
    		      11108075, 11109092, 11109102, 11109105, 11109104, 
    		      11110005, 11110041, 11110042, 11110086};

  }// if( strncmp(energy, "39GeV",5) == 0 )

  if( strncmp(energy, "27GeV",5) == 0 ){

    runIdRange[0] = 12171000;
    runIdRange[1] = 12180000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut = 70.0;
    VtxRcut = 2.0;
    VtxYdelta = 0.0;
    DCAcut = 2.0;
    
    CUTS = false;
    useRunQA = false;
    badRuns ={12172049, 12172056, 12173009, 12173018, 12173026, 
    		      12173053, 12173054, 12173055, 12173056, 12173057, 
    		      12173072, 12174077, 12174096, 12174109, 12175007, 
    		      12175030, 12175089, 12176046, 12176047, 12176067, 
    		      12176069, 12176104, 12178051, 12178093, 12179068, 
    		      12179083, 12179084, 12179085, 12179086};

  }// if( strncmp( energy, "27GeV",5) == 0 )

  if( strncmp(energy, "19GeV",5) == 0) {

    runIdRange[0] = 12110000;
    runIdRange[1] = 12123000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut = 70.0;
    VtxRcut = 2.0;
    VtxYdelta = 0.0;
    DCAcut = 2.0;
    
    CUTS = false;
    useRunQA = false;
    badRuns ={12113081, 12113084, 12114077, 12114079, 12114085,
              12114088, 12114089, 12114091, 12114094, 12114095,
              12114097, 12114098, 12114099, 12114100, 12114101,
              12114102, 12114103, 12114104, 12114110, 12115025,
              12115026, 12116015, 12116016, 12116063, 12116084,
              12117009, 12118045, 12119008, 12119009, 12119011,
              12119015, 12119016, 12119017, 12119019, 12119020,
              12119021, 12119022, 12119023, 12119024, 12119025,
              12119027, 12119028, 12119029, 12119030, 12119032,
              12119035, 12119036, 12119039, 12120018, 12120073};

  }// if( strncmp(energy, "19GeV",5) == 0 )

  if( strncmp(energy, "14GeV",5) == 0 ) {

    runIdRange[0] = 15045000;
    runIdRange[1] = 15075000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut = 70.0;
    VtxRcut = 1.0;
    VtxYdelta = 0.8847;
    DCAcut = 2.0;
    
    CUTS = false; 
    useRunQA = false;
    badRuns ={15053027, 15053028, 15053029, 15053034, 15053035,
              15053048, 15053052, 15053053, 15053054, 15053055,
              15053056, 15053057, 15054019, 15054053, 15054054,
              15055018, 15055131, 15055133, 15055134, 15055135,
              15055136, 15055137, 15055138, 15055139, 15055140,
              15055141, 15056001, 15056002, 15056003, 15056004,
              15056005, 15056006, 15056007, 15056008, 15056009,
              15056113, 15056114, 15056116, 15056117, 15056124,
              15056125, 15057001, 15057003, 15057004, 15057006,
              15057007, 15057010, 15057011, 15057013, 15057014,
              15057018, 15057055, 15057059, 15058006, 15058011,
              15058018, 15060061, 15060069, 15061001, 15061002,
              15062006, 15065012, 15065014, 15066008, 15066013,
              15066017, 15068013, 15068014, 15068016, 15069034,
              15069036, 15070009, 15070010};

  }// if( strncmp(energy, "14GeV",5) == 0)

  if( strncmp(energy, "11GeV",5) == 0 ) {

    runIdRange[0] = 11145000;
    runIdRange[1] = 11165000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut=50.0;
    VtxRcut=2.0;
    VtxYdelta =0.0;
    DCAcut=2.0;

    CUTS = false;
    useRunQA = false;
    badRuns ={11148001, 11148008, 11148009, 11148010, 11148036,
              11148055, 11149017, 11149018, 11149040, 11149043,
              11150017, 11150029, 11151051, 11151057, 11153045,
              11154026, 11154040, 11154059, 11156036, 11156043,
              11156044, 11156045, 11157039};

  }// if( strncmp(energy, "11GeV",5) == 0 )

  if(strncmp(energy, "7GeV",4)==0){

    runIdRange[0] = 11110000;
    runIdRange[1] = 11150000;
    runIdBins = runIdRange[1] - runIdRange[0];

    VtxZcut=50.0;
    VtxRcut=2.0;
    VtxYdelta =0.0;
    DCAcut=2.0;

  }// if(strncmp(energy, "7GeV",4)==0)

  TFile *outFile = new TFile(outFileName, "RECREATE");

  // Histogramming
  // Event
  TH1F *hNeventsVsRunId = new TH1F("hNeventsVsRunId", "<Number of events> vs Run number;RunId;<N_{events}>",
                                    runIdBins, runIdRange[0], runIdRange[1]);
  Int_t mColor = 1;
  hNeventsVsRunId->SetMarkerStyle(20);    // filled circle
  hNeventsVsRunId->SetMarkerColor(mColor);
  hNeventsVsRunId->SetMarkerSize(1.1);
  hNeventsVsRunId->SetLineWidth(2);
  hNeventsVsRunId->SetLineColor(mColor);

  // Reference multiplicity histograms
  // 1D
  TH1D *hRefMult = new TH1D("hRefMult", "Reference multiplicity;RefMult;Entries",
                            600, -0.5, 599.5);
  TH1D *hRefMult2 = new TH1D("hRefMult2","Reference multiplicity in |#eta|<1;RefMult2;Entries",
                             600, -0.5, 599.5);
  TH1D *hGRefMult = new TH1D("hGRefMult","Reference multiplicity of global tracks;gRefMult;Entries",
                             800, -0.5, 799.5);
  // 2D
  TH2D *hRefMultVsAdcZdcE = new TH2D("hRefMultVsAdcZdcE","Reference multiplicity vs Adc_{ZDCe};N_{RefMult};Adc_{ZDCe}",
                                      600, -0.5, 599.5,400,0.,4000.);
  TH2D *hRefMultVsAdcZdcW = new TH2D("hRefMultVsAdcZdcW","Reference multiplicity vs Adc_{ZDCw};N_{RefMult};Adc_{ZDCw}",
                                      600, -0.5, 599.5,400,0.,4000.);
  TH2D *hRefMultVsAdcBbcE = new TH2D("hRefMultVsAdcBbcE","Reference multiplicity vs Adc_{BBCe};N_{RefMult};Adc_{BBCe}",
                                      600, -0.5, 599.5,2500,0.,50000.);
  TH2D *hRefMultVsAdcBbcW = new TH2D("hRefMultVsAdcBbcW","Reference multiplicity vs Adc_{BBCw};N_{RefMult};Adc_{BBCw}",
                                      600, -0.5, 599.5,2500,0.,50000.);


  // ZDC and BBC hist
  // 1D 
  TH1D *hAdcZdcEast = new TH1D("hAdcZdcEast","AdcSum ZDC East;AdcSum_{ZDCe};Counts", 
                                400,0.,4000.);
  TH1D *hAdcZdcWest = new TH1D("hAdcZdcWest","AdcSum ZDC West;AdcSum_{ZDCw};Counts", 
                                400,0.,4000.);
  TH1D *hAdcZdcSum = new TH1D("hAdcZdcSum","AdcSum ZDC;AdcSum_{ZDC};Counts", 
                                400,0.,4000.);
  TH1D *hAdcBbcEast = new TH1D("hAdcBbcEast","AdcSum BBC East;AdcSum_{BBCe};Counts", 
                                2000,0.,40000.);
  TH1D *hAdcBbcWest = new TH1D("hAdcBbcWest","AdcSum BBC West;AdcSum_{BBCw};Counts", 
                                2000,0.,40000.);
  TH1D *hAdcBbcSum = new TH1D("hAdcBbcSum","AdcSum BBC;AdcSum_{BBC};Counts", 
                                2*2000,0.,2*40000.);
  // 2D
  TH2D *hAdcZdcEvsW = new TH2D("hAdcZdcEvsW","AdcSum ZDC East Vs AdcSum ZDC West;AdcSum_{ZDCe};AdcSum_{ZDCw}",
                                400,0.,4000.,400,0.,4000.);
  TH2D *hAdcBbcEvsW = new TH2D("hAdcBbcEvsW","AdcSum BBC East Vs AdcSum BBC West;AdcSum_{BBCe};AdcSum_{BBCw}",
                                2500,0.,50000.,2500,0.,50000.);
  TH2D *hAdcBbcEastVsTile = new TH2D("hAdcBbcEastVsTile","AdcSum BBC East vs Number tile;N_{tile};AdcSum_{BBCe}",
                                      25,0.,25.,4.*400,0.,8.2*4000.);
  TH2D *hAdcBbcWestVsTile = new TH2D("hAdcBbcWestVsTile","AdcSum BBC West vs Number tile;N_{tile};AdcSum_{BBCw}",
                                      25,0.,25.,4.*400,0.,8.2*4000.);

  // Primary Vertex histogram
  // 1D
  TH1D *hVtxX = new TH1D("hVtxX","hVtxX;y [cm]; Entries",
                         100,-3.5,3.5);
  TH1D *hVtxY = new TH1D("hVtxY","hVtxY;y [cm]; Entries",
                         100,-3.5,3.5);
  TH1D *hVtxZ = new TH1D("hVtxZ","hVtxZ;z [cm]; Entries",
                         200, -100., 100.);
  // 2D
  TH2D *hVtxXvsRefMult = new TH2D("hVtxXvsRefMult","hVtxXvsRefMult;x [cm]; RefMult",
                            70,-3.5,3.5,600, -0.5, 599.5);
  TH2D *hVtxYvsRefMult = new TH2D("hVtxYvsRefMult","hVtxYvsRefMult;y [cm]; RefMult",
                            70,-3.5,3.5,600, -0.5, 599.5);
  TH2D *hVtxZvsRefMult = new TH2D("hVtxZvsRefMult","hVtxZvsRefMult;z [cm]; RefMult",
                            200,-100.,100.,600, -0.5, 599.5);
  TH2D *hVtxXvsY = new TH2D("hVtxXvsY", "hVtxXvsY;x [cm];y [cm]",
                            140,-3.5,3.5,140,-3.5,3.5);
  TH2D *hVtxZvsX = new TH2D("hVtxZvsX", "hVtxZvsX;z [cm];x [cm]",
                            200, -100., 100.,140,-3.5,3.5);
  TH2D *hVtxZvsY = new TH2D("hVtxZvsY", "hVtxZvsY;z [cm];y [cm]",
                            200, -100., 100.,140,-3.5,3.5);
  TH2D *hVpdVzVsVtxZ = new TH2D("hVpdVzVsVtxZ","VpdVz Vs Vz Vertex position;z_{VPD} [cm];z [cm]",
                                  400.,-100.,100.,400.,-100.,100.);
  TH2D *hVpdVzDiffVsVz = new TH2D("hVpdVzDiffVsVz","v_{z}(TPC) - v_{z}(VPD) vs. v_{z}(TPC);v_{z}(TPC);v_{z}(TPC) - v_{z}(VPD)",
                                  280, -70., 70., 80, -20., 20.);

  TH1D *hNumberOfPrimaries = new TH1D("hNumberOfPrimaries","Number of primary tracks;Number of primary tracks;Entries",
                                      1000, -0.5, 999.5);
  TH1D *hNumberOfGlobals = new TH1D("hNumberOfGlobals","Number of global tracks;Number of global tracks;Entries",
                                    1500, -0.5, 1499.5);
  TH1D *hCent9 = new TH1D("hCent9","Centralitity;Cent9;Entries",
                          13, -1.5, 11.5);
  TH1D *hCent16 = new TH1D("hCent16","Centralitity;Cent16;Entries",
                          19, -1.5, 17.5);
  TH1D *hBTofHit = new TH1D("hBTofHit","Number of hits in TOF;bTofTrayMult;Entries",
                            800, -0.5, 799.5);
  TH1D *hBTofMatched = new TH1D("hBTofMatched","Number of TOF-matched tracks;bTofMatched;Entries",
                                400, -0.5, 399.5);
  TH1D *hBemcMatched = new TH1D("hBemcMatched","Number of BEMC-matched tracks;bEmcMatched;Entries",
                                400, -0.5, 399.5);
  TH1D *hRanking = new TH1D("hRanking","Primary vertex ranking;Primary vertex ranking;Entries",
                            21, -10.5, 10.5);
  
  TH2D *hBTofTrayMultVsRefMult = new TH2D("hBTofTrayMultVsRefMult","TOF tray multiplicity vs. refMult;refMult;bTofTrayMult",
                                          600, -0.5, 599.5, 1500, -0.5, 1499.5);
  TH2D *hBTofMatchedVsRefMult = new TH2D("hBTofMatchedVsRefMult","TOF-matched tracks vs. refMult;refMult;TOF-matched",
                                          600, -0.5, 599.5, 400, -0.5, 399.5);
  TH1D *hTransSphericity = new TH1D("hTransSphericity","Transverse sphericity;Sphericity;Entries",
                                    10, 0., 1.);
  TH1D *hTransSphericity2 = new TH1D("hTransSphericity2","Transverse sphericity in |#eta|<1;Sphericity;Entries",
                                     10, 0., 1.);
  TH1D *hNumberOfVertices = new TH1D("hNumberOfVertices","Number of primary vertices;Number of primary vertices;Entries",
                                     15, -0.5, 14.5);
  TProfile *hEventProfile[10];
  hEventProfile[0] = new TProfile("hEventProfile_0","Profile of refMult;Run ID;<refMult>",
                                  runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[1] = new TProfile("hEventProfile_1","Profile of TOF tray multiplicity;Run ID;<bTofTrayMultiplicity>",
                                  runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[2] = new TProfile("hEventProfile_2","Profile of TOF-matched tracks;Run ID;<bTofMatched>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[3] = new TProfile("hEventProfile_3","Profile of number of primary tracks;Run ID;<nPrimTracks>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[4] = new TProfile("hEventProfile_4","Profile of number of global tracks;Run ID;<nGlobTracks>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[5] = new TProfile("hEventProfile_5","Profile of ZDC ADC;Run ID; <ADC_{ZDC}>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[6] = new TProfile("hEventProfile_6","Profile of BBC ADC;Run ID; <ADC_{BBC}>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[7] = new TProfile("hEventProfile_7","Profile of primary vertex X position;Run ID; <VtxX> [cm]",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[8] = new TProfile("hEventProfile_8","Profile of primary vertex Y position;Run ID; <VtxY> [cm]",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  hEventProfile[9] = new TProfile("hEventProfile_9","Profile of primary vertex Z position;Run ID; <VtxZ> [cm]",
                                            runIdBins, runIdRange[0], runIdRange[1] );
  for(int iHist=0; iHist<10; iHist++) {
    Int_t mColor = 1;
    hEventProfile[iHist]->SetMarkerStyle(20);    // filled circle
    hEventProfile[iHist]->SetMarkerColor(mColor);
    hEventProfile[iHist]->SetMarkerSize(1.1);
    hEventProfile[iHist]->SetLineWidth(2);
    hEventProfile[iHist]->SetLineColor(mColor);  // black
  }

  

  // Track
  // Momentum histogram
  TH1D *hGlobalPtot = new TH1D("hGlobalPtot","Global track momentum;p (GeV/c);Entries",
                               200, 0., 3. );
  TH1D *hPrimaryPtot = new TH1D("hPrimaryPtot","Primary track momentum;p (GeV/c);Entries",
                                200, 0., 3. );
  TH1D *hGlobalPt = new TH1D("hGlobalPt","Global track transverse momentum;p_{T} (GeV/c)",
                              200, 0., 2.5 );
  TH1D *hPrimaryPt = new TH1D("hPrimaryPt","Primary track transverse momentum;p_{T} (GeV/c)",
                              200, 0., 2.5 );
  TH1D *hPrimaryPx = new TH1D("hPrimaryPx","Primary track momentum p_{x};p_{x} (GeV/c);Entries",
                              400, -2., 2. );
  TH1D *hPrimaryPy = new TH1D("hPrimaryPy","Primary track momentum p_{y};p_{y} (GeV/c);Entries",
                              400, -2., 2. );
  TH1D *hPrimaryPz = new TH1D("hPrimaryPz","Primary track momentum p_{z};p_{z} (GeV/c);Entries",
                              400, -2., 2. );
  TH1D *hGlobalPx = new TH1D("hGlobalPx","Global track momentum p_{x};p_{x} (GeV/c);Entries",
                              400, -2., 2. );
  TH1D *hGlobalPy = new TH1D("hGlobalPy","Global track momentum p_{y};p_{y} (GeV/c);Entries",
                              400, -2., 2. );
  TH1D *hGlobalPz = new TH1D("hGlobalPz","Global track momentum p_{z};p_{z} (GeV/c);Entries",
                              400, -2., 2. );

  TH1D *hNHits = new TH1D("hNHits","Number of hits;nHits;Entries", 80, -0.5, 79.5);
  TH1D *hNHitsRatio = new TH1D("hNHitsRatio","nHitsFit to nHitsPoss ratio;nHitsFit/nHitsRatio;Entries",
                               10, 0., 1. );
  TH1D *hChi2 = new TH1D("hChi2","#chi^{2} of the track;#chi^{2};Entries",
                         200, 0., 20.);
  
  TH1D *hDcaZ = new TH1D("hDcaZ","DCA Z to primary vertex;DCA X (cm);Entries",
                        100, 0., 5.);
  TH1D *hDcaPerp = new TH1D("hDcaPerp","DCA Perp to primary vertex;DCA Perp (cm);Entries",
                        100, 0., 5.);
  TH1D *hDca = new TH1D("hDca","DCA Mag to primary vertex;DCA (cm);Entries",
                        100, 0., 5.);

  TH2D *hDcaVsPt[2];
  hDcaVsPt[0] = new TH2D("hDcaVsPt_0","p_{T} vs. DCA (positive particles);p_{T} (GeV/c);DCA (cm)",
                            350, 0, 3.5, 100, 0., 5.);
  hDcaVsPt[1] = new TH2D("hDcaVsPt_1","p_{T} vs. DCA (negative particles);p_{T} (GeV/c);DCA (cm)",
                            350, 0, 3.5, 100, 0., 5.);

  TH1D *hPhi = new TH1D("hPhi","Azimuthal angle distribution;#phi;Entries",
                        640, -3.2, 3.2 );
  TH1D *hEta = new TH1D("hEta","Track pseudorapidity;#eta;Entries", 220, -1.1, 1.1 );
  TH1D *hEtaG = new TH1D("hEtaG","Track pseudorapidity of global track;#eta;Entires", 220, -1., 1. );
  TH2D *hPtVsEta = new TH2D("hPtVsEta","p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)",
                            240, -1.2, 1.2, 350 , 0., 3.5);
  TH2D *hPrimaryPhiVsPt[2];
  for(int i=0; i<2; i++) {
    hPrimaryPhiVsPt[i] = new TH2D(Form("hPrimaryPhiVsPt_%d",i),
         Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
         350, 0., 3.5, 640, -3.2, 3.2 );
  }
  TH1D* hDedx = new TH1D("hDedx","dE/dx;dE/dx (keV/cm);Entries",
                         125, 0., 12.5);

  TH2D *hDedxVsPt[2], *hNSigmaPionVsPt[2], *hNSigmaElectronVsPt[2], *hNSigmaKaonVsPt[2], *hNSigmaProtonVsPt[2];
  const Char_t *PosNeg[] = {"positive","negative"};
  for( Int_t i = 0; i < 2; i++ ) {
    hDedxVsPt[i] = new TH2D(Form("hDedxVsPt_%i",i),Form("dE/dx vs. p_{T} %s;p_{T} (GeV/c);dE/dx (keV/cm)",PosNeg[i]),
                               420, 0., 3.1, 600, 0., 15.);
    hNSigmaPionVsPt[i] = new TH2D(Form("hNSigmaPionVsPt_%i",i),
                                  Form("n#sigma(#pi) vs. p_{T} %s;p_{T} (GeV/c);n#sigma(#pi)",PosNeg[i]),
                                  420, 0., 3.1, 300, -15., 15.);
    hNSigmaElectronVsPt[i] = new TH2D(Form("hNSigmaElectronVsPt_%i",i),
                                      Form("n#sigma(e) vs. p_{T} %s;p_{T} (GeV/c);n#sigma(e)",PosNeg[i]),
                                      420, 0., 3.1, 300, -15., 15.);
    hNSigmaKaonVsPt[i] = new TH2D(Form("hNSigmaKaonVsPt_%i",i),
                                  Form("n#sigma(K) vs. p_{T} %s;p_{T} (GeV/c);n#sigma(K)",PosNeg[i]),
                                  420, 0., 3.1, 300, -15., 15.);
    hNSigmaProtonVsPt[i] = new TH2D(Form("hNSigmaProtonVsPt_%i",i),
                                    Form("n#sigma(p) vs. p_{T} %s;p_{T} (GeV/c);n#sigma(p)",PosNeg[i]),
                                    420, 0., 3.1, 300, -15., 15.);
  }

  TH2D *hNSigmaElectronVsMassSqrt = new TH2D("hNSigmaElectronVsMassSqrt",
                                             "n#sigma(e) vs. Square mass;n#sigma(e);m^{2} (GeV/c^{2})^{2}",
                                              200, -10., 10., 250, 0., 12.5);
  TH2D *hNSigmaPionVsMassSqrt = new TH2D("hNSigmaPionVsMassSqrt",
                                         "n#sigma(#pi) vs. Square mass;n#sigma(#pi);m^{2} (GeV/c^{2})^{2}",
                                          200, -10., 10., 250, 0., 12.5);
  TH2D *hNSigmaKaonVsMassSqrt = new TH2D("hNSigmaKaonVsMassSqrt",
                                         "n#sigma(K) vs. Square mass;n#sigma(K);m^{2} (GeV/c^{2})^{2}",
                                          200, -10., 10., 250, 0., 12.5);
  TH2D *hNSigmaProtonVsMassSqrt = new TH2D("hNSigmaProtonVsMassSqrt",
                                           "n#sigma(p) vs. Square mass;n#sigma(p);m^{2} (GeV/c^{2})^{2}",
                                           200, -10., 10., 250, 0., 12.5);

  TH2D *hDedxVsPtPID[4];
  for ( int i=0; i<4; i++ ) {
    TString name = "hDedxVsPtPID_";
    name += i;
    TString title = "dE/dx vs. charge*p_{T} ";
    switch (i) {
      case 0: title += "|n#sigma(e)| #leq 2;"; break;
      case 1: title += "|n#sigma(#pi)| #leq 2;"; break;
      case 2: title += "|n#sigma(K)| #leq 2;"; break;
      case 3: title += "|n#sigma(p)| #leq 2;"; break;
      default: title += "unknown PID;";
    }
    title += "charge*p_{T} (GeV/c);dE/dx (keV/cm)";
    hDedxVsPtPID[i] = new TH2D(name.Data(), title.Data(),
                               840, -2.1, 2.1, 600, 0., 12.);
  }

  // TofPidTrait
  TH1D *hTofBeta = new TH1D("hTofBeta","BTofPidTraits #beta;#beta",
                            2000, 0., 2.);
  TH2D *hInvBetaVsPt = new TH2D("hInvBetaVsPt","1/#beta vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta",
                                840, -2.1, 2.1, 200, 0.8, 2.8);
  TH1D *hMassSqr = new TH1D("hMassSqr","m^{2};m^{2} (GeV/c^{2})^{2};dN/dm^{2} (entries)",
                            520, -0.1, 5.1 );

  TH2D *hMassSqrVsPt[2];
  hMassSqrVsPt[0] = new TH2D("hMassSqrVsPt_0","m^{2} vs. p_{T} (positive particles);p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",
                                420, 0., 2.1, 200, -0.2, 1.8);
  hMassSqrVsPt[1] = new TH2D("hMassSqrVsPt_1","m^{2} vs. p_{T} (negative particles);p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",
                                420, 0., 2.1, 200, -0.2, 1.8);

  TH2D *hDedxVsMassSqr[2];
  hDedxVsMassSqr[0] = new TH2D("hDedxVsMassSqr_0","dE/dx vs. mass^{2} charge>0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
             440, -0.4, 1.8, 250, 0., 12.5 );
  hDedxVsMassSqr[1] = new TH2D("hDedxVsMassSqr_1","dE/dx vs. mass^{2} charge<0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
             440, -0.4, 1.8, 250, 0., 12.5 );

  TH2D *hInvBetaDiffElectronVsPt = new TH2D("hInvBetaDiffElectronVsPt","1/#beta - 1/#beta(electron) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(e)",
                                            840, -2.1, 2.1, 200, -0.1, 0.1);
  TH2D *hInvBetaDiffPionVsPt = new TH2D("hInvBetaDiffPionVsPt","1/#beta - 1/#beta(pion) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(#pi)",
                                            840, -2.1, 2.1, 200, -0.1, 0.1);
  TH2D *hInvBetaDiffKaonVsPt = new TH2D("hInvBetaDiffKaonVsPt","1/#beta - 1/#beta(kaon) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(K)",
                                            840, -2.1, 2.1, 200, -0.1, 0.1);
  TH2D *hInvBetaDiffProtonVsPt = new TH2D("hInvBetaDiffProtonVsPt","1/#beta - 1/#beta(p) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(p)",
                                          840, -2.1, 2.1, 200, -0.1, 0.1);
  TProfile *hTrackProfile[6];
  hTrackProfile[0] = new TProfile("hTrackProfile_0","Profile of track #phi;Run ID;<#phi>",
                                           runIdBins, runIdRange[0], runIdRange[1] );
  hTrackProfile[1] = new TProfile("hTrackProfile_1","Profile of track p_{T};Run ID;<p_{T}>",
                                           runIdBins, runIdRange[0], runIdRange[1] );
  hTrackProfile[2] = new TProfile("hTrackProfile_2","Profile of track nHits;Run ID;<nHits>",
                                          runIdBins, runIdRange[0], runIdRange[1] );
  hTrackProfile[3] = new TProfile("hTrackProfile_3","Profile of track DCA;Run ID;<DCA>",
                                          runIdBins, runIdRange[0], runIdRange[1] );
  hTrackProfile[4] = new TProfile("hTrackProfile_4","Profile of track #beta;Run ID;<#beta>",
                                             runIdBins, runIdRange[0], runIdRange[1] );
  hTrackProfile[5] = new TProfile("hTrackProfile_5","Profile of track dE/dx;Run ID;<dE/dx> (keV/cm)",
                                            runIdBins, runIdRange[0], runIdRange[1] );

  for(int iTrk=0; iTrk<6; iTrk++) {
    Int_t mColor = 1;
    hTrackProfile[iTrk]->SetMarkerStyle(20);    // filled circle
    hTrackProfile[iTrk]->SetMarkerColor(mColor);
    hTrackProfile[iTrk]->SetMarkerSize(1.1);
    hTrackProfile[iTrk]->SetLineWidth(2);
    hTrackProfile[iTrk]->SetLineColor(mColor);  // black
  }

  TProfile *hSinPhi[3], *hCosPhi[3];
  for( Int_t i = 0; i < 3; i++) {
    hSinPhi[i] = new TProfile(Form("hSinPhi%i",i+1),
                              Form("Average sin(%i#phi) of Run ID;Run ID; <sin(%i#phi)>",i+1,i+1),
                              runIdBins, runIdRange[0], runIdRange[1]);
    hCosPhi[i] = new TProfile(Form("hCosPhi%i",i+1),
                              Form("Average cos(%i#phi) of Run ID;Run ID; <cos(%i#phi)>",i+1,i+1),
                              runIdBins, runIdRange[0], runIdRange[1]);

  }
  for( Int_t i = 0; i < 3; i++ ) {
    Int_t mColor = 1;
    hSinPhi[i]->SetMarkerStyle(20);    // filled circle
    hSinPhi[i]->SetMarkerColor(mColor);
    hSinPhi[i]->SetMarkerSize(1.1);
    hSinPhi[i]->SetLineWidth(2);
    hSinPhi[i]->SetLineColor(mColor);  // black
    hCosPhi[i]->SetMarkerStyle(20);    // filled circle
    hCosPhi[i]->SetMarkerColor(mColor);
    hCosPhi[i]->SetMarkerSize(1.1);
    hCosPhi[i]->SetLineWidth(2);
    hCosPhi[i]->SetLineColor(mColor);  // black
  }
  


	/*////////////////////////////////////////////////////////////////////////////////////////*/
 /*________________________________START OF EVENT LOOP_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
	for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

  	if ( iEvent % 10000 == 0) {
  		std::cout << "Working on event #[" << (iEvent+1)
     	      	  << "/" << events2read << "]" << std::endl;
    }

    if (iEvent == events2read-1) {
    	std::cout << "Working on event #[" << (events2read)
     	 	      	<< "/" << events2read << "]" << std::endl;
    }
	
		Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
    	std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      	break;
    }

   	// Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
   	StFemtoEvent *event = dst->event();
    if( !event ) {
    	std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    // Simple event cut
    if(CUTS == true && isGoodEvent( event, VtxZcut, VtxRcut, VtxYdelta  ) == false ) continue;

    TVector3 pVtx = event->primaryVertex();
            
    //RunQA
    if( useRunQA == true && std::find(badRuns.begin(), badRuns.end(), event -> runId()) != badRuns.end() ) continue;

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*________________________________FILL HISTOGRAMS_________________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

    RunID = event -> runId();
    Float_t bbcE = 0., bbcW = 0., bbcAdcSum = 0.;
    for( Int_t iTile = 0; iTile < 24; iTile++ ) {
      bbcE += event -> bbcAdcEast(iTile);
      bbcW += event -> bbcAdcWest(iTile);
      bbcAdcSum = bbcE + bbcW;
      hAdcBbcEastVsTile -> Fill( iTile, event -> bbcAdcEast(iTile) ); 
      hAdcBbcWestVsTile -> Fill( iTile, event -> bbcAdcWest(iTile) );
    }

    // Fill event histograms
    hNeventsVsRunId -> Fill( event -> runId() );
    hRefMult->Fill( event->refMult() );
    hRefMult2->Fill( event->refMult2() );
    hGRefMult->Fill( event->gRefMult() );

    hRefMultVsAdcZdcE->Fill( event->refMult(), event->zdcSumAdcEast() );
    hRefMultVsAdcZdcW->Fill( event->refMult(), event->zdcSumAdcWest() );
    hRefMultVsAdcBbcE->Fill( event->refMult(), bbcE );
    hRefMultVsAdcBbcW->Fill( event->refMult(), bbcW );

    hAdcZdcEast->Fill( event->zdcSumAdcEast() );
    hAdcZdcWest->Fill( event->zdcSumAdcWest() );
    hAdcZdcSum->Fill( event->zdcSumAdcEast() + event->zdcSumAdcWest() );
    hAdcBbcEast->Fill( bbcE );
    hAdcBbcWest->Fill( bbcW );
    hAdcBbcSum->Fill( bbcAdcSum );

    hAdcZdcEvsW -> Fill( event->zdcSumAdcEast(), event->zdcSumAdcWest() );
    hAdcBbcEvsW -> Fill( bbcE, bbcW );

    hVtxXvsRefMult->Fill( event->primaryVertex().X(), event->refMult() );
    hVtxYvsRefMult->Fill( event->primaryVertex().Y(), event->refMult() );
    hVtxZvsRefMult->Fill( event->primaryVertex().Z(), event->refMult() );

    hVtxX->Fill( event->primaryVertex().X() );
    hVtxY->Fill( event->primaryVertex().Y() );
    hVtxZ->Fill( event->primaryVertex().Z() );
    hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
    hVtxZvsX->Fill( event->primaryVertex().Z(), event->primaryVertex().X() );
    hVtxZvsY->Fill( event->primaryVertex().Z(), event->primaryVertex().Y() );

    hVpdVzVsVtxZ->Fill(event->vpdVz(), event->primaryVertex().Z());
    hVpdVzDiffVsVz->Fill( event->primaryVertex().Z(),
                          event->primaryVertex().Z() - event->vpdVz() );
    
    
    hNumberOfPrimaries->Fill( event->numberOfPrimaryTracks() );
    hNumberOfGlobals->Fill( event->numberOfGlobalTracks() );
    hCent9->Fill( event->cent9() );
    hCent16->Fill( event->cent16() );
    hBTofHit->Fill( event->numberOfBTofHit() );
    hBTofMatched->Fill( event->numberOfTofMatched() );
    //hBemcMatched->Fill( event->numberOfBEMCMatched() );
    hRanking->Fill( event->ranking() );
    
    hBTofTrayMultVsRefMult->Fill( event->refMult(),
                                  event->numberOfBTofHit() );
    hBTofMatchedVsRefMult->Fill( event->refMult(),
                                 event->numberOfTofMatched() );
    hTransSphericity->Fill( event->transverseSphericity() );
    hTransSphericity2->Fill( event->transverseSphericity2() );
    hNumberOfVertices->Fill( event->numberOfPrimaryVertices() );
    hEventProfile[0]->Fill( event->runId(), event->gRefMult() );
    hEventProfile[1]->Fill( event->runId(), event->numberOfBTofHit() );
    hEventProfile[2]->Fill( event->runId(), event->numberOfTofMatched() );
    hEventProfile[3]->Fill( event->runId(), event->numberOfPrimaryTracks() );
    hEventProfile[4]->Fill( event->runId(), event->numberOfGlobalTracks() );
    hEventProfile[5]->Fill( event->runId(), event->zdcSumAdcEast() + event->zdcSumAdcWest() );
    hEventProfile[6]->Fill( event->runId(), bbcAdcSum );
    hEventProfile[7]->Fill( event->runId(), pVtx.X() );
    hEventProfile[8]->Fill( event->runId(), pVtx.Y() );
    hEventProfile[9]->Fill( event->runId(), pVtx.Z() );

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*________________________________START OF TRACK LOOP_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;
      //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;

      if ( !femtoTrack->isPrimary() ) continue;

      // Simple single-track cut
      if( CUTS == true && isGoodTrack(event, femtoTrack, DCAcut) == false ) continue;

      hGlobalPtot->Fill( femtoTrack->gMom().Mag() );
      hPrimaryPtot->Fill( femtoTrack->pMom().Mag() );
      hGlobalPt->Fill( femtoTrack->gMom().Pt() );
      hPrimaryPt->Fill( femtoTrack->pMom().Pt() );
      hPrimaryPx->Fill( femtoTrack->pMom().X() );
      hPrimaryPy->Fill( femtoTrack->pMom().Y() );
      hPrimaryPz->Fill( femtoTrack->pMom().Z() );
      hGlobalPx->Fill( femtoTrack->gMom().X() );
      hGlobalPy->Fill( femtoTrack->gMom().Y() );
      hGlobalPz->Fill( femtoTrack->gMom().Z() );

      hNHits->Fill( femtoTrack->nHits() );
      hNHitsRatio->Fill( (Float_t)femtoTrack->nHitsFit()/femtoTrack->nHitsPoss() );
      hChi2->Fill( femtoTrack->chi2() );

      hDcaZ->Fill(femtoTrack->gDCAz( pVtx.Z() ) );
      hDcaPerp->Fill( femtoTrack->gDCAxy( pVtx.X(), pVtx.Y() ) );
      hDca->Fill( femtoTrack->gDCA( pVtx.X(), pVtx.Y(), pVtx.Z() ) );

      hPtVsEta->Fill( femtoTrack->pMom().Eta(), femtoTrack->pMom().Pt() );
      hPhi->Fill( femtoTrack->pMom().Phi() );
      hEta->Fill( femtoTrack->pMom().Eta() );
      hEtaG->Fill( femtoTrack->gMom().Eta() );
      hDedx->Fill( femtoTrack->dEdx() * 1e6 );

      // If electron has passed PID nsigma cut
      if ( TMath::Abs( femtoTrack->nSigmaElectron() ) <= 2. ) {
        hDedxVsPtPID[0]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->dEdx() * 1e6 );
      }

      // If pion has passed PID nsigma cut
      if ( TMath::Abs( femtoTrack->nSigmaPion() ) <= 2. ) {
        hDedxVsPtPID[1]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->dEdx() * 1e6 );
      }

      // If kaon has passed PID nsigma cut
      if ( TMath::Abs( femtoTrack->nSigmaKaon() ) <= 2. ) {
        hDedxVsPtPID[2]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->dEdx() * 1e6 );
      }

      // If proton has passed PID nsigma cut
      if ( TMath::Abs( femtoTrack->nSigmaProton() ) <= 2. ) {
        hDedxVsPtPID[3]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->dEdx() * 1e6 );
      }

      if( femtoTrack->charge() > 0 ) {
        hPrimaryPhiVsPt[0]->Fill( femtoTrack->pMom().Pt(),
                                  femtoTrack->pMom().Phi() );
        hDcaVsPt[0]->Fill( femtoTrack->pMom().Pt(), 
                           femtoTrack->gDCA( pVtx.X(), pVtx.Y(), pVtx.Z() ) );
        hDedxVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->dEdx() * 1e6 );
        hNSigmaElectronVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaElectron() );
        hNSigmaPionVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaPion() );
        hNSigmaKaonVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaKaon() );
        hNSigmaProtonVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaProton() );
      }
      else {
        hPrimaryPhiVsPt[1]->Fill( femtoTrack->pMom().Pt(),
  		                            femtoTrack->pMom().Phi() );
        hDcaVsPt[1]->Fill( femtoTrack->pMom().Pt(), 
                           femtoTrack->gDCA( pVtx.X(), pVtx.Y(), pVtx.Z() ) );
        hDedxVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->dEdx() * 1e6 );
        hNSigmaElectronVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaElectron() );
        hNSigmaPionVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaPion() );
        hNSigmaKaonVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaKaon() );
        hNSigmaProtonVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->nSigmaProton() );
      }

      hTrackProfile[0]->Fill( event->runId(), femtoTrack->pMom().Phi() );
      hTrackProfile[1]->Fill( event->runId(), femtoTrack->pMom().Pt() );
      hTrackProfile[2]->Fill( event->runId(), femtoTrack->nHits() );
      hTrackProfile[3]->Fill( event->runId(),
                              femtoTrack->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z() ) );
      hTrackProfile[5]->Fill( event->runId(),
                              femtoTrack->dEdx() * 1e6 );

      for( Int_t i = 0; i < 3; i++) {
        hSinPhi[i]->Fill( event->runId(), sin( (Float_t)(i+1) * femtoTrack->phi() ) );
        hCosPhi[i]->Fill( event->runId(), cos( (Float_t)(i+1) * femtoTrack->phi() ) );
      }

      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {
  	    hTofBeta->Fill( femtoTrack->beta() );
        hInvBetaVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                            femtoTrack->invBeta() );

        hMassSqr->Fill( femtoTrack->massSqr() );
        if( femtoTrack -> charge() > 0) {
          hMassSqrVsPt[0]->Fill( femtoTrack->pMom().Pt(), femtoTrack->massSqr() );
        }
        else {
          hMassSqrVsPt[1]->Fill( femtoTrack->pMom().Pt(), femtoTrack->massSqr() );
        }

        hNSigmaElectronVsMassSqrt->Fill( femtoTrack->nSigmaElectron(), femtoTrack->massSqr() );
        hNSigmaPionVsMassSqrt->Fill( femtoTrack->nSigmaPion(), femtoTrack->massSqr() );
        hNSigmaKaonVsMassSqrt->Fill( femtoTrack->nSigmaKaon(), femtoTrack->massSqr() );
        hNSigmaProtonVsMassSqrt->Fill( femtoTrack->nSigmaProton(), femtoTrack->massSqr() );

        hInvBetaDiffElectronVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                        femtoTrack->invBeta() - TMath::Sqrt(electron_mass_sqr +
                                        femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
        hInvBetaDiffPionVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                    femtoTrack->invBeta() - TMath::Sqrt(pion_mass_sqr +
                                    femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
        hInvBetaDiffKaonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                    femtoTrack->invBeta() - TMath::Sqrt(kaon_mass_sqr +
                                    femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
        hInvBetaDiffProtonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                      femtoTrack->invBeta() - TMath::Sqrt(proton_mass_sqr +
                                      femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
        if ( femtoTrack->charge() > 0 ) {
          hDedxVsMassSqr[0]->Fill( femtoTrack->massSqr(), femtoTrack->dEdx() * 1e6 );
        }
        else {
          hDedxVsMassSqr[1]->Fill( femtoTrack->massSqr(), femtoTrack->dEdx() * 1e6 );
        }
        hTrackProfile[4]->Fill( event->runId(), femtoTrack->beta() );
      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
  }// for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->Write();
  outFile->Close();

  femtoReader->Finish();
  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << std::endl;

}// void FemtoDstAnalyzer()





  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

//********************CHECK EVENT ON GOOD********************//
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZcut, Double_t VtxRcut, Double_t VtxYdelta ) {
  Bool_t check = true; 
  TVector3 pVtx = event->primaryVertex();

  // Reject vertices that are far from the central membrane along the beam
  if( TMath::Abs( pVtx.Z() ) > VtxZcut ) check = false;
  //if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + 0.8847, 2) ) > 1 ) check = false; 14.5
  if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + VtxYdelta , 2) ) > VtxRcut ) check = false;
  
  return check;
}// isGoodEvent(){}

//********************CHECK TRACK ON GOOD********************//
Bool_t isGoodTrack(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCAcut) {
  Bool_t check = true;
  TVector3 pVtx = event->primaryVertex();

  if ( !femtoTrack ) check = false;
  // Must be a primary track
  if ( !femtoTrack->isPrimary() ) check = false;
  if ( ( femtoTrack -> dEdx() ) == 0 ) check = false;
  // Simple single-track cut
  if( femtoTrack -> gMom().Mag() < 0.1 || femtoTrack -> gDCA(pVtx).Mag() > DCAcut ) check = false;    
  if( TMath::Abs( femtoTrack -> eta() ) > 1.0 || femtoTrack -> nHits() < 15 || 
                  femtoTrack -> pt() < 0.2) check = false;
  if(femtoTrack -> pt() > 2.0) check = false; 
  if(  ( (Double_t)femtoTrack -> nHits() )/( (Double_t)femtoTrack -> nHitsPoss() )  < 0.52 ) check = false;

  return check; 
}// isGoodTrack(){}
