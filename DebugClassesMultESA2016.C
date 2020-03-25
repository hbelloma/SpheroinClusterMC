class ESATrack : public TObject
{
	public:
		Float_t   p;
		Float_t   pTPC;
		Float_t   pt;
		//  Float_t   ptcon;
		Float_t   pttrue;
		//  Float_t   tpcchi;
		Float_t   eta;
		Float_t   phi;
		Float_t   dedx;
		Float_t beta;
		Float_t   flength;
		Float_t   dcaxy;
		Float_t   dcaz;
		Int_t     mother; // pdg of mother (can be same particle)
		Short_t   q;
		Short_t   filter;
		Short_t   ncl;
		Float_t   nCr;
		UShort_t  nCl;

		//Short_t   neff;
		Short_t   nclf;
		Short_t   pid;
		Short_t   primary;  
		Short_t   order;

		Int_t     tpcnclS; //number of shared TPC clusters
		//July 2014
		Float_t ncrossrowstpc;//number of crossed rows in the TPC
		Float_t  ratiocrossrowstpcover; // ratio crossed rows over findable clusters in the TPC
		Float_t fraclusterstpcshared; //fracction of TPC shared clusters
		Float_t chi2perclusterits;  // chi2 per ITS cluster
		Float_t chi2perclustertpc;  // chi2 per TPC cluster
		Bool_t iskinkdaughter; //is kink daughter?
		Bool_t istpcrefit; //pass the TPC refit?
		Bool_t isitsrefit; //pass the ITS refit?
		Bool_t hasspdpoint; // has a SPD point?
		//Bool_t dcaxycut;    // Pt-dependent  dca cut
		Float_t TPCsigma_pi; 
		Float_t TPCsigma_ka; 
		Float_t TPCsigma_pr; 
		Float_t chi2tpcconstrainedvsglobal; //the golden cut require chi2tpcconstrainedvsglobal to be smaller than 36



		ESATrack();

		ClassDef(ESATrack, 1);    // Help class
};
//_________________________________________________________
class VZEROCell : public TObject
{
	public:

		Int_t   cellmult;
		Float_t cellindex;
		VZEROCell();

		ClassDef(VZEROCell, 1);    // Help class
};


//_____________________________________________________________________________
class ESAV0 : public TObject
{
	public:
		Float_t   p;
		Float_t   pt;
		Float_t   eta;
		Float_t   phi;
		Float_t   pdca;     // Distance of Closest Approach for positive track
		Float_t   ndca;     // Distance of Closest Approach for positive track
		Float_t   dmassG;
		Float_t   dmassK0;
		Float_t   dmassL;
		Float_t   dmassAL;
		Float_t   alpha;
		Float_t   ptarm;
		Float_t   decayr;
		Float_t   decayl;
		Float_t   chi2;
		Float_t   cospt;
		Float_t   dcav0;
		Float_t   dcadaughters;
		Int_t     pdg;
		Short_t   primary;  
		Short_t   status;  
		ESATrack ptrack;
		ESATrack ntrack;

		ESAV0();

		ClassDef(ESAV0, 1);    // Help class
};


//_____________________________________________________________________________
class ESATrackMC : public TObject
{
	public:
		Float_t pMC;
		Float_t ptMC;
		Float_t etaMC;
		Float_t phiMC;
		Float_t yMC;
		Short_t qMC;
		Short_t pidMC;
		Short_t orderMC;
		Int_t   pdgMC;

		ESATrackMC();

		ClassDef(ESATrackMC, 1);    // Help class for MC track debug info
};

//_____________________________________________________________________________
class ESAEvent : public TObject
{
	public:
		ULong64_t eventid;     // unique event id
		Int_t     run;         // run number
		UInt_t    time;        // time of event
		Float_t   v0Mperc;     // centrality V0A+V0C, default
		Float_t   v0Aperc;     // centrality V0A
		Float_t   adMperc;     // centrality ADA+ADC
	        Float_t   mag;         // magnetic field
		Float_t   zvtxTrk;        // rec vertex tracks
		Float_t   zvtxSpd;        // rec vertex SPD
		Float_t   zvtxMC;      // MC true vertes
		Float_t   ptmax;       // Max pt of tracks for this event
		Float_t   ptmaxMC;     // Max pt of MC tracks

		Short_t   vtxstatus;   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)

		Short_t   vtxstatus2015pp;   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)

		Int_t   trackmult05;   // Track mult (eta smaller than 05)
		Int_t   trackmult08;   // Track mult (eta smaller than 08)
		Short_t   n;           // Number of added tracks 
		Short_t   trackmultMC; // MC track mult (primary tracks)
		Short_t   nMC;         // MC number of added tracks 
		Short_t   process;     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
		Short_t   trig;        // 0=untriggered, &1 = MB, &2=V0 AND

		Bool_t   pileup;      // Is the event marked as pileup?


		Bool_t	bgreject;


		Bool_t	isincompletedaq;


		ESAEvent();

		ClassDef(ESAEvent, 1);    // Help class
};

//_____________________________________________________________________________
ClassImp(ESATrack)

	ESATrack::ESATrack():
		TObject(),
		p(-1),
		pTPC(-1),
		pt(-1),
		//  ptcon(-1),
		pttrue(-1),
		//  tpcchi(-1),
		eta(-999),
		phi(-999),
		dedx(-999),
		beta(-99),
		flength(-9999),

		dcaxy(-999),
		dcaz(-999),
		mother(0),
		q(-999),
		filter(-999),
		ncl(-999),
		nCr(-999),
		nCl(-999),
		//neff(-999),
		nclf(-999),
		pid(-999),
		primary(-999),
		order(-1),
		//filterset1(0),
		// filterset2(0),
		//filterset3(0),
		tpcnclS(0),
		//July 2014
		ncrossrowstpc(-999),//number of crossed rows in the TPC
		ratiocrossrowstpcover(-999), // ratio crossed rows over findable clusters in the TPC
		fraclusterstpcshared(-999), //fracction of TPC shared clusters
		chi2perclusterits(-999),  // chi2 per ITS cluster
		chi2perclustertpc(-999),  // chi2 per TPC cluster
		iskinkdaughter(0), //is kink daughter?
		istpcrefit(0), //pass the TPC refit?
		isitsrefit(0), //pass the ITS refit?
		hasspdpoint(0), // has a SPD point?
		TPCsigma_pi(-999),
		TPCsigma_ka(-999),
		TPCsigma_pr(-999),
		chi2tpcconstrainedvsglobal(-999) //the golden cut require chi2tpcconstrainedvsglobal to be smaller than 36




{
	// default constructor
}
//_____________________________________________________________________________
ClassImp(VZEROCell)

	VZEROCell::VZEROCell():
		TObject(),
		cellmult(-1),
		cellindex(-999)

{
	// default constructor
}

//_____________________________________________________________________________
ClassImp(ESAV0)

	ESAV0::ESAV0():
		TObject(),
		p(-1),
		pt(-1),
		eta(-999),
		phi(-999),
		pdca(-1),
		ndca(-1),
		dmassG(-1),
		dmassK0(-1),
		dmassL(-1),
		dmassAL(-1),
		alpha(-999),
		ptarm(-999),
		decayr(-999),
		decayl(-999),
		chi2(-1),
		cospt(-999),
		dcav0(999),
		dcadaughters(999),
		pdg(0),
		primary(-1),  
		status(0),  
		ptrack(),
		ntrack()
{
	// default constructor
}

//_____________________________________________________________________________
ClassImp(ESATrackMC)

	ESATrackMC::ESATrackMC():
		TObject(),
		pMC(-1),
		ptMC(-1),
		etaMC(-99),
		phiMC(-99),
		yMC(-99),
		qMC(-99),
		pidMC(-99),
		orderMC(-1),
		pdgMC(0)
{
	// default constructor
}

//_____________________________________________________________________________
ClassImp(ESAEvent)

	ESAEvent::ESAEvent():
		TObject(),
		eventid(0),      // unique event id
		run(-1),         // run number
		time(-1),        // time of event
		v0Mperc(1000),   // VZERO
		v0Aperc(1000),   // VZEROA
		adMperc(1000),   // AD

		mag(+999),       // magnetic field
		zvtxTrk(+999),   // rec vertex tracks
		zvtxSpd(+999),   // rec vertex SPD

		zvtxMC(+999),    // MC true vertes
		ptmax(-1),       // Max pt of tracks for this event
		ptmaxMC(-1),     // Max pt of MC tracks

		vtxstatus(-2),   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)
		vtxstatus2015pp(-999),   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)

		trackmult05(-1),   // Track mult (no cuts)
		trackmult08(-1),   // Track mult (no cuts)
		n(-1),           // Number of added tracks 
		trackmultMC(-1), // MC track mult (primary tracks)
		nMC(-1),         // MC number of added tracks 
		process(-2),     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
		trig(-1),        // Was the event triggered

		pileup(kFALSE),       // Is the event marked as pileup?

		bgreject(-1),


		isincompletedaq(kFALSE)//




{
	// default constructor
}
