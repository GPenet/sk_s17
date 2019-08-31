struct XINDEX3 {
	uint32_t cellsbf; //cells bit field in band  mode
	uint32_t ideb;// index in main table
	uint32_t tcells[3];
	void Open(uint32_t filter, int ind,int * tc) {
		cellsbf= filter;
		ideb  = ind;
		memcpy(tcells, tc, sizeof tcells);
	}
};
struct X_EXPAND_3_5 {
	uint32_t cellsbf,d; //0;1;2 cells
	inline void Add0() { cellsbf = 0;	d = 0; }
	inline void Add1(int i1) { cellsbf = 1 << i1; d = (i1<<8 )| 1;
	}
	inline void Add2(int i1,int i2){ cellsbf = (1 << i1) | (1 << i2);
	d = (i1 << 16) | (i2 << 8) | 2;
	}
};

struct SPOT17_NOUAS2_3{//hosting a classical search with stack limit control
	int known, active;
	GINT64  stacks;
	int * tua, nua;
	void newspot(int * oua, int onua, SPOT17_NOUAS2_3 * sp, int cell, int bit = 0);
};

struct G17TMORE{// FIFO table of more for bands 1+2
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init(){ maxt = G17MORESIZE; nt = 0; }
	inline void Add(uint64_t v){//add a new more in FIFO 
		if (nt < maxt){// use next location
			curt = nt;
			t[nt++] = v;
		}
		else{// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}
	inline void Add_If_New(uint64_t v){// check oldest first
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// from curt to 0
		if ((*Rt) == V)return;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) goto add;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)if ((*Rt) == V)return;
	add:
		Add(v);
	}
	inline int Check(uint64_t v){// check oldest first
		if (!nt) return 0;
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}
	void Print(int modegua){
		register uint64_t *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return;

		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}

	}


	void PrintUasDirect(){
		for (int i = 0; i < nt; i++){
			register uint64_t w = t[i];
			cout << Char2Xout(w) << endl;
		}
	}

}g17more, g17morebig;

struct MORE32 {// FIFO table of more for band b
	uint32_t  t[32];
	int nt, maxt, curt;
	inline void Init() { maxt = 32; nt = 0; }
	inline void Add(uint32_t v) {//add a new more in FIFO 
		if (nt < maxt) {// use next location
			curt = nt;
			t[nt++] = v;
		}
		else {// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}

	inline int Check(uint32_t v) {// check oldest first
		if (!nt) return 0;
		register uint32_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt > Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}

};

struct G17TB3GO{
	int ib3, minirows_triplets;
	GINT64 pairs, // 27 bits  pairs holes 27 bits pairs cells(later pairs+triplets)
		triplets,
		count, // 4 items 9 bits/9 minirows 1;2;3 pairs any number of pairs
		countsum,// 4 values (4x16 bits) min clues total and per stack
		countstack;// count of clues per stack {bands 12 + countsum)
	void Debug();
}g17tb3go[512];

//========== indexstep max number of blocs 
#define G17BLOCSUA 4
#define G17BLOCGSUA 10
//640 uas 320 guas ,,, step size

#define G17CHKX 256
#define G17CHKY 256
#define MAXNIND6 20000
#define NVUAS128 4
#define NVGUAS128 15
struct V256_UAS { 
	BF128 v[NVUAS128];
	inline void operator&= (const V256_UAS &r) {
		v[0] &= r.v[0]; v[1] &= r.v[1];
		// assuming NVUAS128 4
		v[2] &= r.v[2]; v[3] &= r.v[3];

	}
	void Debug(const char * lib,int mirror=0);
	void Fout(const char * lib);
	void Cout();
};
struct V256_GUAS { 
	BF128 v[NVGUAS128]; 
	inline void operator&= (const V256_GUAS &r) {
		for (int i = 0; i < NVGUAS128; i++)v[i] &= r.v[i];
	}
};

//#include "g17_debuggingcode.cpp"
// standard first band (or unique band)
struct STD_B416 {
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i);
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas()
		;
	void InitC10(int i);
	void InitG12(int i);
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	XINDEX3 xindex3[2000];
	X_EXPAND_3_5 x_expand_3_5[70000]; //51520 for 5 alone
	int nxindex3,n3_5;
	// row solution pattern in digit
	int mini_digs[9], mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27], nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb);
	void DoExpandBand();	
	void DebugExpand();
	void PrintShortStatus();
};

struct STD_B3 :STD_B416 {// data specific to bands 3
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket4;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81];
		int ua2_i27[81];
	}guas;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	//BF128 tbands_UA4_6s, tbands_pairs, tbands_triplets;
	//int tuas46[81];
	//_______________________
	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	int IsGua(int i81);
	int IsGua3(int i81);
	int GetI81_2(int bf) {
		for (int i = 0; i < 81; i++) if (guas.ua_pair[i] == bf) return i;
		return -1;
	}
	int GetI81_3(int bf) {
		for (int i = 0; i < 81; i++) if (guas.ua_triplet[i] == bf) return i;
		return -1;
	}
	void PrintB3Status();
};

//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
		revised_gangbf[9],// same revised UA2s UA3s ***
		mini_digs[9], mini_pairs[27], // UA2s UA3  ***
		//valid_pairs, //  27 bits valid sockets UA2s ***
		nfloors, limstep, map[9], cptdebug, modemore;
	BF128 valid_sockets;

	//=============== uas collector 
	int limsize, floors;
	uint64_t  tuaold[1000],// previous non hit uas infinal table of uas for bands 1+2
		tua[TUA64_12SIZE],// 
		tuab1b2[200];// collecting bands uas in 2x mode
	uint32_t nuaold, nua, nuab1b2,
		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	int ib, digp;
	uint64_t w0, ua;
	//_____________________ functions collect UAs bands 1+2
	int Initgen();
	void BuildFloorsAndCollectOlds(int fl);
	//int AddUA64(uint64_t * t, uint32_t & nt);
	inline void AddUA(uint64_t v) {
		ua = v; AddUA64(tua, nua, ua);
	}
	inline void AddUACheck(uint64_t v) {
		if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
		ua = v; AddUA64(tua, nua, ua);
	}
	void BuilOldUAs(uint32_t r0);
	int CheckOld();
	int CheckMain(uint64_t wua);
	void CollectMore();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	//_____________________ functions collect UA2s UAs3 socket 

	void ProcessSocket2(int i81);
	int DebugUas();
};

#define SIZETGUA 150
#define GUAREDSIZE 100
struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	int   *gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// actice cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t tua[SIZETGUA];
		int col1, col2;// columns of the socket
		int i_81; // index 0_80 for this 
		int i9;// complementary column in minirow
		int id1, id2; // index of digits in gang 27 
		// Current band1+2 data
		int digs, dig1, dig2;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		int gangcols[9];// revised gangster
		uint32_t nua;// nua_start, nua_end;
		void Debug(const char * lib);

	}tsgua2[81];
	struct SGUA3 {// 81 possible UA3 sockets
		// permanent data
		uint64_t tua[SIZETGUA];
		int col1;// first columns 0-9 
		int i_81, imini; // index 0_80 for this 
		int id1, id2, id3; // index of digits in gang 27 
		// Current band1+2 data
		int  dig1, dig2, dig3;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		uint32_t nua;// nua_start, nua_end, nua;
		void Debug(const char * lib);
	}tsgua3[81];
	// __________________________  primary UAs tables and creation of such tables
	uint64_t  // tua3x[3000],// dynamic sub tables
		*ptua2;// pointer to current table cycle search 2/3
	uint32_t  ntua2, ntua3, nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int tactive2[81], nactive2, tactive3[81], nactive3;
	int   tcolok[2], ncolok;

	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		gang27 = gang[0];
		InitialSockets2Setup();
		InitialSockets3Setup();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void Build_CheckUAs_Subsets_Table();
	void Build_tuacheck(int fl);
	int Have_tuacheck_subset();
	void SecondSockets2Setup();// band1+2 level
	void SecondSockets2MoreUAs();// band1+2 level
	void GuaCollectMore();
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl, int diag = 0);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	int DebugFreshUA(uint64_t ua);
	int Debug17(SGUA2 & w);
	//int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3


	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);
};

struct G17B3HANDLER {
	int known_b3, rknown_b3, active_b3,   ib3, nb3,
		active_sub, ndead, wactive0, nmiss, ncritical,
		irloop, wua, stack;
	uint32_t mini_bf1, mini_bf2, mini_bf3,pairsbf,  pairs27, mini_triplet,
		*uasb3if, nuasb3if, *uasb3of, nuasb3of,andoutf;
	GINT64 stack_count;
	int diagh;
	// ================== entry in the proces
	void Init(int i);
	uint32_t IsMultiple(int bf);
	int BuildIfShortB3();
	int ShrinkUas1();
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs();
	void Go_Critical(uint32_t * wua=0);
	void CriticalLoop();
	void CriticalExitLoop();
	void Critical_0_UA();
	void CriticalFinalCheck();
	//===================== process not critical
	void Go_Not_Critical_missn();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini( int M, int mask);
	void Go_Subcritical();
	void Go_SubcriticalMiniRow();
	//===============  debugging 
	void PrintStatus();

};

struct BANDS_AB {// handling bands 12 in A B mode
	struct BANDB { // when band a is filled
		uint32_t tua[512],nua;
		uint32_t pairs27, critbf, mini1, mini2, mini3, triplet,
			mini_all,ncrit,nmiss,ncb2;
		// uas in field outfiel
		uint32_t *tuaif, *tuaof, nuaif, nuaof, andoutf;
		// expansion 
		int known_bb, rknown_bb, active_bb,
			active_sub, ndead,  wactive0,
			irloop,diag,diagbug;
		uint32_t  wua;


		void AddIF(uint32_t ua) {
			if (nuaif >= 256)nuaif = 255;
			ua &= BIT_SET_27;
			ua |= _popcnt32(ua) << 27;
			AddUA32(tuaif, nuaif,ua);
		}
		void AddOF(uint32_t ua) {
			if (nuaof >= 256)nuaof = 255;
			ua &= BIT_SET_27;
			ua |= _popcnt32(ua) << 27;
			AddUA32(tuaof, nuaof,ua);
		}
		inline void Ncrit() {
			mini_all = mini1 | mini2 | mini3 | triplet;
			ncrit= _popcnt32(mini_all)+ _popcnt32(mini3);
			nmiss = ncb2 - ncrit;
		}
		//void Init(int i);
		int BuildIF_short();
		int ShrinkUas1();
		void Go();
		//=============== process critical
		void CriticalAssignCell(int Ru);
		void Go_Critical(uint32_t * wua = 0);
		void CriticalLoop();
		void CriticalExitLoop();
		void Critical_0_UA();
		//===================== process not critical
		void Go_Not_Critical_missn();
		//==================== process subcritical no cell added outside the GUAs field
		void SubMini(int M, int mask);
		void Go_Subcritical();
		void Go_SubcriticalMiniRow();
		int IsFinalOrMultiple(uint32_t * wua=0);
		void Status();
		void DebugIfOf();
	}sbb;

	uint32_t ni3, mode_ab, ia, ib,myuab,
		indd,indf,filt32,ncluesbandb,stack_filter;
	GINT64 stack_countba, stack_count, stack_countf;
	XINDEX3 * myi3,wi3;
	X_EXPAND_3_5 * myt3_5,wi3_5;
	STD_B1_2 * mybb;
	//======= band B initial infield outfield and more outfield table
	uint32_t btuaif[256], btuaof[3000], tuaif[3000],
		nbif,nbof;
	uint32_t more_of[128], nmoreof, more_if[128], nmoreif;
	//============== reduction of UAs GUAs
	uint32_t tuasmini[36][100], ntuasmini[36], ntua,
		activemini[36],nactivemini;
	uint64_t tua[512];
	//__________________________ secondary guas table 
	uint32_t tuar2[81][GUAREDSIZE], tuar3[81][GUAREDSIZE],
		ntuar2[81], ntuar3[81];
	uint32_t guar2i81[81], guar3i81[81], nguared_2, nguared_3;
	BF128 forced81_2, forced81_3, final81_2, final81_3;
	//========== tclues for valid XY 
	uint32_t tclues[40];// mini 25+band a
	int ncluesa, nclues;
	//==================== current band 3 to process
	uint32_t tcluesb12[20], ncluesb3;
	int  ntb3, nmiss;
	uint32_t mini_bf1, mini_bf2, mini_bf3, pairsbf, pairs27, mini_triplet;
	uint32_t all_used_minis, mincount;
	uint32_t uasb3_1[2000], uasb3_2[2000],uas_in[2000], nuasb3_1, nuasb3_2,nuas_in;

	void Go(STD_B1_2 & ba, STD_B1_2 & bb, int i, int mode);
	void AddMini(int imini, uint32_t ua);
	void Init3clues();
	int Init3_5clues();
	int IsMultiple(int bf,int diag=0);
	void CriticalFinalCheck(int bf);
	void GuasCollect(int bf);
	int SetUpGuas2_3(int ib3);
	int EndCollectBand3(int ib3);
	void GoBand3(int ib3);
	void ExpandBand3();
	int BuildUasB3_in(uint32_t known, uint32_t field);
	void DebugInit();
	void Status();
};

  
struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	BF128 p17diag;// known 17 pattern for tests
	int b3lim,debug17, diag,diagbug,aigstop,
		npuz, a_17_found_here;
	BANDS_AB bands_ab;
	G17B3HANDLER g17hh0;
	MORE32 moreb;
	//______sockets common to  all bands 3  
	BF128 isguasocket2all, isguasocket3all;
	//=====================process
	void GoM10();// end preparation for a given band1+band2+ table band3 pass 656 566
	void Go();
	inline void SetUp(BANDS_AB::BANDB * bb){
		bb->tuaof = bands_ab.btuaof;
		bands_ab.ncluesb3 = 6;
	}
		
	//================ debugging code
	int DebugK17M10();
	void GodebugInit(int mode);
	int GodebugFindKnown17();
	int GodebugCheckUas(const char * lib);
};