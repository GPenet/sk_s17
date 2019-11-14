struct XINDEX3 {
	uint32_t cellsbf; //cells bit field in band  mode
	uint32_t ideb5,ideb;// index in main table
	uint32_t tcells[3];
	void Open(uint32_t filter, int ind,int * tc) {
		cellsbf= filter;
		ideb  = ind;
		memcpy(tcells, tc, sizeof tcells);
	}
	void Open5_6(uint32_t filter, int ind5, int ind6) {
		cellsbf = filter;
		ideb = ind6;
		ideb5 = ind5;
	}
	void status( int i) {
		cout << Char27out(cellsbf)
			<< " d5=" << ideb5 << "\td6=" << ideb
			<< " i=" << i<< endl;
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
struct VECT256 {// vector for other uas after 3 clues
	BF128 v[2];
	inline int IsEmpty(VECT256 &v2) {
		BF128 w = v[0] & v2.v[0];
		w |= v[1] & v2.v[1];
		return w.isEmpty();
	}
	inline void Init(int n) {
		if (n > 128) { v[0].SetAll_1(); v[1] = maskLSB[n - 128]; }
		else { v[1].SetAll_0(); v[0] = maskLSB[n]; }
	}
	inline void Set(int i, int j) { v[i].Set(j); }
	inline void Setx(int x) { Set((x >> 7), (x & 127)); }
	inline void Clear(int i, int j) {	v[i].clearBit(j);	}
	inline void Clearx(int x) { Clear((x >> 7), (x & 127)); }
	inline void And(VECT256 &v2) {	v[0] &= v2.v[0]; v[1] &= v2.v[1];}
	inline int On(int i,int j){return  v[i].On(j);	}
	inline void Table(int * t, int & n) {
		n = v[0].Table128(t);
		if (v[1].isNotEmpty()) {
			int tu2[128], ntu2 = v[1].Table128(tu2);
			for (int i = 0; i < ntu2; i++)
				t[n++] = tu2[i] + 128;
		}
	}
	inline void Not(VECT256 &v2) {	v[0] -= v2.v[0]; v[1] -= v2.v[1];	}
	void Print(const char * lib) {
		cout << "V256 status for " << lib<<endl;
		uint64_t *w = v[0].bf.u64, id = 0;
		for (int i = 0; i < 4; i++, id += 64, w++)
			if (*w)cout << Char64out(*w) << " " << id << endl;
	}
};
struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount,minplus;
	void SetMincount() {// after direct setting minis
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs

		// set up pair + triplet bitfield
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					critbf |= field;
		}
		minplus = mincount;
	}
	GINT64_t  Count_per_stack() {
		GINT64 cc; cc.u64 = 0;
		for (int i = 0, st = 0111; i < 3; i++, st <<= 1) {
			cc.u16[i]= _popcnt32(all_used_minis&st) + 
				_popcnt32(mini_bf3&st);
		}
		return cc;
	}
	void Status(const char * lib) {
		cout<<lib << "critical Status mincount ="<<mincount<< " minplus=" <<minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
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
	void Print(int modegua) {
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
		while (--Rt > Rtl) {
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}

	}


	void PrintUasDirect() {
		for (int i = 0; i < nt; i++) {
			register uint64_t w = t[i];
			cout << Char2Xout(w) << endl;
		}
	}

};
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
struct GUAN {// one occurrence of the active guan 
	uint64_t *tua, killer;
	uint32_t colbf, digsbf, ncol, // 2-4 cols
		 nua,  nfree;
	int i81;// i81 or pattern (if 4 columns)
	void Enter(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		i81 = ind;// index 0-80 if gua2 gua3
		tua = t; nua = n; colbf = cbf; digsbf = dbf;
		killer = BIT_SET_2X;
		for (uint32_t i = 0; i < nua; i++)killer &= tua[i];
		ncol = _popcnt32(colbf);
	}
	void Debug1Guan(int i,int all=0) {
		cout << Char2Xout(killer);
		cout << " i=" << i << "\tcols" << Char9out(colbf);
		cout << "\tdigs" << Char9out(digsbf) << " ncol=" << ncol
			<< " nua=" << nua << " i81=" << i81
			<< endl;
		if (all) {
			cout << Char2Xout(killer) << " kill " << endl;
			for (uint32_t i = 0; i < nua; i++)
				cout << Char2Xout(tua[i]) << endl;
		}
	}

};

// standard first band (or unique band)
struct STD_B416 {
	XINDEX3 *myi3;// index
	uint32_t *mybv5, *mybv6;// valid bands
	int nmyi3, nmybv5, nmybv6;

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
	void ExpandBand();
	void DebugExp(){
		cout << band<<" band id " << i416;
		cout << "\tn5=" << nmybv5;
		cout << "\tn6=" << nmybv6 << endl;
	}
	void DebugIndex() {
		cout << "debug index bande" << endl;
		for (int i = 0; i < nmyi3; i++)
			myi3[i].status(i);
	}
	void Debug5() {
		int indf5 = myi3[nmyi3 - 1].ideb5;
		cout << "5 clues table indf=" << indf5 << endl;
		for (int i = 0; i < indf5; i++) {
			cout << Char27out(mybv5[i]) << endl;
			if(_popcnt32(mybv5[i])!=5)
				cout << Char27out(mybv5[i])
				<< " bug in count"<< endl;
		}

	}
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
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
	void PrintShortStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket2_46;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81];
	}guas;
	struct G3X3X {// active guas at 3x3y
		BF128 g46, g23, gall, g3x4; 
	}g3x3yx;
	struct G3X {// active guas at 3x3y
		BF128 g46, g23, gall;
	}g3x;
	BF128  g3x23;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	MINCOUNT smin;
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
	void SetUpMincountxy(BF128 & validsockets);
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
	uint32_t patb, ib, digp,colb, cola;
	uint64_t w0, ua,p12;
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
	int BuilOldUAs(uint32_t r0);
	int CheckOld();
	int CheckMain(uint64_t wua);
	void CollectMore2digits();
	void Collect2digits2_4_cols();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	//_____________________ functions collect UA2s UAs3 socket 

	void ProcessSocket2(int i81);
	int DebugUas();
};

#define SIZETGUA 150
struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore,diagbug,
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
	int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t * tua;
		int col1, col2;// columns of the socket
		int i_81,iguan; // index 0_80 for this 
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
		uint64_t * tua,killer;
		int col1;// first columns 0-9 
		int i_81, imini,iguan; // index 0_80 for this 
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
	MINCOUNT smin;
	int known_b3, rknown_b3, active_b3,   ib3, nb3,
		active_sub, ndead, wactive0, nmiss, //ncritical,
		irloop, wua, stack;
	uint32_t *uasb3if, nuasb3if, *uasb3of, nuasb3of,andoutf;
	GINT64 stack_count;
	int diagh;
	// ================== entry in the proces
	void Init( );
	inline int AddCell_Of(uint32_t cell, int bit) {
		register int s= C_stack[cell];
		if (stack_count.u16[s] > 5) return 0;
		stack_count.u16[s]++;
		if (stack_count.u16[s] > 5) {
			s =~( 07007007 << (3 * s));// mask
			wua &= s;
			wactive0 &= s;
		}
		nmiss--;
		known_b3 |= bit;
		return 1;
	}
	uint32_t IsMultiple(int bf);
	int ShrinkUas1();
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs();
	void Go_Critical();
	void CriticalLoop();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini(int M, int mask);
	void Go_Subcritical();
	void Go_SubcriticalMiniRow();
	void Go_SubcriticalMiniRow_End(int stack);
	//===================== process not critical
	void ShrinkUasOfB3();
	void Go_miss1_b3();
	void Go_miss2_b3();
	void Go_miss3_b3();

};
/*
entry 92
maxindex= 983
maxn5= 51516
maxn6= 237762
maxdet5= 261
maxdet6= 2004
*/

struct BANDS_AB {// handling bands 12 in A B mode

	G17TMORE moreuas_AB, moreuas_AB_small, moreuas_AB_big;
	MORE32 moreuas_b3;
	uint32_t  mode_ab, ia, ib,myuab,
		indd,indf,ncluesbandb,
		nxy_filt1,nvalid;// bufferstoring AB/XY passing filt1
	GINT64  stack_count, stack_countf;
	uint32_t *mybv5;
	STD_B1_2 * mybb,*myba;
	//=============== bands B order at start AB
	uint64_t tuaB[TUA64_12SIZE];// valid uas mode AB ordered
	uint32_t tiBc_ideb[257],
		tiBc_pat[256],// pattern ib B
		tiBc_kill[256],// killer in A (& of all uas)
		ntuaB, ntiBc;
	BF128 v128B, v128_mult, vc128A[27], vc128B[27],
		v128B_y6[MAXNIND6],v128B_x5[MAXNIND5];
	VECT256 v256B, v256_mult, vc256A[27], vc256B[27],
		v256B_y6[MAXNIND6], v256B_x5[MAXNIND5];

	//========== bands AB uas at index 3
	uint32_t tybf[MAXNIND6],// sub table 6 clues band B 3Y
		txbf[MAXNIND5];// same for 3X 5 clues
	uint32_t i3, i5, iy3, // index in main loops  
		nx3, ny3;// number of valid 5/6 clues in the 3X3Y
	uint32_t  bf3, bf5,ybf3;
	uint64_t b1b2_3x3y_bf,b1b2_xy_bf;
	uint32_t   ntib3c;

	int diag, diagbug;

	MINCOUNT smin;
	//========== tclues for valid XY 
	uint32_t tclues[40];// mini 25+band a
	int ncluesa, nclues;
	uint32_t  bfA,bfB,bf1,bf2;
	BF128 bands_active_pairs, bands_active_triplets,
		valid_vect;
	//==================== current band 3 to process
	int cur_ib;
	uint32_t tcluesb12[20], ncluesb3x;
	uint32_t   nmiss;
	uint32_t uasb3_1[2000], uasb3_2[2000],uas_in[2000], 
		nuasb3_1, nuasb3_2,nuas_in,b3_andout;

	void Go(STD_B1_2 & ba, STD_B1_2 & bb, int i, int mode);
	void Go3X3Y();
	void Go3X3Y128();
	void Go3X3Y256();
	void DebugBuildUasAB(int mode=0);
	void DebugAdd12(uint32_t itemp, uint64_t uab12, GINT64  w);
	void CleanTempXY();
	void CleanTempXY_2(int itemp);
	int CTXY_Unique(int itemp);
	void CTXY_SetupValidVect();
	int CTXY_SetupB3MinplusCurib();
	void CTXY_Gob3Curib();
	void B12InTclues() {
		nclues = 0;
		stack_count.u64 = 0;
		uint64_t w= b1b2_xy_bf;
		uint32_t xcell, cell;
		while (bitscanforward64(xcell, w)) {
			w ^= (uint64_t)1 << xcell;
			cell = From_128_To_81[xcell];
			stack_count.u16[C_stack[cell]]++;
			tclues[nclues++] = cell;
		}
	}
	void CleanValid();
	void FinalCheckB3(uint32_t bfb3);
	void Fout(uint32_t bfb3);
	void MergeUasExpand();
	void ExpandBand3(uint32_t *tua, uint32_t nua);
	void Debug_If_Of_b3();
};
struct TU_LOCK {// current pointer to buffer for bands expansion
	XINDEX3 * px_index;// pointer to next index 3 (b1,b2,dummy b3)
	uint32_t * px_5, *px_6;// pointer to valid band5/6 next
	// store index end band 1 to avoid rebuilding
	XINDEX3 * rpx_index;// pointer to next index 3 (b1,b2,dummy b3)
	uint32_t * rpx_5, *rpx_6;// pointer to valid band5/6 next

	void LockExpand(int nindex, int n5, int n6) {
		px_index = &px_index[nindex];
		px_5 =&px_5[ n5];
		px_6 = &px_6[n6];
	}
	void InitBuf();
	inline void Store1() {
		rpx_index= px_index;
		rpx_5 = px_5;
		rpx_6 = px_6;
	}
	inline void Restore1() {
		px_index = rpx_index;
		px_5 = rpx_5;
		px_6 = rpx_6;

	}

};

struct TU_GUAN {// GUAN process (used GUA all kinds)
	struct GUAN3X3Y  {// pointing to the guar sub vector id,if
		int i, 	indd, indf;//  i=id/128  id% 128  if%128 
	}tuaindex[128];
	GUAN tguan[256], guanw, tgua3x[128],tgua3x3y[128];
	//GUAN tguan4[81], tguan4_3x[81];
	uint32_t nguan,ngua3x3y, ngua3x;
	uint32_t ng2, ng3;// debugging only
	uint64_t  guabuf[15000], *pguabuf;
	uint64_t  guabuf3x[12000], *pguabuf3x;
	uint64_t  guabufr[10000], *pguabufr;
//	uint32_t nuar,nuar23x;
	BF128 v3x3y, vcells3x3y[54], vxy, vxy_raw,xy23sockets,
		v4_3x; 
	void AddGuan(uint64_t *t, uint32_t n, uint32_t cbf,
		int32_t dbf, uint32_t ind) {
		if (nguan < 256) {
			pguabuf = &pguabuf[n];// lock space
			tguan[nguan++].Enter(t, n, cbf, dbf, ind);
		}
		//tguan[nguan - 1].Debug1Guan(nguan - 1);
	}

	void Init(){
		pguabuf = guabuf;// reinit gua buffer use
		nguan = 0;//and guan table
	}
	void Build_Guas3X();// after 3 clues A reduce tables

	void Build_Guas3X3Y();// after 3 clues A reduce tables
	void Debug3X3Y();
	void Debug1() {
		for (uint32_t i = 0; i < nguan; i++)
			tguan[i].Debug1Guan(i);
	}
	void Debug3X() {
		cout << "test Debug3X nguan=" << ngua3x << " n3x=" << ngua3x << endl;
		for (uint32_t i = 0; i < ngua3x; i++)
			tgua3x[i].Debug1Guan(i);
	}
	void Debugvxy() {
		cout << " guan actifs en vB" << endl;
		for (uint32_t i = 0; i < 128; i++)
			if(vxy.On(i))
				tgua3x3y[i].Debug1Guan(i);
	}
};
  
struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	BF128 p17diag;// known 17 pattern for tests
	int b3lim, debug17,debug17_check,
		diag, diagbug,debugb3, aigstop,
		npuz, a_17_found_here;
	uint32_t	iband1,iband2;
	G17B3HANDLER g17hh0;

	MORE32 moreb;
	//______sockets common to  all bands 3  
	BF128 isguasocket2all, isguasocket3all;
	

	//=====================process
	void GoM10();// end preparation for a given band1+band2+ table band3 pass 656 566
	void GoM10Uas();// end preparation for a given band1+band2+ table band3 pass 656 566
	void GoM10_guas_four_columns();// end preparation for a given band1+band2+ table band3 pass 656 566
	void Go();
			
	//================ debugging code
	inline uint32_t k17x(int ix) {	return p17diag.bf.u32[ix];	}
	void DebugGetPuz(const char * p) {
		p17diag.SetAll_0();
		for (int i = 0; i < 81; i++)
			if (p[i] != '.')p17diag.Set_c(i);

		cout <<"this is a "
			<<_popcnt32(p17diag.bf.u32[0])
			<< _popcnt32(p17diag.bf.u32[1])
			<< _popcnt32(p17diag.bf.u32[2]) 
			<<" pattern for the expected puzzle"<< endl;
	}
	int DebugK17M10();
	void GodebugInit(int mode);
	int GodebugFindKnown17();
	int GodebugCheckUas(const char * lib);
};