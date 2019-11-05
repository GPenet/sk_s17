//#define VALc15 13438
//#define VALc17 1032
//#define TESTXY 0
//2536435903110145  i1=2 i2=35
//#define TESTXY2 0
//22521159232789761
//536435903110145
//#define DEBUGLEVEL 10  nb12=3461507
//123456789456789123789132564268591437341627895597843216634278951815964372972315648;
//1....6...4......2.....3.5.....59........2...........16.....89.1..5.......72......; 189; 339; 341; 6ua < 12 to add   clean


//___ start process expand bands collect uas guas ...
void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	if (aigstop)return;
	const char * diagband = "281364975374591862569872341"; // la bonne
	//const char * diagband =   "349628517586371492712945863"; // le bug
	const char * diagband3 = "712648593845913627936725418";
	const char * diagpuz = ".2..........1...36..8..7......3.......4...86......2.......485...........93.....1.";
	diag = diagbug = 0;
	//if (!strcmp(diagband, myband2.band)) {
	//if (genb12.nb12 != 2694160) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//cout << myband2.band << "band2 id=" << myband2.i416<<" nb12="<< genb12.nb12 << endl;
	if (0){//(!strcmp(diagband, myband2.band)) {
		if (strcmp(diagband, myband2.band)) return;
		cout << myband2.band << "go band2 id=" << myband2.i416 << " nb12="<< genb12.nb12 << endl;
		DebugGetPuz(diagpuz);
		diagbug = 2;
	}
	else {
		//if (1) return;
		if (genb12.nb12 == sgo.vx[6]) diagbug = 1;
		//cout << myband2.band << "go band2 id=" << myband2.i416 << " nb12=" << genb12.nb12 << endl;
	}
	p_cpt2g[0] ++;
	p_cpt2g[1] +=genb12.nband3;
	if (genb12.nband3 > p_cpt2g[23])	p_cpt2g[23] = genb12.nband3;

	if (g17b.debug17) diag = diagbug = g17b.debug17;
	//______ true start band 2 expand
	tulock.Restore1();	
	myband2.ExpandBand();// expand band2
	if (!(myband1.nmybv5 | myband2.nmybv5)) return; // no 656 no 566
	tulock.LockExpand(myband2.nmyi3, myband2.nmybv5, myband2.nmybv6);
	if (diagbug == 2) {
		GodebugFindKnown17();
		//myband1.DebugIndex();
	}
	//myband1.DebugIndex();
	//for (int i = 0; i < 100; i++)
	//	cout << Char27out(myband1.mybv5[i]) << " v5 i=" << i << endl;
	//if (1) return;
	
	/*
	//================= expand 6 for bands 3

	nb3 = genb12.nband3;
	debugb3 = -1;
	for (int ib3 = 0; ib3 < nb3; ib3++) {// lock only 6 clues
		genb12.bands3[ib3].ExpandBand();// expand band3
		tulock.LockExpand(genb12.bands3[ib3].nmyi3, 0, genb12.bands3[ib3].nmybv6);
		if (diagbug) {
			if (strcmp(diagband3, genb12.bands3[ib3].band) == 0) {
				genb12.bands3[ib3].DebugExp();
				cout << "this is the expected band 3 ib3="<<ib3 << " out of nb3="<<nb3 << endl;
				debugb3 = ib3;
			}
		}
	}
	*/

	GoM10Uas();//expand bands 3  collect UAs 
	if (g17b.debug17) if (DebugK17M10()) return;

	p_cpt2g[18] += genuasb12.nua;
	p_cpt2g[19] += genb12.ntua2;
	p_cpt2g[20] += genb12.ntua3;
	if (g17b.debug17) 			GodebugInit(0);
	Go();// standard entry point for all 
}
void G17B::GoM10Uas() {
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	//cout << " call init gen test sockest 2 4 bands 12" << endl;
	if (genuasb12.Initgen()) return;
	genb12.BuildGang9x3();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	tuguan.Init();
	genb12.SecondSockets2Setup();// collect GUA2s 
	tuguan.ng2 = tuguan.nguan;

	//tuguan.tguan[13].Debug1Guan(13, 1);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	genb12.SecondSockets3Setup();// collect GUA3s 
	tuguan.ng3 = tuguan.nguan- tuguan.ng2;
	// setupsockets common to all band3
	isguasocket2all = genb12.bands3[0].guas.isguasocket2;
	isguasocket3all = genb12.bands3[0].guas.isguasocket3;
	for (int ib3 = 1; ib3 < genb12.nband3; ib3++) {
		isguasocket2all &= genb12.bands3[ib3].guas.isguasocket2;
		isguasocket3all &= genb12.bands3[ib3].guas.isguasocket3;
	}
	// ===== add new guas four cells four columns 2 boxes 
	//tuguan.BuildCellKillersVector();
}

//____ launch {band A ; band B}  needed
void G17B::Go(){// start loops on 3_5 clues 
	if (aigstop)return;
	if (g17b.debug17)cout << "entry g17b.go" << endl;
	if ((!sgo.vx[4]) && myband1.nmybv5) {// A is band1
		bands_ab.Go(myband1, myband2, 0, 1);
	}
	if (aigstop)return;
	if (myband2.nmybv5) {
		// prepare zh2b_g with bands 2 first
		//memcpy(zh2b_g.puz0, myband2.band0, sizeof myband2.band0);
		//memcpy(&zh2b_g.puz0[27], myband1.band0, sizeof myband1.band0);
		//zh2b_g.GetBands(myband2.gangster, myband1.gangster);// set sol/pm
		//_______________________________________________________
		bands_ab.Go(myband2, myband1, 1, 1);
	}
}
#define XCHUNK 100
#define YCHUNK 200

//_______________ processing bandA bandB (5clues 6 clues)
void BANDS_AB::Go(STD_B1_2 & ba, STD_B1_2 & bb, int i, int mode) {
	if (g17b.debug17)cout << "entry BANDS_AB.go i=" <<i << endl;
	ia = i; ib = 1 - i;
	p_cpt2g[4] += ba.nmybv5*bb.nmybv6/100;
	if (g17b.debug17) {
		if (_popcnt32(g17b.p17diag.bf.u32[ia]) != 5) return;
		cout << "\t\t\tgo a b  ia=" << ia << endl;
	}
	mode_ab = mode;// 1 if must be 5 clues 
	mybb = &bb;
	myba = &ba;
	uint32_t ni3 = ba.nmyi3 - 1, nyi3 = bb.nmyi3 - 1;
	mybv5 = ba.mybv5;
	//myi3 = ba.myi3;
	for (i3 = 0; i3 < ni3; i3++) {
		// store the x 5 clues in a fixed place 
		indd = myba->myi3[i3].ideb5;
		indf = myba->myi3[i3 + 1].ideb5;
		if (indd == indf)continue;// no X 5 clues here
		nx3 = indf - indd;
		memcpy(txbf, &myba->mybv5[indd], 4 * nx3);
		XINDEX3 wi3 = myba->myi3[i3];
		bf3 = wi3.cellsbf;
		tuguan.Build_Guas3X();
		for (iy3 = 0; iy3 < nyi3; iy3++) {
			// store the y 6 clues in a fixed place 
			uint32_t yindd = mybb->myi3[iy3].ideb,
				yindf = mybb->myi3[iy3 + 1].ideb;
			ny3 = yindf - yindd;
			//cout << "ny3=" << ny3 << endl;
			memcpy(tybf, &mybb->mybv6[yindd], 4 * ny3);
			XINDEX3 wiy3 = mybb->myi3[iy3];
			ybf3 = wiy3.cellsbf;
			Go3X3Y();
		}
	}

}
void TU_GUAN::Build_Guas3X() {// guas filter => already killed/forced plus sub table
	pguabuf3x = guabuf3x;
	register uint64_t F = bands_ab.bf3;
	if (bands_ab.ia)F <<= 32;
	ngua3x =  0;
	{
		// go through the first table to extract still valid
		for (uint32_t ig = 0; ig < nguan; ig++) {
			GUAN & g = tguan[ig];
			if (g.killer&F)continue;
			GUAN gn = g;// catch permanent data
			gn.tua = pguabuf3x;
			gn.nua = 0;
			register uint64_t killer = BIT_SET_2X;
			for (uint32_t iua = 0; iua < g.nua; iua++) {
				register uint64_t ua = g.tua[iua];
				if (ua & F)continue;
				killer &= ua;
				pguabuf3x[gn.nua++] = ua;
				if (gn.nua > 20) break;// optimization  
			}
			if (!gn.nua)continue;
			gn.killer = killer;
			tgua3x[ngua3x++] = gn;
			pguabuf3x += gn.nua;
			if (ngua3x > 127) break; // limit 128 here
		}
	}

	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		int blocib3 = ib3 >> 7, irib3 = ib3 & 127;
		STD_B3 & b = genb12.bands3[ib3];
		STD_B3::GUAs &g = b.guas;
		memset(&b.g3x, 0, sizeof b.g3x);
		for (uint32_t ig = 0; ig < ngua3x; ig++) {
			GUAN gu = tgua3x[ig];
			int i81 = gu.i81, stack = i81 / 27;
			if ((gu.ncol == 2 && g.isguasocket2.On_c(i81)) ||
				(gu.ncol == 3 && g.isguasocket3.On_c(i81))) {
				b.g3x.g23.Set(ig);
			}
			else if (gu.ncol == 2 && g.isguasocket2_46.On_c(i81)) {
				b.g3x.g46.Set(ig);
			}
			b.g3x.gall = b.g3x.g46 | b.g3x.g23;
		}
	}
}
void TU_GUAN::Debug3X() {
	cout << "test Debug3X nguan=" << nguan << " n3x=" << ngua3x << endl;
	//for (uint32_t i = 0; i < ngua3x ; i++)tgua3x[i].Debug1Guan(i);
	v4_3x.Print("v4_3x status");
}

//____________ processing a 3X3Y set
void BANDS_AB::Go3X3Y() {
	moreuas_AB.Init();
	moreuas_AB_small.Init();
	moreuas_AB_big.Init();
	//indd = myba->myi3[i3].ideb5;
	//indf = myba->myi3[i3 + 1].ideb5;
	//if (indd == indf)return;// no X 5 clues here
	//wi3 = myba->myi3[i3];
	p_cpt2g[2]++;
	p_cpt2g[3] += nx3 * ny3;
	  
	if (ia) b1b2_3x3y_bf = ((uint64_t)bf3 << 32) | ybf3;
	else b1b2_3x3y_bf = ((uint64_t)ybf3 << 32) | bf3;
	if (g17b.debug17) {
		if(b1b2_3x3y_bf & (~g17b.p17diag.bf.u64[0]))return;
		cout <<Char2Xout(b1b2_3x3y_bf)<<" this is a valid 3X3Y"<<endl;
		tuguan.Debug3X();

	}

	tuguan.Build_Guas3X3Y();//<<<<<<<<<<<<<<<<<<<
	//tuguan.Debug3X3Y();

	{	// Sort  on B increasing number of cells
		register uint64_t F = b1b2_3x3y_bf;
		ntuaB = 0;
		for (uint32_t i = 0; i < genuasb12.nua; i++) {
			uint64_t R = genuasb12.tua[i] & BIT_SET_2X;
			if (R&F) continue;// hit by 3X3Y
			uint64_t Ra = R & BIT_SET_27, Rb = R >> 32;
			if (ia) { Rb = R & BIT_SET_27; Ra = R >> 32; }
			Rb <<= 32;// now high bits
			uint64_t cc = _popcnt64(Rb);
			R = Ra | Rb | (cc << 59);
			for (uint64_t j = 0; j < ntuaB; j++) {
			// Insert in tt ordered on B
				if (R < tuaB[j]) {// insert here
					for (uint64_t k = ntuaB; k > j; k--)	tuaB[k] = tuaB[k - 1];
					tuaB[j] = R;
					ntuaB++;
					goto nextuaold;
				}
			}
			tuaB[ntuaB++] = R;// added end of table
		nextuaold:;
		}
	}
	//p_cpt2g[32] += ntuaB;
	{// build index of the first 256 entries
		ntiBc = 0;
		uint32_t oldb = 0;
		uint64_t i = 0;
		for (; i < ntuaB; i++) {
			uint64_t R = tuaB[i];
			uint32_t Ra = (uint32_t)(R & BIT_SET_27),
				Rb = (uint32_t)(R >> 32)& BIT_SET_27;
			if (Rb == oldb) {
				tiBc_kill[ntiBc - 1] &= Ra;
				continue;
			}
			if (ntiBc >= 256) break;// limit for matrix test
			oldb = Rb;
			tiBc_pat[ntiBc] = Rb;
			tiBc_kill[ntiBc] = Ra;
			tiBc_ideb[ntiBc++] = (uint32_t)i;
		}
		tiBc_ideb[ntiBc] = (uint32_t)i;// last index
	}
	p_cpt2g[24 + (ntiBc + 127) / 128]++;
	if (!ntiBc) return; //just in case
	nxy_filt1 = 0;
	//cout << "i3=" << i3 << " iy3=" << iy3 << " ntiBc="
	//	<< " nx3="<<nx3<<" ntiBc="<< ntiBc << endl;
	if (ntiBc <= 128)Go3X3Y128();
	else Go3X3Y256();
	if (nxy_filt1)CleanTempXY();// final clean 3X3Y

}
void TU_GUAN::Build_Guas3X3Y() {// guas filter => already killed/forced plus sub table
	pguabufr = guabufr;
	ngua3x3y = ngua3x;
	v3x3y.SetAll_0();
	{
		register uint64_t F = bands_ab.b1b2_3x3y_bf;
		// go through the 3x table to extract still valid
		for (uint32_t ig = 0; ig < ngua3x; ig++) {
			GUAN & g = tgua3x[ig];
			if (g.killer&F)continue;
			GUAN gn = g;// catch permanent data
			gn.tua = pguabufr;
			gn.nua = 0;
			register uint64_t killer = BIT_SET_2X;
			for (uint32_t iua = 0; iua < g.nua; iua++) {
				register uint64_t ua = g.tua[iua];
				if (ua & F)continue;
				killer &= ua;
				pguabufr[gn.nua++] = ua;
				if (gn.nua > 15) break;// optimization  
			}
			if (!gn.nua)continue;
			gn.killer = killer;
			tgua3x3y[ig] = gn;
			v3x3y.Set(ig);
			pguabufr += gn.nua;
		}
	}
}
void TU_GUAN::Debug3X3Y(){
	cout << "test build nguan=" << nguan<< " n3x3y=" << ngua3x3y<< endl;
	cout << Char2Xout(bands_ab.b1b2_3x3y_bf) << " 3x3ybf" << endl;
	//Debug1();
	v3x3y.Print(" v3x3y");
	//for(uint32_t i=0;i< ngua3x3y;i++)tgua3x3y[i].Debug1Guan(i);
	//for (int i = 0; i < 54; i++)
		//cout << Char64out(vcells3x3y[i].bf.u64[0]) << " vc i=" << i << endl;
	cout << "\n\nfirst band status" << endl;
	STD_B3 & b = genb12.bands3[0];
	b.g3x.g23.Print("b.g3x.g23");
	b.g3x.g46.Print("b.g3x.g46");
	b.g3x.gall.Print("b.g3x.gall");
}
void BANDS_AB::Go3X3Y128() {
	v128B = maskLSB[ntiBc];
	memset(vc128A, 255, sizeof vc128A);// all bits to 1
	memset(vc128B, 255, sizeof vc128B);// all bits to 1
	{	// build vectors A (killer filter) and vector B (matrix filter)
		uint32_t cc;
		for (uint32_t i = 0; i < ntiBc; i++) {
			register uint32_t Rw = tiBc_kill[i];
			while (bitscanforward(cc, Rw)) {// look for  possible cells
				register uint32_t bit2 = 1 << cc;
				Rw ^= bit2;// clear bit
				vc128A[cc].Clear(i);
			}
			Rw = tiBc_pat[i];
			while (bitscanforward(cc, Rw)) {// look for  possible cells
				register uint32_t bit2 = 1 << cc;
				Rw ^= bit2;// clear bit
				vc128B[cc].Clear(i);
			}
		}
	}

	{// _____ apply vcellsb to valid band B 6
		register uint32_t* R6 = tybf;
		register BF128 * Rv = v128B_y6;
		for (uint32_t i6 = 0; i6 < ny3; i6++, R6++, Rv++) {
			*Rv = v128B;
			register uint32_t bf = (*R6)^ ybf3;// 3 more clues
			int ccb;
			while (bitscanforward(ccb, bf)) {
				bf ^= 1 << ccb; //clear bit
				(*Rv)&=vc128B[ccb];
			}
		}
	}
	//============== setup the v128_mult vector
	v128_mult.SetAll_0();
	for (uint32_t i = 0; i < ntiBc; i++) {
		// set bit more than one A ua
		if ((tiBc_ideb[i + 1] - tiBc_ideb[i]) > 1)
			v128_mult.Set(i);
	}
	// _________________loop on remaining 2 clues
	for (i5 = 0; i5 < nx3; i5++) {
		bf5 = txbf[i5];
		ncluesa = 0;
		BitsInTable32((int *)tclues, ncluesa, bf5);
		{ //build the band B X vector
			BF128 vw = v128B;
			uint32_t bfplus = bf5 ^ bf3;
			register uint32_t w = bfplus;
			{// find the final 128 small vector
				int cca;// this is bandA
				while (bitscanforward(cca, w)) {
					w ^= 1 << cca; //clear bit
					vw &= vc128A[cca];
				}
			}// now vw has possible b pats not killed
			{// check index multiple uas
				BF128  w = vw;
				register uint32_t F = bf5;
				w &= v128_mult; // patterns to chek
				int t[128], nt = w.Table128(t);
				for (int i = 0; i < nt; i++) {
					int ind = t[i];// index of the pattern in tiBc
					uint32_t iud = tiBc_ideb[ind], iuf = tiBc_ideb[ind + 1];
					for (uint32_t iua = iud; iua < iuf; iua++) {
						register uint32_t Ra = tuaB[iua] & BIT_SET_27;
						if (!(Ra & F)) {//
							goto nextindex; // first is ok
						}
					}
					vw.Clear(ind);
				nextindex:;
				}
			}// now vw all valid entries in tiBc
			v128B_x5[i5]=vw;
		}

		if (g17b.aigstop)return;
	}// end lop 3 to 5 clues in band A


	// now the matrix process with a vector size 128
	uint32_t xch_bf[XCHUNK], ych_bf[YCHUNK];
	BF128 vxch[XCHUNK], vych[YCHUNK];
	uint32_t  ideby = 0, iendy = YCHUNK;
	if (iendy > ny3)iendy = ny3;
	while (ideby < ny3) { //Y chunk
		register int ny = iendy - ideby;
		// collect date for Y chunk
		memcpy(ych_bf, &tybf[ideby], 4 * ny);
		memcpy(vych, &v128B_y6[ideby], 16 * ny);
		uint32_t idebx = 0, iendx = XCHUNK;
		if (iendx > nx3)iendx = nx3;

		while (idebx < nx3) {// X chunk
			register int nx = iendx - idebx;
			// collect date forXchunk
			memcpy(xch_bf, &txbf[idebx], 4 * nx);
			memcpy(vxch, &v128B_x5[idebx], 16 * nx);
			register BF128 * Ry = vych;
			for (int iy = 0; iy < ny; iy++, Ry++) {
				BF128 y128 = *Ry;
				register BF128 * Ra = vxch;
				for (int ix = 0; ix < nx; ix++, Ra++) {
					BF128 x128 = *Ra;
					if ((x128 & y128).isNotEmpty()) continue;
					GINT64 & w = tempXY[nxy_filt1++];
					w.u32[0] = xch_bf[ix];
					w.u32[1] = ych_bf[iy];
				}
			}// end of a matrix XY
			if (nxy_filt1 > 1000) {
				//cout << "call clean nfilt1=" << nxy_filt1 << endl;
				//nxy_filt1 = 0;
				CleanTempXY();
			}
			idebx = iendx; iendx += XCHUNK;
			if (iendx > nx3)iendx = nx3;
		}
		ideby = iendy; iendy += YCHUNK;
		if (iendy > ny3)iendy = ny3;
	}

}
void BANDS_AB::Go3X3Y256() {
	v256B.Init(ntiBc);
	memset(vc256A, 255, sizeof vc256A);// all bits to 1
	memset(vc256B, 255, sizeof vc256B);// all bits to 1
	{	// build vectors A (killer filter) and vector B (matrix filter)
		uint32_t cc;
		for (uint32_t i = 0; i < ntiBc; i++) {
			uint32_t bloc = i >> 7, ir = i & 127;
			register uint32_t Rw = tiBc_kill[i];
			while (bitscanforward(cc, Rw)) {// look for  possible cells
				register uint32_t bit2 = 1 << cc;
				Rw ^= bit2;// clear bit
				vc256A[cc].Clear(bloc,ir);
			}
			Rw = tiBc_pat[i];
			while (bitscanforward(cc, Rw)) {// look for  possible cells
				register uint32_t bit2 = 1 << cc;
				Rw ^= bit2;// clear bit
				vc256B[cc].Clear(bloc,ir);
			}
		}
	}

	{// _____ apply vcellsb to valid band B 6
		register uint32_t* R6 = tybf;
		register VECT256 * Rv = v256B_y6;
		for (uint32_t i6 = 0; i6 < ny3; i6++, R6++, Rv++) {
			*Rv = v256B;
			register uint32_t bf = (*R6) ^ ybf3;// 3 more clues
			int ccb;
			while (bitscanforward(ccb, bf)) {
				bf ^= 1 << ccb; //clear bit
				(*Rv).And( vc256B[ccb]);
			}
		}
	}
	//============== setup the v256_mult vector
	memset(&v256_mult, 0, sizeof v256_mult);
	for (uint32_t i = 0; i < ntiBc; i++) {
		// set bit more than one A ua
		if ((tiBc_ideb[i + 1] - tiBc_ideb[i]) > 1)
			v256_mult.Setx(i);
	}
	// _________________loop on remaining 2 clues
	for (i5 = 0; i5 < nx3; i5++) {
		bf5 = txbf[i5];
		ncluesa = 0;
		BitsInTable32((int *)tclues, ncluesa, bf5);
		{ //build the band B X vector
			VECT256 vw = v256B;
			uint32_t bfplus = bf5 ^ bf3;
			register uint32_t w = bfplus;
			{// find the final 128 small vector
				int cca;// this is bandA
				while (bitscanforward(cca, w)) {
					w ^= 1 << cca; //clear bit
					vw.And(vc256A[cca]);
				}
			}// now vw has possible b pats not killed
			{// check index multiple uas
				VECT256  w = vw;
				register uint32_t F = bf5;
				w.And( v256_mult); // patterns to chek
				int t[256], nt;			w.Table(t,nt);
				for (int i = 0; i < nt; i++) {
					int ind = t[i];// index of the pattern in tiBc
					uint32_t iud = tiBc_ideb[ind], iuf = tiBc_ideb[ind + 1];
					for (uint32_t iua = iud; iua < iuf; iua++) {
						register uint32_t Ra = tuaB[iua] & BIT_SET_27;
						if (!(Ra & F)) {//
							goto nextindex; // first is ok
						}
					}
					vw.Clearx(ind);
				nextindex:;
				}
			}// now vw all valid entries in tiBc
			v256B_x5[i5] = vw;
		}

		if (g17b.aigstop)return;
	}// end lop 3 to 5 clues in band A


	// now the matrix process with a vector size 256
	uint32_t xch_bf[XCHUNK], ych_bf[YCHUNK];
	VECT256 vxch[XCHUNK], vych[YCHUNK];
	uint32_t  ideby = 0, iendy = YCHUNK;
	if (iendy > ny3)iendy = ny3;
	while (ideby < ny3) { //Y chunk
		register int ny = iendy - ideby;
		// collect date for Y chunk
		memcpy(ych_bf, &tybf[ideby], 4 * ny);
		memcpy(vych, &v256B_y6[ideby], 32 * ny);
		uint32_t idebx = 0, iendx = XCHUNK;
		if (iendx > nx3)iendx = nx3;

		while (idebx < nx3) {// X chunk
			register int nx = iendx - idebx;
			// collect date forXchunk
			memcpy(xch_bf, &txbf[idebx], 4 * nx);
			memcpy(vxch, &v256B_x5[idebx], 32 * nx);
			register VECT256 * Ry = vych;
			for (int iy = 0; iy < ny; iy++, Ry++) {
				VECT256 y256 = *Ry;
				register VECT256 * Ra = vxch;
				for (int ix = 0; ix < nx; ix++, Ra++) {
					VECT256 x256 = *Ra;
					if ((x256.v[0] & y256.v[0]).isNotEmpty()) continue;
					if ((x256.v[1] & y256.v[1]).isNotEmpty()) continue;
					GINT64 & w = tempXY[nxy_filt1++];
					w.u32[0] = xch_bf[ix];
					w.u32[1] = ych_bf[iy];
				}
			}// end of a matrix XY
			if (nxy_filt1 > 1000) {
				//cout << "call clean nfilt1=" << nxy_filt1 << endl;
				//nxy_filt1 = 0;
				CleanTempXY();
			}
			idebx = iendx; iendx += XCHUNK;
			if (iendx > nx3)iendx = nx3;
		}
		ideby = iendy; iendy += YCHUNK;
		if (iendy > ny3)iendy = ny3;
	}

}
	//__________________________ cleaning after first filter >= 5 bands
void BANDS_AB::CleanTempXY() {//GINT64 tempXY[15000];
	p_cpt2g[5]+= nxy_filt1;
	nvalid = 0;
	for (uint32_t itemp = 0; itemp < nxy_filt1; itemp++) {
		GINT64  w = tempXY[itemp];
		bfA = w.u32[0]; bfB = w.u32[1];
		if (ia) { bf1 = bfB; bf2 = bfA; }// work now in bf1 bf2 mode;
		else { bf1 = bfA; bf2 = bfB; }
		b1b2_xy_bf = ((uint64_t)bf2 << 32) | bf1;

		// check the "more" table if available
		if (moreuas_AB_small.Check(b1b2_xy_bf))continue;
		if (moreuas_AB.Check(b1b2_xy_bf)) continue;
		if (moreuas_AB_big.Check(b1b2_xy_bf)) continue;

		if (g17b.debug17) {//skip if not ok
			if (b1b2_xy_bf != g17b.p17diag.bf.u64[0]) continue;
			cout << Char2Xout(b1b2_xy_bf) << "we have the right XY" << endl;
		}
		B12InTclues();
		if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
			stack_count.u16[2] > 6) continue;
		if (g17b.diagbug == 2) {
			if (w.u32[0] == g17b.k17x(ia) && w.u32[1] == g17b.k17x(ib))
				cout << Char27out(w.u32[1]) << " right band B try check valid" << endl;
			//else continue;//<< 
		}
		if (genb12.nband3 < 5)CleanTempXY_2(itemp);
		else 		{
			p_cpt2g[6]++;
			if (CTXY_Unique(itemp)) continue;;
			p_cpt2g[7]++;
			valid_b12[nvalid++] = b1b2_xy_bf;
		}
	}
	nxy_filt1 = 0;
	if (nvalid)CleanValid();
}
void BANDS_AB::CleanValid() {
	for (uint32_t iv = 0; iv < nvalid; iv++) {
		b1b2_xy_bf = valid_b12[iv];
		B12InTclues();
		if (zhou[0].PartialInitSearch17(tclues, nclues))return;// would be  bug
		CTXY_SetupValidVect();
		for (cur_ib = 0; cur_ib < genb12.nband3; cur_ib++) {
			if (CTXY_SetupB3MinplusCurib()) continue;
			CTXY_Gob3Curib();
		}
	}
}
	//__________________________ cleaning after first filter < 5 bands
void BANDS_AB::CleanTempXY_2(int itemp) {
	p_cpt2g[32]++;
	CTXY_SetupValidVect();
	int aig = 1;
	for (cur_ib = 0; cur_ib < genb12.nband3; cur_ib++) {
		if (CTXY_SetupB3MinplusCurib()) continue;
		if (aig) {
			p_cpt2g[6]++;
			if (CTXY_Unique(itemp)) return;// only once
			p_cpt2g[7]++;
			if (zhou[0].PartialInitSearch17(tclues, nclues))return;// would be  bug
		}
		aig = 0;
		CTXY_Gob3Curib();
	}
	

}int BANDS_AB::CTXY_Unique(int itemp) {
	//_____ check for unique solution bands A;B
	register uint64_t myua = zh2b[0].ValidXY(tclues, nclues, 0);
	if (myua) {//___ if not do the best with the UA
		//if (g17b.debug17 ) 		cout << Char2Xout(myua) << "uaret band b multiple" << endl;
		uint64_t cc64 = _popcnt64(myua&BIT_SET_2X);
		if (cc64 < 12) {// this should never be check for a bug
			DebugAdd12(itemp, myua, tempXY[itemp]);
			g17b.aigstop = 1;
			return 1;
		}
		if (cc64 < 18)			moreuas_AB_small.Add(myua);
		else if (cc64 < 21)			moreuas_AB.Add(myua);
		else moreuas_AB_big.Add(myua);
		if (cc64 < 18) {//  update the UA table (for next 3 clues band A)
			//cout << Char2Xout(uab12) << "to add" << endl;
			genuasb12.AddUACheck(myua | (cc64 << 59));
			p_cpt2g[31]++;
		}
		return 1;
	}// end not a unique solution

	return 0;
}
		// clean tasks (not in the same sequence in both branches)
void BANDS_AB::CTXY_SetupValidVect() {
	valid_vect.SetAll_0();
	{		//_____________ build valid guas
		register uint64_t F = b1b2_xy_bf;
		for (uint32_t ig = 0; ig < tuguan.ngua3x3y; ig++) {
			if (tuguan.v3x3y.Off(ig))continue; // killed in 3x3y
			GUAN gu = tuguan.tgua3x3y[ig];
			if (gu.killer & F) continue;
			if (gu.nua < 2) {
				valid_vect.Set(ig);
				continue;// if 1, killer is enough
			}
			for (uint32_t iua = 0; iua < gu.nua; iua++) {
				if (gu.tua[iua] & F) continue;
				valid_vect.Set(ig);
				break; // first not hit is enough
			}
		}
	}

}
int BANDS_AB::CTXY_SetupB3MinplusCurib() {
	STD_B3 & b = genb12.bands3[cur_ib];
	b.SetUpMincountxy(valid_vect);
	if (b.smin.mincount > 6)return 1;
	// build and check start count per stack
	stack_countf.u64 = stack_count.u64 + b.smin.Count_per_stack().u64;
	if (stack_countf.u16[0] > 6 || stack_countf.u16[1] > 6 ||
		stack_countf.u16[2] > 6) return 1; // not ok

	// build free staks pattern for outfield uas
	uint32_t fstk = BIT_SET_27;
	if (stack_countf.u16[0] == 6) fstk ^= 07007007;
	if (stack_countf.u16[1] == 6) fstk ^= 070070070;
	if (stack_countf.u16[2] == 6) fstk ^= 0700700700;

	// Build in field out field and check minplus
	nuasb3_1 = nuasb3_2 = 0;
	BF128 w = valid_vect & b.g3x.gall;// use all 
	register int  Rfilt = b.smin.critbf;
	register uint32_t Rfstk = fstk;
	uint32_t  ua, nfree = 6 - b.smin.mincount;
	int iguan;
	b3_andout = fstk;

	while ((iguan = w.getFirst128()) >= 0) {
		w.Clear(iguan);
		GUAN & gu = tuguan.tgua3x3y[iguan];
		if (gu.ncol == 2)ua = b.guas.ua_pair[gu.i81];
		else if (gu.ncol == 3)ua = b.guas.ua_triplet[gu.i81];
		register uint32_t Ru = ua & BIT_SET_27;
		if (Ru & Rfilt) {
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(uasb3_1, nuasb3_1, Ru);
		}
		else {
			Ru &= Rfstk;
			if (!Ru) return 1; //not ok
			b3_andout &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(uasb3_2, nuasb3_2, Ru);
		}
	}

	if (nuasb3_2) {// check the limit 6
		if (!nfree) return 1;//not ok
		if ((!b3_andout) && nfree == 1) return 1;//not ok
	}
	{// add now band 3 uas ready to be loaded
		uint32_t * to = b.tua;
		for (uint32_t i = 0; i < b.nua; i++) {
			register uint32_t Ru = to[i];
			if (Ru & Rfilt)AddUA32(uasb3_1, nuasb3_1, Ru);
			else {
				Ru &= Rfstk;// loosing count
				if (!Ru)return 1;//not ok  
				b3_andout &= Ru;
				Ru |= _popcnt32(Ru) << 27;// restore count
				AddUA32(uasb3_2, nuasb3_2, Ru);
			}
		}
	}
	if (g17b.debug17) {
		valid_vect.Print("valid_vect");
		b.smin.Status( " b status CTXY_SetupB3MinplusCurib()");
		cout << nuasb3_1 << " " << nuasb3_2 << " if of final" << endl;
		void Debug_If_Of_b3();
	}
	if (nuasb3_2) {// check again the limit 6
		if (!nfree) return 1;//not ok
		b.smin.minplus++;
		if ((!b3_andout)) {
			if (nfree == 1) return 1;//not ok
			b.smin.minplus++;
		}
	}
	return 0;
}
void STD_B3::SetUpMincountxy(BF128 & validsockets) {
	BF128 ws = validsockets & g3x.g23;
	memset(&smin, 0, sizeof smin);
	int iguan;
	while ((iguan = ws.getFirst128()) >= 0) {
		ws.Clear(iguan);
		GUAN &gu = tuguan.tgua3x3y[iguan];
		int i81 = gu.i81, ncol = gu.ncol;
		if (ncol == 2) {
			int imini = guas.ua2_imini[i81], mask = (7 << (3 * imini));
			int  bit = 1 << imini, pat = guas.ua_pair[i81];
			if (smin.mini_bf2&bit) smin.mini_bf3 |= bit;
			if (smin.mini_bf1&bit) smin.mini_bf2 |= bit;
			smin.mini_bf1 |= bit;
			smin.critbf |= pat;
			smin.pairs27 |= mask ^ pat;
		}
		else smin.mini_triplet |= 1 << guas.ua3_imini[i81];
	}
	smin.SetMincount();
}
void BANDS_AB::CTXY_Gob3Curib() {
	p_cpt2g[10]++;
	//if (1) continue;//<<<<<<<<<<<<<test phase 1
	STD_B3 & b = genb12.bands3[cur_ib];
	//_________________________ this is a band3 to process 
	if (g17b.debug17) cout << "go band3 ib3=" << cur_ib << endl;
	moreuas_b3_small.Init();
	moreuas_b3.Init();
	memcpy(&genb12.grid0[54], b.band0, 4 * 27);
	smin = b.smin;
	//____________________ call the relevant band 3 process
	if (smin.minplus < 5) {// apply a matrix filter
		p_cpt2g[11]++;
		if((smin.minplus-smin.mincount)==2)
			if(smin.mincount==1)p_cpt2g[34]++;
			else if(smin.mincount == 2)p_cpt2g[33]++;
		MergeUasExpand();
		return;
	}
	{
		int nmiss = 6 - smin.mincount;
		if (g17b.debug17) {
			cout << " nmiss=" << nmiss << " nuas 1=" << nuasb3_1 << " nuas 2=" << nuasb3_2 << endl;
		}
		p_cpt2g[12 + nmiss] ++;
		G17B3HANDLER hh0;
		hh0.diagh = 0;// sbb.diag;

		//if ((!nmiss)  && p_cpt2g[12] < 100) hh0.diagh = 1;
		if (hh0.diagh) {
			smin.Status("b3 goband3 status before init");
		}

		hh0.Init();
		if (hh0.diagh) {
			cout << "stack count " << hh0.stack_count.u16[0] << hh0.stack_count.u16[1]
				<< hh0.stack_count.u16[2] << endl;
			hh0.smin.Status("b3 goband3 after init ");
		}
		if (!nmiss)hh0.Go_Critical();
		else 		if (nmiss == 1)	hh0.Go_miss1_b3();
		else if (nmiss == 2)	hh0.Go_miss2_b3();
		else hh0.Go_miss3_b3();
	}
}
void BANDS_AB::DebugAdd12(uint32_t itemp, uint64_t uab12, GINT64  w) {

	cerr << "ua < 12 to add clean" << endl;

	cout << endl << endl << Char2Xout(uab12) << " ua < 12 to add   clean" << endl;
	cout << "bug location band 2 id=" << genb12.nb12 << endl;
	cout << myband1.band << endl;
	cout << myband2.band << endl;
	cout << Char2Xout(w.u64) << " XY at call" << endl;
	cout << "ia=" << ia << " i3=" << i3 << " iy3=" << iy3 << " itemp=" << itemp << endl << endl;
	g17b.GodebugInit(0);
	DebugBuildUasAB(7);
	//myua = zh2b[0].ValidXY(tclues, nclues, 1);
}
			//___________________ direct expansion in band 3

void BANDS_AB::MergeUasExpand() {
	STD_B3 & b = genb12.bands3[cur_ib];
		if (!nuasb3_1) {// will be if no mini row ua
			ExpandBand3(uasb3_2, nuasb3_2);// direct on outfiel
			return;
		}
		if (!nuasb3_2) {// should be processed in another way???
			ExpandBand3(uasb3_1, nuasb3_1);// direct in field
			return;
		}
	//}
	uint32_t tua[500], nua=0,aig=0;
	//int ig4;
	for (uint32_t i1 = 0; i1 < nuasb3_1; i1++) {
		register uint32_t ua1 = uasb3_1[i1];
		if (aig ||_popcnt32(ua1) < 5) {
			tua[nua++] = ua1;
		}
		else {// insert gua4 and outfield in bloc
			aig = 1;// do it only once
			//while ((ig4 = w.getFirst128()) >= 0) {// first guas4
			//	w.Clear(ig4);
			//	tua[nua++] = b.g4pat[ig4];
			//}
			for (uint32_t i2 = 0; i2 < nuasb3_2; i2++)
				tua[nua++] = uasb3_2[i2];
			tua[nua++] = ua1;
		}
	}
	ExpandBand3(tua, nua);// now expand the merged table
}
void BANDS_AB::ExpandBand3(uint32_t *tua, uint32_t nua) {// find all 5 and 6 clues solutions
	if (g17b.debug17 >2) {
		cout << "expandb3 nua=" << nua << endl;
		for (uint32_t i = 0; i < nua; i++) cout << Char27out(tua[i]) << endl;
	}
	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3,
			stack[3];
	}spb3[7], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	// init the stack status
	for (int i = 0; i < 3; i++) {
		s3->stack[i] = stack_count.u16[i];
		if (s3->stack[i] == 6)s3->active_cells &= ~(07007007 << (3 * i));
	}
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tua[0] & s3->active_cells;
	int tcells[10];
	//____________ here start the search 6 clues
next:
	uint64_t ispot = s3 - spb3;
	if (g17b.debug17 > 2 && g17b.npuz == 15) {
		cout << Char27out(s3->all_previous_cells) << " ispot=" << ispot << " iuab3=" << s3->iuab3;
		cout << Char27out(tua[s3->iuab3]) << " ua";
		cout << Char27out(s3->possible_cells) << endl;
	}
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3+1; *sn3 = *s3; // (copy the stack count)
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		{// if the stack is limit update sn3 active
			int st = C_stack[cell];
			sn3->stack[st]++;
			if(sn3->stack[st]>5)
				sn3->active_cells &= ~(07007007 << (3 * st));
		}
		// nextspot:take the next available ua to loop		
		for (uint32_t i = s3->iuab3 + 1; i < nua; i++) {
			if (tua[i] & filter)continue;
			if (ispot >= 5) 	goto next;//passing the limit	
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & sn3->active_cells;
			s3 = sn3; // switch to next spot
			goto next;
		}

	}
	if (ispot < 5) {// use active as possible
		sn3->possible_cells =  sn3->active_cells;
		if (g17b.debug17 > 2 && g17b.npuz == 15) {
			cout << "\t\t" << Char27out(sn3->possible_cells) << " ispot=" << ispot << " continue with active cells" << endl;
		}
		sn3->iuab3 = nua;
		s3 = sn3; // switch to next spot
		goto next;
	}
	if (g17b.debug17 > 2 && g17b.npuz == 15) {
		cout << "\t\t" << Char27out(sn3->all_previous_cells)  << " possible solution to check" << endl;
	}
	// this is a possible 17 do final check
	bands_ab.FinalCheckB3(sn3->all_previous_cells);  
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
	//DebugExpand();
}

//_________ final called by all branches 
void BANDS_AB::FinalCheckB3(uint32_t bfb3) {
	int diagloc = 0;
	if (moreuas_b3_small.Check(bfb3))return;
	if (moreuas_b3.Check(bfb3))return;
	if (diagloc)cout << Char27out(bfb3) << " b3 to check" << endl;
	register uint32_t ir = zhou[1].CallMultipleB3(zhou[0], bfb3, 0);
	if (ir) {
		register uint32_t ua = zh_g2.cells_assigned.bf.u32[2];
		if (_popcnt32(ua) < 10)moreuas_b3_small.Add(ua);
		else moreuas_b3.Add(ua);
		if (diagloc)cout << Char27out(ua) << " b3 ua to add" << endl;
		return;
	}

	cout << Char32out(bfb3) << "one sol to print final check " << endl;
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < nclues; i++) {
		int cell = tclues[i];
		ws[cell] = genb12.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[54 + i] = genb12.grid0[54 + i] + '1';
	fout1 << ws << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
	g17b.a_17_found_here++;

}

//================ part 2  band 3 processing using guas2/3

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	if (diag) cout << Char27out(bf) << " call multipleb3" << endl;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{	
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = genb12.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) {
		if (diag) {
			cout << "not eight digits" << endl;
			ImageCandidats();
		}
		return 1;// can not be one solution
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (diag) ImageCandidats();
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (diag) {
		cout << "after update" << endl;
		ImageCandidats();
	}
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0,diag);

	return zh_g.nsol;  
}
int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4= R3 & FD[3][0]; 
		R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	BF128 R5 = R4 & FD[4][0]; R4 |= R3 & FD[4][0];
		R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R5 |= R4 & FD[5][0]; R4 |= R3 & FD[5][0];
		R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R5 |= R4 & FD[5][6]; R4 |= R3 & FD[6][0];
		R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R5 |= R4 & FD[7][0]; R4 |= R3 & FD[7][0];
		R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R5 |= R4 & FD[8][0]; R4 |= R3 & FD[8][0];
		R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		zh_g.pairs = R2 - R3;
		zh_g2.triplets = R3 - R4;
		zh_g2.quads = R4 - R5;
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diag) {
	if (diag) {
		char ws[82];
		cout << zh_g.pairs.String3X(ws)<< " pairs "<< endl;
		cout << zh_g2.triplets.String3X(ws) << " triplets " << endl;
		cout << zh_g2.quads.String3X(ws) << " quads " << endl;
	}
	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	// here the target is to have ua band 3 as small as possible
	
	if (w.bf.u32[2]) w.bf.u32[0] = w.bf.u32[1] = 0;// select band 3 in priority
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];
	if (diag) 
		cout << "guess17 index=" << index << " cell " << cellsFixedData[cell].pt << endl;

	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		if (diag)cout << "guess true" 
			<< digit + 1 << cellsFixedData[cell].pt << endl;
		mynext->Compute17Next(index + 1, diag);
		if (zh_g.go_back) return;

	}
	// if first step try first false
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			if (diag)cout << "guess false "
				<< idig + 1 << cellsFixedData[cell].pt << endl;
			mynext->Compute17Next(index + 1, diag);
			if (zh_g.go_back) return;
		}
	}
}
void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (diag>1) { cout << "index=" << index << endl; ImageCandidats(); }
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = genb12.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {// ignore true solution
				if (diag) ImageCandidats();
				zh_g.nsol++;
				zh_g.go_back = 1;// closed anyway
			}
		}
		return;
	}
	Guess17(index , diag);// continue the process
}

void G17B3HANDLER::Init( ) {
	BANDS_AB & bab = bands_ab;
	smin = bab.smin;
	uasb3of = bab.uasb3_2;	nuasb3of = bab.nuasb3_2;
	uasb3if = bab.uasb3_1;	nuasb3if = bab.nuasb3_1;
	andoutf = bab.b3_andout;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 =smin.critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = 6 - smin.mincount;
	nb3 = 6;
	stack_count = bab.stack_countf;
	for (int istack = 0, stp = 0111; istack < 3; istack++, stp <<= 1)
		if (stack_count.u16[istack] > 5) {// critical stack
			wactive0 &= ~(07007007 << (3 * istack));// clear outfield
			register int m2stack = stp & smin.mini_bf2, shrink = TblShrinkMask[m2stack];
			if (m2stack) {// common cell(s) to assign
				register int Mask = tbitsrows[shrink] << (3 * istack);
				//adjust count and known
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
				smin.mini_bf2 &= ~stp; // clear the 2pairs bit(s) in stack
				active_b3 &= (~Mask);// clear the  field bf
				smin.critbf &= (~Mask);
				smin.pairs27 &= (~Mask);
				smin.mincount -= _popcnt32(shrink);
			}
		}

}

uint32_t G17B3HANDLER::IsMultiple(int bf) {
	if (bf == rknown_b3) return 0;
	if (_popcnt32(bf) > 25) return 0;
	uint32_t ua = 0;
	rknown_b3 = bf;
	BANDS_AB & bab = bands_ab;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, 0);
	if (g17b.diag >= 2)	cout<<Char27out(ir) << "ir retour multiple ir="<<ir  << endl;
	return ua;
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (smin.mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[smin.mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~smin.pairs27);// and set the common cell as assigned
			}
		}
		smin.mini_bf2 = 0;
	}
}
//================= critical process
void G17B3HANDLER::CriticalAssignCell(int Ru){// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & smin.mini_bf3){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		smin.mini_bf3 ^= bit; // now only a pairto hit
		smin.mini_bf1 |= bit;
	}
	else{// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		smin.mini_bf1 &= ~bit;
		smin.mini_triplet &= ~bit;
	}
}


void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	p_cpt2g[27]++;
	if (g17b.diag >= 2|| diagh)	cout << Char27out(known_b3) << "entry critical nb3if= "<<
		nuasb3if << endl;
	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	CriticalLoop();
}
void G17B3HANDLER::CriticalLoop(){// after optional assignment  
	while (1) {// first shrink uas in field
		irloop = 0;
		uint32_t * tn = &uasb3if[nuasb3if], n = 0;
		register uint32_t Ra = active_b3,
			Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuasb3if; iua++) {
			register int Ru = uasb3if[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
			if (!Ru) return ;// dead branch
			if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
				CriticalAssignCell(Ru);
				Ra = active_b3; //can be  modified
				irloop = 1;// should loop for new singles
			}
			else tn[n++] = Ru;
		}
		uasb3if = tn;
		nuasb3if = n;
		if (!n) irloop = 0;// no need to loop again
		if (!irloop) break;
	}
	if (_popcnt32(known_b3) > 6) return; 
	if (!active_b3) {// must be here expected number of clues
		if (nuasb3if) return ; //can not be valid
		bands_ab.FinalCheckB3(known_b3);
		return ; // branch closed 
	}
	p_cpt2g[28]++;

	if (0) {
		if (diagh) {
			cout << Char27out(known_b3) << "end shrink  nb3if= "
				<< nuasb3if << " p_cpt2g[27] " << p_cpt2g[27] << endl;
			for (uint32_t i = 0; i < nuasb3if; i++)
				cout << Char27out(uasb3if[i]) << endl;
		}
		return;
	}
	int wua = uasb3if[0] & active_b3,cell;
	while (bitscanforward(cell, wua)) {
		register int bit = 1 << cell;
		wua ^= bit;// clear bit

		// clean the bit in active_b3, this is now a dead cell downstream
		active_b3 ^= bit;

		G17B3HANDLER hn = *this;
		hn.CriticalAssignCell(bit);
		hn.CriticalLoop();
	}
}

//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow() {
	if (g17b.debug17)cout << "entry Go_SubcriticalMiniRow() ndead=" << ndead << endl;
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++,  bit <<= 1, mask <<= 3) {
		stack = i % 3;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (g17b.debug17)cout << Char27out(M) << " mini row to process i=" << i << endl;
		if (bit & smin.mini_bf1) {// it was a gua2 pair assign both 
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf1 ^= bit;
			hn.SubMini( M, mask);
		}
		else if (bit & smin.mini_bf2)// it was 2 gua2 pair assign 2 out of 3 
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_bf2 ^= bit;
				hn.SubMini(M, mask);
			}
		else if (bit & smin.mini_bf3) {// it was 3 gua2 pair assign 3 out of 3 
			G17B3HANDLER hn = *this;
			hn.smin.mini_bf3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & smin.mini_triplet)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.smin.mini_triplet ^= bit;
				hn.SubMini(M, mask);
			}
		else {// second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void G17B3HANDLER::SubMini( int M, int mask) {
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added 
	active_b3 &= ~mask;
	active_sub ^= M;
	
	// now adjust the stack count
	stack_count.u16[stack]++;
	if (stack_count.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else {	// leave sub critical mode and enter the critical mode 
		Critical2pairs();// assign 2 pairs in minirow to common cell
		CriticalLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
	if (diagh) {
		cout << Char27out(known_b3) << " entry Go_Subcritical()" << endl;
		cout << Char27out(active_sub) << " active_sub " << endl;
	}
	for (int ist = 0; ist < 3; ist++) {// check stacks 
		if (stack_count.u16[ist] > 5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	if (0) {// likely not very productive now
		int cct = _popcnt32(smin.critbf) - smin.mincount;
		if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
		if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	}
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}

//======================================================================= not critical sequence
void G17B3HANDLER::ShrinkUasOfB3() {
	if (known_b3) {// shrink the out field table
		uint32_t * tn = &uasb3of[nuasb3of], n = 0;
		register uint32_t Ra = wactive0,
			Rfilt = known_b3;
		andoutf = BIT_SET_27;
		for (uint32_t iua = 0; iua < nuasb3of; iua++) {
			register int Ru = uasb3of[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		uasb3of = tn;
		nuasb3of = n;
	}
}

void G17B3HANDLER::Go_miss1_b3() {// not called if more than 1 needed
	if (diagh) cout << Char27out(known_b3) << " entry Go_miss1_b3()" << endl;
	wua = wactive0;
	ShrinkUasOfB3();// if known from up stream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				if (diagh) cout << Char27out(uabr) << " nmiss1 first ua of found" << endl;
				andoutf = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical();
			}
		}
	}
	if (nuasb3of) wua &= andoutf;
	if (diagh) cout << Char27out(wua) << " nmiss1 b3 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	
			if (hn.AddCell_Of(res, bit))hn.Go_Critical();
		}
	}
}
void G17B3HANDLER::Go_miss2_b3() {
	ShrinkUasOfB3();// if known from upstream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				if (diagh) cout << Char27out(uabr) << " nmiss2 b3 first ua of found" << endl;
				uasb3of[0] = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical(); 
			}
		}
	}
	wua = wactive0;
	if (nuasb3of)wua &= uasb3of[0];// use first ua  
	if (diagh) cout << Char27out(wua) << " nmiss2 b3 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;
			if(hn.AddCell_Of(res,bit))hn.Go_miss1_b3();
		}
	}
}
void G17B3HANDLER::Go_miss3_b3() {// always first direct entry
	if (diagh) cout << Char27out(wua) << " nmiss3 b3 wua" << endl;
	ShrinkUasOfB3();// if known from upstream
	if (smin.mincount) {
		if (!nuasb3of) {// subcritical in hn if solved
			int uabr = IsMultiple(known_b3 | active_b3);
			if (uabr) {// on ua outfield seen
				if (diagh) cout << Char27out(uabr) << " nmiss3b3 first ua of found" << endl;
				uasb3of[0] = uabr;
				nuasb3of = 1;
			}
			else {// confirmed subcritical possible
				G17B3HANDLER hn = *this;
				hn.Go_Subcritical();
			}
		}
	}
	 wua = wactive0;	
	 if (nuasb3of)wua &= uasb3of[0];// use first ua  
	 if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;
			if (hn.AddCell_Of(res, bit))hn.Go_miss2_b3();
		}
	}
}


//=============================== debugging sequences

int G17B::DebugK17M10() {
	cout << "b1 i416=" << myband1.i416 << " b2 i416=" << myband2.i416
		<< " n5b1=" << myband1.nmybv5 << " n5b2=" << myband2.nmybv5 << endl;
	GodebugInit(0);
	if (GodebugCheckUas("check uas")) return 1;
	if (GodebugFindKnown17()) {
		cout << "failed to find the expected 17 in expanded bands" << endl;
		return 1;
	}
	else cout << iband1 << " " << iband2 << " seen solution" << endl;
	if (0) {
		if (debug17 > 2) {
			cout << "index band 1" << endl;
			myband1.DebugExp();
			cout << "index band 2" << endl;
			myband2.DebugExp();

		}
	}
	return 0;
}
void G17B::GodebugInit(int mode) {// called from DebugK17M10
	cout << "band1 n5=" << myband1.nmybv5 << " n6=" << myband1.nmybv6 << endl;
	cout << "band2 n5=" << myband2.nmybv5 << " n6=" << myband2.nmybv6 << endl;
	cout << "n bands3      \t" << genb12.nband3  << "\tua bands1+2   \t" << genuasb12.nua;
	cout << "\nnguas socket2  \t" << genb12.ntua2 << endl;
	cout << "nguas socket3  \t" << genb12.ntua3<<endl;
	cout << "guas number sockets  2;3;all\t"
		<<tuguan.ng2<<"\t" << tuguan.ng3 << "\t" << tuguan.nguan << endl;

	if (debug17 > 1) {
		STD_B3 &b3 = genb12.bands3[0];
		cout << "band3 permanent data" << endl;
		b3.guas.isguasocket2.Print3("guasocket2");
		b3.guas.isguasocket3.Print3("guasocket3");
		b3.guas.isguasocket2_46.Print3("guasocket2_46");
	}


	if (mode & 1)tuguan.Debug1();// guas list and killers

	if (mode & 2) {
		cout << "table uas" << endl;
		uint64_t *t = genuasb12.tua;
		uint32_t n = genuasb12.nua;
		for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << endl;

	}
	if (mode & 4) {// status of guas2 3 first band
		cout << "debug status of guas2 3 first band" << endl;
		STD_B3::GUAs & bguas = genb12.bands3[0].guas;
		BF128 ws = bguas.isguasocket2 | bguas.isguasocket3;
		for (uint32_t i = 0; i < 128; i++) {
			GUAN wg = tuguan.tguan[i];
			if (wg.ncol > 3) continue;
			int i81 = wg.i81;
			if (wg.ncol == 2 && ws.On_c(i81)) {// here 3X 27
				if (bguas.isguasocket2.On_c(i81))
					cout << Char27out(bguas.ua_pair[i81]) << " sock2";
				else cout << Char27out(bguas.ua_triplet[i81]) << " sock3";

				cout <<"i = " << i	<< " i81=" << i81 << " "
					<< Char2Xout(wg.killer) << "killer" << endl;
				uint64_t *tua = wg.tua;
				uint32_t nua = wg.nua;
				for (uint32_t iu = 0; iu < nua; iu++)
					cout << Char2Xout(tua[iu]) << endl;
			}

		}
		tuguan.Debug1();

	}
	

}
int G17B::GodebugFindKnown17() {// locate the known 17 
	uint32_t *t1, nt1, *t2, nt2;
	if (_popcnt32(p17diag.bf.u32[0]) == 5) {
		cout << "this is a 5 6 6 known 17" << endl;
		t1 = myband1.mybv5; nt1= myband1.nmybv5;
		t2 = myband2.mybv6; nt2 = myband2.nmybv6;
	}
	else {
		cout << "this is a 6 5 6 known 17" << endl;
		t1 = myband1.mybv6; nt1 = myband1.nmybv6;
		t2 = myband2.mybv5; nt2 = myband2.nmybv5;
	}
	register uint32_t R1 = p17diag.bf.u32[0];
	for (iband1 = 0; iband1 < nt1; iband1++) {
		if (R1 == t1[iband1])goto check2;
	}
	cout << "fault band1" << endl;
	return 1;// fault band1
check2:;
	R1 = p17diag.bf.u32[1];
	for (iband2 = 0; iband2 < nt2; iband2++) {
		if (R1 == t2[iband2]) {
			cout << "found iband 1=" << iband1 << " iband 2=" << iband2 << endl;
			return 0;
		}
	}	
	cout << "fault band2" << endl;
	return 1; //fault band 2
}
int G17B::GodebugCheckUas(const char * lib) {
	uint32_t nua = genuasb12.nua;
	uint64_t * tua= genuasb12.tua;
	for (uint32_t i = 0; i < nua; i++) {
		if (tua[i] & g17b.p17diag.bf.u64[0]) continue;
		cout << lib << "check ua failed" << endl;
		cout << Char2Xout(tua[i]) << " not hit by the known 17" << endl;
		return 1;
	}
	return 0;
}

void BANDS_AB::DebugBuildUasAB(int mode) {
	cout << Char27out(bf3) << " debug builduasAB nindex=" << ntiBc
		<< " nuas=" << ntuaB << endl;
	//b384_bf3.Print(" 384_bf3");
	if(mode &1)
	for (uint32_t i = 0; i < ntuaB; i++)
		cout << Char2Xout(tuaB[i])<<" i="<<i << endl;
	if (mode & 2)
		for (uint32_t i = 0; i < ntiBc; i++) {
		cout << Char27out(tiBc_pat[i])  << " i=" << i
			<< "\t ideb=" << tiBc_ideb[i];
		//if (tiBc_cc[i] < 4)cout << "\timini=" << tiBc_imini[i];
		cout << "\t" << Char27out(tiBc_kill[i]) << " kill" << endl;
	}
	//cout << Char64out(b64B_pat23) << " pairs triplets status" << endl;
	//v512_mult.Print(" v512_mult status ");

}

void BANDS_AB::Debug_If_Of_b3() {
	cout << "nif=" << nuasb3_1 << endl;
	for (uint32_t i = 0; i < nuasb3_1; i++)
		cout << Char27out(uasb3_1[i] )<< endl;
	cout << "nof=" << nuasb3_2 << endl;
	for (uint32_t i = 0; i < nuasb3_2; i++)
		cout << Char27out(uasb3_2[i]) << endl;
}
