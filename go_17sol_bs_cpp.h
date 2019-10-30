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
		cout << myband2.band << "go band2 id=" << myband2.i416 << " nb12=" << genb12.nb12 << endl;
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
	GoM10_guas_four_columns();
	//tuguan.BuildCellKillersVector();
}
void G17B::GoM10_guas_four_columns() {// add guas 4 columns 2 digits 2 rows 2 boxes
	int nb3 = genb12.nband3, ntguan4 = 0,limit=256-tuguan.nguan;
	for (int ib1 = 0; ib1 < 2; ib1++) {
		for (int i = 0; i < 27; i++) {
			int i8a = 27 * ib1 + i;// first i81
			GEN_BANDES_12::SGUA2  wa = genb12.tsgua2[i8a];
			uint32_t cbf1 = (Zhoucol << wa.col1) | (Zhoucol << wa.col2),
				cols1 = (1 << wa.col1) | (1 << wa.col2);
			for (int ib2 = ib1 + 1; ib2 < 3; ib2++) {
				for (int j = 0; j < 27; j++) {
					int i8b = 27 * ib2 + j;// second  i81
					GEN_BANDES_12::SGUA2  wb = genb12.tsgua2[i8b];
					if (wa.digs != wb.digs) continue; // same digits needed
					uint32_t cbf2 = cbf1 | (Zhoucol << wb.col1) | (Zhoucol << wb.col2);
					TEMPGUAN4 & tg = tempguan4[ntguan4++];
					tg.Init(wa.digs, cols1 | (1 << wb.col1) | (1 << wb.col2));
					for (int ib3 = 0; ib3 < nb3; ib3++) {
						STD_B416 & b3 = genb12.bands3[ib3];
						uint32_t patd1d2 = b3.fd_sols[0][wa.dig1] | b3.fd_sols[0][wa.dig2];
						patd1d2 &= cbf2; // limit to the four columns 
						int nr = 0;
						for (int ish = 0; ish <= 18; ish += 9)
							if ((patd1d2 >> ish) & 0x1ff) nr++;
						if (nr != 2) continue; // must be 2 rows 
						int nmini = 0;
						for (int imini = 0, mask = 7; imini < 9; imini++, mask <<= 3)
							if (patd1d2&mask)nmini++;
						if (nmini != 4)continue;// forget 2 mini
						tg.AddBand(ib3, patd1d2);
					}
					if (ntguan4 >= limit) goto exit_collect_sockets;
				}
			}
		}
	}
exit_collect_sockets:;

	if (0) {
		cout << "_collect_sockets n= " << ntguan4 <<" limit="<<limit<< endl;
		//return;
	}
	// Find GUAs bands 1+2 for the socket 4
	for (int itg = 0; itg < ntguan4; itg++) {
		uint32_t patcols = tempguan4[itg].colbf,
			d1d2bf= tempguan4[itg].digsbf;
		// build revised gangster
		int gangcols[9];
		memcpy(gangcols, genb12.gangb12, sizeof gangcols);
		for (int icol = 0, bit = 1; icol < 9; icol++, bit <<= 1)
			if (patcols&bit)gangcols[icol] ^= d1d2bf;
		// find guas of the socket
		zh2b_g.InitGangster(genb12.gangb12, gangcols);
		zh2b5_g.sizef5 = GUALIMSIZE;
		zh2b5_g.modevalid = 0;
		genb12.ptua2 = tuguan.pguabuf;
		uint64_t *tua = tuguan.pguabuf;
		genb12.nua2 = 0;
		//================== GUA collector 2 bands
		genb12.GuaCollect(d1d2bf);
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if ((fl & d1d2bf) == d1d2bf)	genb12.GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = floors_4d[i];
			if ((fl & d1d2bf) == d1d2bf) genb12.GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 5digits
			int fl = 0x1ff ^ floors_4d[i];
			if ((fl & d1d2bf) == d1d2bf) genb12.GuaCollect(fl);
		}
		if (genb12.nua2) {
			if (genb12.nua2 > 40)genb12.nua2 = 40;
			if (0) {
				cout << "ua collector got " << genb12.nua2 << " uas " << endl;
				uint64_t killer = BIT_SET_2X;
				for (uint32_t i = 0; i < genb12.nua2; i++) 	killer &= tua[i];
				cout << Char2Xout(killer) << " killer" << endl;
			}
			tuguan.AddGuan(tua, genb12.nua2, patcols, d1d2bf, itg);
		}
	}
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
		memcpy(zh2b_g.puz0, myband2.band0, sizeof myband2.band0);
		memcpy(&zh2b_g.puz0[27], myband1.band0, sizeof myband1.band0);
		zh2b_g.GetBands(myband2.gangster, myband1.gangster);// set sol/pm
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

	for (iy3 = 0; iy3 < nyi3; iy3++) {
		// store the y 6 clues in a ficed place 
		uint32_t yindd = mybb->myi3[iy3].ideb,
			yindf = mybb->myi3[iy3 + 1].ideb;
		ny3 = yindf - yindd;
		//cout << "ny3=" << ny3 << endl;
		memcpy(tybf, &mybb->mybv6[yindd], 4 * ny3);
		wiy3 = mybb->myi3[iy3];
		for (i3 = 0; i3 < ni3; i3++) {
			indd = myba->myi3[i3].ideb5;
			indf = myba->myi3[i3 + 1].ideb5;
			if (indd == indf)continue;// no X 5 clues here
			nx3 = indf - indd;
			memcpy(txbf, &myba->mybv5[indd], 4 * nx3);
			wi3 = myba->myi3[i3];
			Go3X3Y();
		}
	}
}
void BANDS_AB::Go3X3Y() {
	moreuas_AB.Init();
	moreuas_AB_small.Init();
	moreuas_AB_big.Init();
	indd = myba->myi3[i3].ideb5;
	indf = myba->myi3[i3 + 1].ideb5;
	if (indd == indf)return;// no X 5 clues here
	wi3 = myba->myi3[i3];
	p_cpt2g[2]++;
	p_cpt2g[3] += nx3 * ny3;
	bf3 = wi3.cellsbf;
	ybf3 = wiy3.cellsbf;
	if (ia) b1b2_3x3y_bf = ((uint64_t)bf3 << 32) | ybf3;
	else b1b2_3x3y_bf = ((uint64_t)ybf3 << 32) | bf3;
	  
	if (g17b.debug17) {
		if(b1b2_3x3y_bf & (~g17b.p17diag.bf.u64[0]))return;
		cout <<Char2Xout(b1b2_3x3y_bf)<<" this is a valid 3X3Y"<<endl;
	}

	tuguan.Build_Guas3X3Y();
	if (p_cpt2g[35] < tuguan.nuar)p_cpt2g[35] = tuguan.nuar;
	//tuguan.Debug3X3Y();
	//if (1) return;
	//if (p_cpt2g[2] < 429) 	return;	
	//if (p_cpt2g[2] > 434) 	return;
	//  if (p_cpt2g[2] == 429) tuguan.Debug3X3Y();

	{	// Sort  on B increasing number of cells
		register uint64_t F = (ia) ? bf3 : ybf3;
		F <<= 32;// band 2 high bits
		F |= (ia) ? ybf3 : bf3; // filter for uas
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
	p_cpt2g[32] += ntuaB;
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
	if (ntiBc <= 128)Go3X3Y128();
	else Go3X3Y256();
	if (nxy_filt1)CleanTempXY();// final clean 3X3Y

}

void TU_GUAN::Build_Guas3X3Y() {// guas filter => already killed/forced plus sub table
	pguabufr = guabufr;
	ngua3x3y =nuar23= 0;
	{
		register uint64_t F = bands_ab.b1b2_3x3y_bf;
		vmult.SetAll_0();
		// go through the first table to extract still valid
		for (uint32_t ig = 0; ig < nguan; ig++) {
			GUAN & g = tguan[ig];
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
				if (gn.nua > 25) break;
			}
			if (!gn.nua)continue;
			gn.killer = killer;
			if (gn.nua > 1)vmult.setBit(ngua3x3y);
			if (gn.ncol < 4) 	nuar23++;
			tgua3x3y[ngua3x3y++] = gn;
			pguabufr += gn.nua;
		if(ngua3x3y>127) break; // limit 128 here
		}
	}
	nuar =(uint32_t)( pguabufr - guabufr);

	if (nuar > 128 * NGUARBLOCS) nuar = 128 * NGUARBLOCS;
	{ //Setup tuaindex for quick check of a given gua
		int indd = 0, indf;
		for (uint32_t i = 0; i < ngua3x3y; i++) {
			GUAN3X3Y & ggi = tuaindex[i];
			indf = indd + tgua3x3y[i].nua;
			ggi.i = indd / 128;
			ggi.indd = indd % 128;
			ggi.indf = indf % 128;
			indd = indf;// for the next tgua3x3y
		}
	}
	{// setup "cells to guan" vectors
		v3x3y = maskLSB[ngua3x3y];
		memset(vcells3x3y, 255, sizeof vcells3x3y);// all bits to 1
		uint32_t cc;
		for (uint32_t i = 0; i < ngua3x3y; i++) {// set uas
			register uint64_t Rw = tgua3x3y[i].killer;
			Rw &= BIT_SET_2X;
			while (bitscanforward64(cc, Rw)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= bit2;// clear bit
				vcells3x3y[From_128_To_81[cc]].Clear(i);
			}
		}
	}
	{// setup "cells to guan uas" vectors
		vua3x3y.Init((uint32_t)nuar);
		memset(vuac3x3y, 255, sizeof vuac3x3y);
		uint32_t cc;
		for (uint32_t iua = 0; iua < nuar; iua++) {
			register int bloc = iua >> 7, ir = iua & 127;
			register uint64_t Rw = guabufr[iua], biti = (uint64_t)1 << ir;
			Rw &= BIT_SET_2X;
			while (bitscanforward64(cc, Rw)) {// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= bit2;// clear bit
				vuac3x3y[From_128_To_81[cc]].v[bloc].clearBit(ir);
			}
		}
	}
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		int blocib3 = ib3 >> 7, irib3 = ib3 & 127;
		STD_B3 & b = genb12.bands3[ib3];
		STD_B3::GUAs &g = b.guas;
		uint32_t * td = b.fd_sols[0];// per digit status
		memset(&b.g3x3y, 0, sizeof b.g3x3y);
		for (uint32_t ig = 0; ig < ngua3x3y; ig++) {
			GUAN gu = tgua3x3y[ig];
			int i81=gu.i81,stack=i81/27;
			if ((gu.ncol == 2 && g.isguasocket2.On_c(i81)) ||
				(gu.ncol == 3 && g.isguasocket3.On_c(i81))) {
				b.g3x3y.g23.Set(ig);
				b.g3x3y.g23s[stack].Set(ig);
			}
		
			else if (gu.ncol == 2 && g.isguasocket2_46.On_c(i81)) {
				b.g3x3y.g46.Set(ig);
			}
			else if(gu.ncol == 4 ) {// expected  4cols
				TEMPGUAN4 &tp = tempguan4[gu.i81];
				if (tp.b3bf.On(blocib3, irib3)) {//ok for this band
					b.g3x3y.g46.Set(ig);
				}
			}
		}
	}
}
void TU_GUAN::Debug3X3Y(){
	cout << "test build nguan=" << nguan<< " n3x3y=" << ngua3x3y
		<<" nuar="<<nuar<< endl;
	cout << Char2Xout(bands_ab.b1b2_3x3y_bf) << " 3x3ybf" << endl;
	//Debug1();
	v3x3y.Print(" v3x3y");
	vmult.Print(" vmult");
	for(uint32_t i=0;i< ngua3x3y;i++)tgua3x3y[i].Debug1Guan(i);
	for (uint32_t i = 0; i < nuar; i++)
		cout << Char2Xout(guabufr[i]) << " uai=" << i << endl;
	//vua3x3y.Print(" vua3x3y vector");
	//vuac3x3y[10].Print(" vuav3x3y[10] vector");


	//for (int i = 0; i < 54; i++)
		//cout << Char64out(vcells3x3y[i].bf.u64[0]) << " vc i=" << i << endl;
	cout << "\n\nfirst band status" << endl;
	STD_B3 & b = genb12.bands3[0];
	b.g3x3y.g23.Print("b.g3x3y.g23");
	b.g3x3y.g23s[0].Print("b.g3x3y.g23s[0]");
	b.g3x3y.g46.Print3("b.g3x3y.g46");
}

void TU_GUAN::GUA23Collect() {// only guas 23 here
	xy23sockets.SetAll_0();
	//cout << "entry gua23collect" << endl;
	vxy_raw = v3x3y;
	vuaxy = vua3x3y;
	uint64_t dbf = bands_ab.b1b2_xy_bf^ bands_ab.b1b2_3x3y_bf;
	uint32_t clue, cell;
	//cout << Char2Xout(dbf) << "dbf" << endl;
	while (dbf) {
		bitscanforward64(clue, dbf);
		dbf ^= (uint64_t)1 << clue;
		cell = From_128_To_81[clue];
		vxy_raw &= vcells3x3y[cell];
		vuaxy.And(vuac3x3y[cell]);
	}// now vxy_raw vuaxy free of erased through "kill"
	//vxy_raw.Print("vxy_raw");
	BF128 wvxy = vxy_raw & maskLSB[nuar23]; // here only potential mini row sockets
	int iguan;
	while ((iguan = wvxy.getFirst128()) >= 0) {
		wvxy.Clear(iguan);
		GUAN3X3Y & gi = tuaindex[iguan];// always below 128 uas 
		if (gi.indf>gi.indd) {// same bloc
			BF128 w= maskLSB[gi.indf];
			w-= maskLSB[gi.indd];
			w &= vuaxy.v[gi.i];
			if (w.isNotEmpty()) goto nextiguan;
		}
		else {// this is a 2 blocs vector
			BF128 w; w.SetAll_1();
			w -= maskLSB[gi.indd];
			w &= vuaxy.v[gi.i];
			if (w.isNotEmpty()) goto nextiguan;
			w= maskLSB[gi.indf];
			w &= vuaxy.v[gi.i+1];
			if (w.isNotEmpty()) goto nextiguan;
		}
		vxy_raw.Clear(iguan);
		continue;
	nextiguan:// set xy2sockets,xy3sockets;
		GUAN & g= tgua3x3y[iguan];
		if (g.ncol <4 )xy23sockets.Set_c(g.i81);
	}
	//vxy_raw.Print(" vxy  clean detail limited to 23 ");

}
void TU_GUAN::GUAEndCollect() {// guas vector after 6 clues B
	// first check reduced vector to get potential active reduced guan
	//cout << "entry InitC_Guas()"  << endl;
	vxy = vxy_raw;// restart from vector valid for mini rows

	// for each valid not mini row, check if a ua not hit is there

	BF128 wvxy = vxy_raw;
	wvxy-=maskLSB[nuar23]; // already checked
	int iguan;
	while ((iguan = wvxy.getFirst128()) >= 0) {
		wvxy.Clear(iguan);
		GUAN3X3Y & gi = tuaindex[iguan];// always below 128 uas 
		if (gi.indf > gi.indd) {// same bloc
			BF128 w = maskLSB[gi.indf];
			w -= maskLSB[gi.indd];
			w &= vuaxy.v[gi.i];
			if (w.isNotEmpty()) continue;
		}
		else {// this is a 2 blocs vector
			BF128 w; w.SetAll_1();
			w -= maskLSB[gi.indd];
			w &= vuaxy.v[gi.i];
			if (w.isNotEmpty()) continue;;
			w = maskLSB[gi.indf];
			w &= vuaxy.v[gi.i + 1];
			if (w.isNotEmpty()) continue;
		}
		vxy.Clear(iguan);
	}
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

void BANDS_AB::CleanTempXY() {//GINT64 tempXY[15000];
	//if (g17b.debug17>1)//|| g17b.diagbug == 2)
		//cout << "entry CleanTempXY() nxy_filt1=" << nxy_filt1
		//	<<" p_cpt2g[2]="<< p_cpt2g[2] << endl;
		//if (1) return;
	//if(i3==57 && ia)cout << "entry CleanTempXY() ntempx=" << nxy_filt1 << endl;
	p_cpt2g[5]+= nxy_filt1;
	//if (1) { nxy_filt1 = 0; return; }
	for (uint32_t itemp = 0; itemp < nxy_filt1; itemp++) {
		GINT64  w = tempXY[itemp];
		//if (g17b.diagbug == 2)cout << Char2Xout(w.u64) << "XY" << endl;
		/*
		if (g17b.diagbug == 2) {
			if (w.u32[0] != g17b.k17x(ia)) continue;
			if (w.u32[1] != g17b.k17x(ib)) continue;
			cout << Char2Xout(w.u64) << "XY we have the right band B" << endl;
		}
		if (g17b.diagbug == 2) {
			if (w.u32[0] == g17b.k17x(ia) && w.u32[1] == g17b.k17x(ib))
				cout << Char27out(w.u32[1]) << "we have the right band B" << endl;
		}
		*/
		bfA = w.u32[0]; bfB = w.u32[1];
		// check the "more" table if available
		if (moreuas_AB_small.Check(w.u64))continue;
		if (moreuas_AB.Check(w.u64)) continue;
		if (moreuas_AB_big.Check(w.u64)) continue;

		// rebuild table of clues and build stack count for A,B
		if (ia) b1b2_xy_bf = ((uint64_t)bfA << 32) | bfB;
		else b1b2_xy_bf = ((uint64_t)bfB << 32) | bfA;
		if (g17b.debug17) {//skip if not ok
			if (b1b2_xy_bf != g17b.p17diag.bf.u64[0]) continue;
			cout << Char2Xout(b1b2_xy_bf) << "we have the right XY" << endl;
		}
		nclues = 0;
		uint32_t cell, wbf = bfA;
		stack_count.u64 = 0;
		while (wbf) {
			bitscanforward(cell, wbf);
			wbf ^= 1 << cell;
			stack_count.u16[C_stack[cell]]++;
			tclues[nclues++] = cell;
		}
		wbf = bfB;
		while (wbf) {
			bitscanforward(cell, wbf);
			wbf ^= 1 << cell;
			stack_count.u16[C_stack[cell]]++;
			tclues[nclues++] = cell + 27;
		}
		if (stack_count.u16[0] > 6 || stack_count.u16[1] > 6 ||
			stack_count.u16[2] > 6) continue;

		//______ quick check of the critical stack
		tuguan.GUA23Collect();
		if (FirstCheckGuaStack()) continue; //stack>6 in quick mode
		if (FirstCheckGuaAll()) continue; //mincount>6 without minplus
		//if (1) return;
		p_cpt2g[6]++;
		if (g17b.diagbug == 2) {
			if (w.u32[0] == g17b.k17x(ia) && w.u32[1] == g17b.k17x(ib))
				cout << Char27out(w.u32[1]) << " right band B try check valid" << endl;
			//else continue;//<< 
		}
		//_____ check for unique solution bands A;B
		register uint64_t myua = zh2b[0].ValidXY(tclues, nclues, 0);
		if (myua) {//___ if not do the best with the UA
			//if (g17b.debug17 ) 		cout << Char2Xout(myua) << "uaret band b multiple" << endl;
			register uint64_t uab = myua >> 32, uaa = myua & BIT_SET_27, uab12 = myua;
			if (ia)uab12 = uab | (uaa << 32);
			uint64_t cc64 = _popcnt64(uab12);
			if (cc64 < 12) {// this should never be check for a bug
				DebugAdd12(itemp,uab12,w);
				g17b.aigstop = 1;
				return;
			}
			if (cc64 < 18)			moreuas_AB_small.Add(myua);
			else if(cc64<21)			moreuas_AB.Add(myua);
			else moreuas_AB_big.Add(myua);
			if (cc64 < 18) {//  update the UA table (for next 3 clues band A)
				//cout << Char2Xout(uab12) << "to add" << endl;
				genuasb12.AddUACheck(uab12 | ((uint64_t)cc64 << 59));
				p_cpt2g[31]++;
			}
			continue;
		}// end not a unique solution
		if (g17b.diagbug == 2) cout  << " right band B  valid working now on band 3" << endl;
		memcpy(tcluesb12, tclues, 11 * 4);
		if (ia) for (int i = 0; i < 11; i++)
			if (i > 4)tcluesb12[i] -= 27; else tcluesb12[i] += 27;
		if (zhou[0].PartialInitSearch17(tcluesb12, nclues))// A;B order ok to check
			return;// would be  bug

		//________________________ after a valid bands 1+2 seen
		p_cpt2g[7]++;
		tuguan.GUAEndCollect();;// get final gangster valid vector 
		if (EndCheckCurBand3()) continue;
		if (0) {
			cout << "minplus=" << genb12.bands3[cur_ib].smin.minplus << endl;
			Debug_If_Of_b3();
		}
		GoBand3();
		while (++cur_ib< genb12.nband3) {
			if (FirstCheckGuaAll()) continue; //quickcheck till ok
			if (EndCheckCurBand3()) continue; // end check cur ib
			GoBand3();// passing minplus filter 
		}
	}
	nxy_filt1 = 0;
}

void BANDS_AB::GoBand3() {
	if (g17b.debug17) cout << "go band3 ib3=" << cur_ib << endl;
	STD_B3 & b = genb12.bands3[cur_ib];
	moreuas_b3_small.Init();
	moreuas_b3.Init();
	memcpy(&genb12.grid0[54],b.band0, 4 * 27);
	smin = b.smin;
	//____________________ call the relevant band 3 process
	if (smin.minplus < 5|| smin.mincount < 3	) {// apply a matrix filter
		p_cpt2g[10]++;
		MergeUasExpand();
		return;
	}
	p_cpt2g[11]++;




	int nmiss = 6 - smin.mincount;
	if (g17b.debug17) {
		cout << " nmiss=" << nmiss << " nuas 1=" << nuasb3_1 << " nuas 2=" << nuasb3_2 << endl;
		if (g17b.debug17) {
			cout << "table uab3_1" << endl;
			for (uint32_t i = 0; i < nuasb3_1; i++)
				cout << Char27out(uasb3_1[i]) << " i=" << i << endl;
			cout << "table uab3_2" << endl;
			for (uint32_t i = 0; i < nuasb3_2; i++)
				cout << Char27out(uasb3_2[i]) << " i=" << i << endl;
		}
	}
	p_cpt2g[12 + nmiss] ++;
	G17B3HANDLER hh0; hh0.Init();
	hh0.diagh = 0;// sbb.diag;
	if (hh0.diagh) {
		cout << "status after init" << endl;
		hh0.smin.Status("b3 goband3 ");
	}
	if (!nmiss)hh0.Go_Critical();
	else if (nmiss == 1)	hh0.Go_miss1_b3();
	else if (nmiss == 2)	hh0.Go_miss2_b3();
	else hh0.Go_miss3_b3();
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

int BANDS_AB::FirstCheckGuaStack() {// quick check of critical situations
	// consider critical stack (usually one)
	int maxcount = 0, stack;
	for (int i = 0; i < 3; i++) {
		int st_count = stack_count.u16[i];
		if (st_count > maxcount) { maxcount = st_count; stack = i; }
	}
	if (maxcount < 5)return 0;// 443 small chances to exit false
	if (maxcount >6)return 1; // no band ok
	cur_ib = 0;
	for (; cur_ib < genb12.nband3; cur_ib++) {
		STD_B3 & b = genb12.bands3[cur_ib];
		BF128 ws = tuguan.vxy_raw & b.g3x3y.g23s[stack];
		int count = ws.Count();
		if (!count) return 0;// one non critical band found
		if(maxcount==6 && count)return 0;
		// need 2 minimum clues in the stack to exclude the band
		if (count >= 4)continue; // (2 minirows needed or 3 pairs in a minirow)
		p_cpt2g[8]++;
		if (0) {
			cout << "entry first check need 2 cluesstack=" << stack << " maxcount=" << maxcount
				<< " ws count=" << count << endl;
			tuguan.vxy_raw.Print(" vxy_raw");
			b.g3x3y.g23s[stack].Print(" b.g3x3y.g23s[stack]");
			ws.Print(" ws");
		}
		MINCOUNT sm;	 memset(&sm, 0, sizeof sm);
		int iguan,bf1_old=0;
		while ((iguan = ws.getFirst128()) >= 0) {
			ws.Clear(iguan);
			GUAN &gu = tuguan.tgua3x3y[iguan];
			int i81 = gu.i81,ncol=gu.ncol;
			if (ncol == 2) {
				int  bit = 1 << b.guas.ua2_imini[i81];
				if (sm.mini_bf2&bit) sm.mini_bf3 |= bit;
				if (sm.mini_bf1&bit) sm.mini_bf2 |= bit;
				sm.mini_bf1 |= bit;
			}
			else sm.mini_bf1 |= 1 << b.guas.ua3_imini[i81];// can use here bf1 pairs closed
			if (bf1_old) {	if (sm.mini_bf1 != bf1_old) goto nextb3;}
			else bf1_old = sm.mini_bf1;// stop as soon as you get 2 mini rows
		}
		sm.mincount = _popcnt32(sm.mini_bf1) + _popcnt32(sm.mini_bf3) ;
		if (0)sm.Status("exit first check ");
		if (sm.mincount < 2) return 0; // non critical (quick) band found
		nextb3:;
	}
	return 1;// all bands are in critical mode skip this xy

}
void STD_B3::SetUpMincountxy() {
	p_cpt2g[9]++;
	BF128 ws = tuguan.vxy_raw & g3x3y.g23;
	if (p_cpt2g[9]<10) {
		cout << "entry second quick check ="<< endl;
		tuguan.vxy_raw.Print(" vxy_raw");
		g3x3y.g23.Print(" g3x3y.g23");
		ws.Print(" ws");
	}
	memset(&smin, 0, sizeof smin);
	int iguan;
	while ((iguan = ws.getFirst128()) >= 0) {
		ws.Clear(iguan);
		GUAN &gu = tuguan.tgua3x3y[iguan];
		int i81 = gu.i81, ncol = gu.ncol;
		if (ncol == 2) {
			int imini = guas.ua2_imini[i81], mask = (7 << (3 * imini));
			int  bit = 1 << imini,pat= guas.ua_pair[i81];
			if (smin.mini_bf2&bit) smin.mini_bf3 |= bit;
			if (smin.mini_bf1&bit) smin.mini_bf2 |= bit;
			smin.mini_bf1 |= bit;
			smin.critbf |= pat;
			smin.pairs27 |=mask^pat;
		}
		else smin.mini_triplet |= 1 << guas.ua3_imini[i81];
	}
	smin.SetMincount();
	if (p_cpt2g[9] < 10)smin.Status("exit second quick");
}
int BANDS_AB::FirstCheckGuaAll() {// quick check of critical situations
	for (; cur_ib < genb12.nband3; cur_ib++) {
		STD_B3 & b = genb12.bands3[cur_ib];
		b.SetUpMincountxy();
		if (b.smin.mincount < 7) return 0; // non critical (quick) band found
	}
	return 1;// all bands are in critical mode skip this xy
}
int BANDS_AB::EndCheckCurBand3() {// setup min count done 
	STD_B3 & b = genb12.bands3[cur_ib];

	// build and check start count per stack
	stack_countf.u64 = stack_count.u64 +b.smin.Count_per_stack().u64;
	if (stack_countf.u16[0] > 6 || stack_countf.u16[1] > 6 ||
		stack_countf.u16[2] > 6) return 1; // not ok

	// build free staks pattern for outfield uas
	uint32_t fstk = BIT_SET_27;
	if (stack_countf.u16[0] == 6) fstk ^= 07007007;
	if (stack_countf.u16[1] == 6) fstk ^= 070070070;
	if (stack_countf.u16[2] == 6) fstk ^= 0700700700;

	// Build in field out field and check minplus
	nuasb3_1 = nuasb3_2 = 0;
	int b3bloc = cur_ib >> 7, b3ir = cur_ib & 127;
	register int  Rfilt = b.smin.critbf;
	register uint32_t Rfstk = fstk;
	uint32_t  ua, nfree = 6 - b.smin.mincount;
	int iguan;
	b3_andout = fstk;

	BF128 vw = (b.g3x3y.g46 | b.g3x3y.g23)&tuguan.vxy;
	if (1) {
		b.g3x3y.g46.Print("active 46 for the band");
		vw.Print("active vw");
	}
	{// catch first all uas b3 from guas
		while ((iguan = vw.getFirst128()) >= 0) {
			vw.Clear(iguan);
			GUAN & gu = tuguan.tgua3x3y[iguan];
			if (gu.ncol == 2)ua = b.guas.ua_pair[gu.i81];
			else if (gu.ncol == 3)ua = b.guas.ua_triplet[gu.i81];
			else {
				TEMPGUAN4 & gu4 = tempguan4[gu.i81];
				if (!gu4.b3bf.On(b3bloc, b3ir))continue;
				ua = gu4.b3pat[cur_ib];
			}
			register uint32_t Ru = ua & BIT_SET_27;
			if (Ru & Rfilt) {
				Ru |= _popcnt32(Ru) << 27;
				AddUA32(uasb3_1, nuasb3_1, Ru);
			}
			else {
				Ru &= Rfstk;
				if (!Ru) return 1;//not ok  <<<<< vérifier 
				b3_andout &= Ru;
				Ru |= _popcnt32(Ru) << 27;
				AddUA32(uasb3_2, nuasb3_2, Ru);
			}
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
				if (!Ru)return 1;//not ok  <<<<< vérifier 
				b3_andout &= Ru;
				Ru |= _popcnt32(Ru) << 27;// restore count
				AddUA32(uasb3_2, nuasb3_2, Ru);
			}
		}
	}
	if (g17b.debug17)cout << nuasb3_1 << " " << nuasb3_2 << " if of final" << endl;
	if (nuasb3_2) {// check again the limit 6
		if (!nfree) return 1;//not ok
		b.smin.minplus++;
		if ((!b3_andout)) {
			if (nfree == 1) return 1;//not ok
			b.smin.minplus++;
		}
	}

	return 0; // ok to go
}

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
		int cell = tcluesb12[i];
		ws[cell] = genb12.grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[54 + i] = genb12.grid0[54 + i] + '1';
	fout1 << ws << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
	g17b.a_17_found_here++;

}
int  BANDS_AB::BuildUasB3_in(uint32_t known, uint32_t field) {
	nuas_in = 0;
	for (uint32_t iua = 0; iua < nuasb3_1; iua++) {
		register uint32_t  Ru = uasb3_1[iua];
		if (Ru & known) continue;// already hit, forget it
		Ru &= field;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(uas_in, nuas_in, Ru);
	}
	return 0;
}
void BANDS_AB::MergeUasExpand() {
	if (!nuasb3_1) {// will be if no mini row ua
		ExpandBand3(uasb3_2, nuasb3_2);// direct on outfiel
		return;
	}
	if (!nuasb3_2) {// should be processed in another way???
		ExpandBand3(uasb3_1, nuasb3_1);// direct in field
		return;
	}
	uint32_t tua[500], nua=0;
	struct UAT {// merging structue for both tables
		uint32_t *t, n, cc, ua,ind;
		void Init(uint32_t *te, uint32_t ne) {
			t = te; n = ne; ind = 0;
			GetNext();
		}
		void GetNext() {
			if (ind == n)cc = 100;
			else {
				uint32_t w = t[ind++];
				ua = w & BIT_SET_27;
				cc = w >> 27;
			}
		}
	}ua1,ua2;
	ua1.Init(uasb3_1, nuasb3_1);
	ua2.Init(uasb3_2, nuasb3_2);
	while (1) {
		if (ua1.cc < ua2.cc) {
			tua[nua++] = ua1.ua;
			ua1.GetNext();
		}
		else if (ua2.cc == 100) break;// end of the merge process
		else {
			tua[nua++] = ua2.ua;
			ua2.GetNext();
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


//================ part 2  band 3 processing

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
/*
void ZHOU::Guess17(int index, int diag) {
	if (diag) {
		char ws[82];
		cout << zh_g.pairs.String3X(ws)<< " pairs "<< endl;
		cout << zh_g2.triplets.String3X(ws) << " triplets " << endl;
		cout << zh_g2.quads.String3X(ws) << " quads " << endl;
	}
	if (zh_g2.s17_b3_mini) {// look once for mini rows
		zh_g2.s17_b3_mini = 0;
		uint32_t b3 = cells_unsolved.bf.u32[2], aig = 0;
		if (!b3)return; // if b3 solved, all is solved
		for (uint32_t i = 0, mask = 7; i < 9; i++, mask <<= 3) {
			uint32_t mini = b3 & mask;
			if (_popcnt32(mini) < 2) continue;
			aig = 1;
			//cout << Char27out(mini) << " mini " << endl;
			// try this mini row as unsolved
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			for (uint32_t j = 0, bit = 1; j < 27; j++, bit <<= 1) {
				if (bit&mini)continue;
				if (b3&bit) {
					int cell = j + 54, digit = zh_g2.grid0[cell];
					mynext->Seta_c(digit, cell);
				}
			}
			//cout << Char27out(mynext->cells_unsolved.bf.u32[2]) << " unsolved b3 " << endl;
			//cout << "appel suite zh_g.go_back="<< zh_g.go_back << endl;
			int ir = mynext->Full17Update();// solve as much as possible
			if (ir == 2) continue;
			uint32_t b3_n = mynext->cells_unsolved.bf.u32[2];
			if (!b3_n) continue; // now solved
			mynext->Guess17(0, 0);
			if (zh_g.nsol) return;
			zh_g.go_back = 0;// see why it is 1
		}
		if (aig) {
			int ir = Apply17SingleOrEmptyCells();// restore zh_g
			zh_g.go_back = 0;// see why it is 1
		}
	}
	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	{ // select band with more unsolved cells
		uint32_t nfreecells = 0, nw;
		if (w.bf.u32[0]) {
			nfreecells = _popcnt32(cells_unsolved.bf.u32[0]);
		}
		if (w.bf.u32[1]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[1]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
				}
			}
			else	nfreecells = _popcnt32(cells_unsolved.bf.u32[1]);
		}
		if (w.bf.u32[2]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[2]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
					w.bf.u32[1] = 0;
				}
			}
		}
	}
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell],
		tdig[10], ndig = 0;
	if (diag) {
		cout << "guess17 index=" << index << " cell " << cellsFixedData[cell].pt << endl;
	}
	// if first step try first false
	if(!index)	ClearCandidate_c(digit, cell);// force false
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig][0].On(xcell))tdig[ndig++] = idig;
	for (int idig = 0; idig < ndig; idig++) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(tdig[idig], cell, xcell);
		if (diag)cout <<"guess index="<<index<<" "
			<< tdig[idig]+1<< cellsFixedData[cell].pt << endl;
		mynext->Compute17Next(index+1,diag);
		if (zh_g.go_back) return;

	}
	if (!index) {
		FD[digit]->Set_c(cell);// restore the candidate
		SetaCom(digit, cell, xcell);
		if (diag)cout << "guess last index=" << index << " "
			<< digit + 1 << cellsFixedData[cell].pt << endl;
		Compute17Next(index, diag);

	}
}
*/
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
	diagh = 0;
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
	if (ir) {//consider store the fresh ua b3
		/*
		STD_B3 & myb3 = genb12.bands3[ib3];
		BF128 wua = zh_g2.cells_assigned;
		ua = wua.bf.u32[2];
		int cc = _popcnt32(ua);
		if (g17b.diag >= 2|| diagh) cout << Char27out(ua) << " ua b3 retour" << endl;
		//if (1) cout << Char27out(ua) << " ua b3 retour" << endl;
		if ( cc < 4) {// needs more tests on performance evolution
			if (cc == 2) {// fresh gua2
				int i81 = myb3.GetI81_2(ua);
				if (i81 >= 0) {
					bab.final81_2.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar2[i81]) bab.guar2i81[bab.nguared_2++] = i81;
					if (bab.ntuar2[i81] < GUAREDSIZE)
						bab.tuar2[i81][bab.ntuar2[i81]++] = wua.bf.u32[bab.ib];
					GEN_BANDES_12::SGUA2 & sg = genb12.tsgua2[i81];
					if (sg.nua < SIZETGUA)
						AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
			if (cc == 3) {// fresh gua2
				int i81 = myb3.GetI81_3(ua);
				if (i81 >= 0) {
					bab.final81_3.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar3[i81]) bab.guar3i81[bab.nguared_3++] = i81;
					if (bab.ntuar3[i81] < GUAREDSIZE)
						bab.tuar3[i81][bab.ntuar3[i81]++] = wua.bf.u32[bab.ib];
					GEN_BANDES_12::SGUA3 & sg = genb12.tsgua3[i81];
					if (sg.nua < SIZETGUA)
						AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
		}
		*/
	}
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
int G17B3HANDLER::BuildIfShortB3() {
	if (g17b.diag >= 2)	cout << Char27out(known_b3) << "entry BuildIfShortB3 " << endl;
	if (bands_ab.BuildUasB3_in(known_b3, active_b3))return 1;
	uasb3if = bands_ab.uas_in;
	nuasb3if = bands_ab.nuas_in;
	if (g17b.diag >= 2) {
		cout << "uas_in afer buildshort" << endl;
		for(uint32_t i=0;i< nuasb3if;i++)
			cout << Char27out(uasb3if[i]) << endl;
	}
	return 0;
}

int G17B3HANDLER::ShrinkUas1() {
	irloop = 0;
	uint32_t * tn = &uasb3if[nuasb3if], n = 0;
	register uint32_t Ra = active_b3,
		Rfilt = known_b3;
	for (uint32_t iua = 0; iua < nuasb3if; iua++) {
		register int Ru = uasb3if[iua];
		if (Ru & known_b3) continue;// already hit, forget it
		Ru &= active_b3;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else {
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
	}
	uasb3if = tn;
	nuasb3if = n;
	if (!n ) irloop = 0;// no need to loop again
	return 0;

}
void G17B3HANDLER::Go_Critical(uint32_t * wua){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	if (g17b.diag >= 2|| diagh)	cout << Char27out(known_b3) << "entry critical " << endl;
	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	if (!active_b3){
		bands_ab.FinalCheckB3(known_b3);
		return; // should be seen earlier if possible
	}
	//cout<< Char27out(active_b3) <<" active after 2 pairs"<<endl;smin.Status();
	int ua = IsMultiple(known_b3 | active_b3);
	if (ua) {
		if (wua) *wua &= ua;
		return;
	}
	if(BuildIfShortB3())return;
	if (ShrinkUas1()) return;// dead branch
	//cout << Char27out(active_b3) << " active aftershrinkuas" << endl; smin.Status();
	if (!active_b3) {
		//cout << Char27out(known_b3) << " call final check" << endl;  
		bands_ab.FinalCheckB3(known_b3);
		return; // should be seen earlier if possible
	}
	if (IsMultiple(known_b3 | active_b3))	return;// not valid using all cells
	if (irloop)		CriticalLoop();
	else CriticalExitLoop();
}

void G17B3HANDLER::CriticalLoop(){
	if (ShrinkUas1()) return;
	if (irloop)CriticalLoop();
	else CriticalExitLoop();
}
void G17B3HANDLER::CriticalExitLoop(){
	int nmissb =6 - _popcnt32(known_b3);// missing clues
	if (g17b.diag >= 2|| diagh) cout<<Char27out(known_b3) << "critical exit loop missb="<< nmissb 
		<<"nuasb3of="<< nuasb3of << endl;
	if (nmissb < 0)return;
	if (!active_b3){// nothing more to assign 
		if (nmissb)return;// dead cell in a mini row 3 pairs
		bands_ab.FinalCheckB3(known_b3);
		return;
	}
	// check known + active with brute force
	int wknown = known_b3 | active_b3;
	if (IsMultiple(known_b3 | active_b3))	return;// not valid using all cells
	if (g17b.diag > 2|| diagh) cout << "critical exit loop moreto do nuasb3=" << nuasb3if << endl;

	if (nuasb3if){		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_b3;
			for (uint32_t i = 0; i < nuasb3if; i++) {
				and_uas &= uasb3if[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (smin.mini_bf1) {	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, smin.mini_bf1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3 & mask;// catch the minirow
		}
		else {
			for (uint32_t i = 0; i < nuasb3if; i++) {
				register int ua = uasb3if[i]& active_b3,
					cc = _popcnt32(ua);
				if (cc < sizeua) { wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && smin.mini_triplet) {// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, smin.mini_triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 & mask;// catch the minirow

			}
		}

		while (bitscanforward(cell, wua)){
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop();
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void G17B3HANDLER::Critical_0_UA(){
	int nmissb = 6 - _popcnt32(known_b3);// missing clues
	if (g17b.diag >= 2) {
		cout << Char27out(known_b3) << "critical_0ua missb=" << nmissb << endl;
		smin.Status("b3 0 ua ");
	}
	if (nmissb < 0)return;
	if (!nmissb) {// nothing more to assign (granted at first call in a branch)
		bands_ab.FinalCheckB3(known_b3);
		return;
	}
	if (smin.mini_bf3) {// in active minirows with 3 pairs, assign 2
		while (smin.mini_bf3) {
			uint32_t mini;
			bitscanforward(mini, smin.mini_bf3);
			int shift = 3 * mini, bit = 1 << shift;
			smin.mini_bf3 ^= 1 << mini; //clear bit the mini row is always killed
			active_b3 &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++) {
				int * tpi = tp[i];
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (smin.mini_bf1) {// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, smin.mini_bf1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_b3 & mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		smin.mini_bf1 ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			if (x&bb) {
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (smin.mini_triplet) {// safety control should always be
		uint32_t mini;
		bitscanforward(mini, smin.mini_triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_b3 &= ~mask;// clear the minirow
		smin.mini_triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
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
		if (ShrinkUas1()) return;// dead branch
		if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
		if (irloop)		CriticalLoop();
		else CriticalExitLoop();
	}
}
void G17B3HANDLER::Go_Subcritical() {// nmiss to select in the critical field
	if (diagh) cout << Char27out(known_b3) << " entry Go_Subcritical()" << endl;
	active_b3 = active_sub = smin.critbf;
	// check first if a global solution  is still possible
	if (diagh) cout << Char27out(active_sub) << " active_sub Go_Subcritical()" << endl;
	if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	int cct = _popcnt32(smin.critbf) - smin.mincount;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	for (int ist = 0; ist < 3; ist++) {// check stacks 
		if (stack_count.u16[ist] > 5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
	if (diagh) cout << Char27out(active_sub) << " active_sub Go_Subcritical()" << endl;
	//if (BuildIfShortB3())return;
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
	uint32_t wua = wactive0;
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
				hn.Go_Subcritical();// to test later
			}
		}
	}
	if (nuasb3of) wua &= andoutf;
	if (diagh) cout << Char27out(wua) << " nmiss1 b3 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	hn.nmiss--;	hn.known_b3 |= bit;
			hn.Go_Critical();
		}
	}
}
void G17B3HANDLER::Go_miss2_b3() {
	ShrinkUasOfB3();// if known from up stream
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
				hn.Go_Subcritical();// to test later
			}
		}
	}
	uint32_t wua = wactive0;
	if (nuasb3of)wua &= uasb3of[0];// use first ua  
	if (diagh) cout << Char27out(wua) << " nmiss2 b3 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	hn.nmiss--;	hn.known_b3 |= bit;
			hn.Go_miss1_b3();
		}
	}
}
void G17B3HANDLER::Go_miss3_b3() {// always first direct entry
	uint32_t wua = wactive0 & uasb3of[0];// use first ua 
	if (diagh) cout << Char27out(wua) << " nmiss3 b3 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this;	hn.nmiss--;	hn.known_b3 |= bit;
			hn.Go_miss2_b3();
		}
	}
}



//================ debugging code
void G17Debug_CoutEntry(int * g) {
for (int i = 0; i < 54; i++) cout << g[i] + 1;
cout << "; entry canonical" << endl;
}
const char * g17tl1[10] = { "entry", "steps1", "steps2", "steps3", "steps4",
"common", "f1", "dead", "n51", "n62" };


const char * g17tl[40] = { "0__", "falluas", "fstack1", "UasGuas", "fb12v", "  b3v\t\t", "fnov", //0-6
"h3go", "critb", "ncritb", "m1go", "subgo", "subvv ",   //12
"gcrit", "gcvv", "gc_exit", "m2", "m2vv", "m1_ss", "19_", "fin check",//20
"21_", "22_", "23_", "24_", "25_", "26__", "gocheck", "check0", "nocheck0",//29
"30nua", "nua2", "nua3", "nb3", "__34", "vt1", "vt2", "", "T1", "T2",//39
};

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



void G17TB3GO::Debug() {
	cout << "status G17TB3GO ib3=" << ib3 << endl;
	//char ws[82];
	//cout << g17xy.bands_active_triplets.String3X(ws)
	//	<< "g17xy active triplets" << endl;
	cout << genb12.bands3[ib3].band << " \tbande3" << endl;
	cout << Char27out(pairs.u32[1]) << " in field bf" << endl << endl;
	cout << Char27out(pairs.u32[0]) << " pairs bf" << endl<<endl;
	cout << Char9out(count.u16[3]) << " all minis2" << endl;
	cout << Char9out(count.u16[0]) << "     minis2/1" << endl;
	cout << Char9out(count.u16[1]) << "     minis2/2" << endl;
	cout << Char9out(count.u16[2]) << "     minis2/3" << endl<<endl;
	cout << Char9out(minirows_triplets) << " mini triplets" << endl << endl;

	cout << countsum.u16[0] << "\t" << countsum.u16[1] << "\t" << countsum.u16[2]
		<< "\t total " << countsum.u16[3] << endl;
	cout << countstack.u16[0] << "\t" << countstack.u16[1] << "\t" << countstack.u16[2]
		<< "\t total " << countstack.u16[3] << endl;
}