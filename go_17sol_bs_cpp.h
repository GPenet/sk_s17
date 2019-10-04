//#define VALc15 13438
//#define VALc17 1032
//#define TESTXY 0
//2536435903110145  i1=2 i2=35
//#define TESTXY2 0
//22521159232789761
//536435903110145
//#define DEBUGLEVEL 10  nb12=3461507
//123456789456789123789132564268591437341627895597843216634278951815964372972315648;
//1....6...4......2.....3.5.....59........2...........16.....89.1..5.......72......; 189; 339; 341; 6


//___ start process expand bands collect uas guas ...
void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	if (aigstop)return;
	const char * diagband ="268591437341627895597843216";
	const char * diagpuz = "1....6...4......2.....3.5.....59........2...........16.....89.1..5.......72......";
	diag = diagbug = 0;
	if (diag )cout << "entry m10 nb12=" << genb12.nb12 << " nbands3=" << genb12.nband3 << " p_cpt2g[0]=" << p_cpt2g[0] << endl;
	if (genb12.nb12 == sgo.vx[6]) diagbug = 1;
	p_cpt2g[0] ++;
	p_cpt2g[1] +=genb12.nband3;
	if (genb12.nband3 > p_cpt2g[23])	p_cpt2g[23] = genb12.nband3;
	if (diag>1) {
		if (strcmp(diagband, myband2.band)) return;
		cout << "entry m10 nb12=" << genb12.nb12 << endl;
		cout << "this is the band in diag" << endl;
		if (diag == 2) {
			p17diag.SetAll_0();
			cout << diagpuz << " puz known" << endl;
			for (int i = 0; i < 81; i++) if (diagpuz[i] != '.')
				p17diag.Set_c(i);
		}
	}
	if (g17b.debug17) diag = diagbug = g17b.debug17;
	//______ true start band 2 expand
	tulock.Restore1();	
	myband2.ExpandBand();// expand band2
	if (!(myband1.nmybv5 | myband2.nmybv5)) return; // no 656 no 566
	tulock.LockExpand(myband2.nmyi3, myband2.nmybv5, myband2.nmybv6);
	if (g17b.debug17) {
		cout <<"b1 i416="<< myband1.i416 << " b2 i416=" << myband2.i416
			<< " n5b1=" << myband1.nmybv5 << " n5b2=" << myband2.nmybv5 << endl;
		if (0 &&g17b.npuz == 1025) {
			myband1.DebugExp();
			//myband1.DebugIndex();
			if (1) {
				cout << "band 1 partiel  6clues" << endl;
				uint32_t * t = myband1.mybv6,
					nd=7285,nf=7715;
				for (uint32_t i = nd; i < nf;i++)
					cout << Char27out(t[i]) << endl;
			}
			//myband1.Debug6();
		}
		if (GodebugFindKnown17()) {
			cout << "failed to find the expected 17 in expanded bands" << endl;
			return;
		}
		else cout << iband1 << " " << iband2 << " seen solution" << endl;

	}
	GoM10Uas();//expand bands 3  collect UAs 


	if (g17b.debug17 > 1) {
		if (DebugK17M10()) return;

		//g17b.a_17_found_here = 1;
		//return;
	}
	p_cpt2g[18] += genuasb12.nua;
	p_cpt2g[19] += genb12.ntua2;
	p_cpt2g[20] += genb12.ntua3;
	Go();// standard entry point for all 
}
void G17B::GoM10Uas() {
	int nb3 = genb12.nband3;
	//==== expand 6 for bands 3
	for (int ib3 = 0; ib3 < nb3; ib3++) {// lock only 6 clues 
		genb12.bands3[ib3].ExpandBand();// expand band3
		tulock.LockExpand(genb12.bands3[ib3].nmyi3, 0, genb12.bands3[ib3].nmybv6);
		genb12.bands3[ib3].DebugExp();
	}
	//=========================== collect UAs  old process 
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if (genuasb12.Initgen()) return;
	// extract small uas 
	tusmall.BuildSmallUas();
	genb12.BuildGang9x3();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	tuguan.Init();
	genb12.SecondSockets2Setup();// collect GUA2s 
	tuguan.ng2 = tuguan.nguan;
	genb12.SecondSockets3Setup();// collect GUA3s 
	tuguan.ng3 = tuguan.nguan- tuguan.ng2;

	// setupsockets common to all band3
	isguasocket2all = genb12.bands3[0].guas.isguasocket2;
	isguasocket3all = genb12.bands3[0].guas.isguasocket3;
	for (int ib3 = 1; ib3 < nb3; ib3++) {
		isguasocket2all &= genb12.bands3[ib3].guas.isguasocket2;
		isguasocket3all &= genb12.bands3[ib3].guas.isguasocket3;
	}
	// ===== add new guas four cells four columns 2 boxes 
	GoM10_guas_four_columns();
	tuguan.BuildCellKillersVector();
	//=== prepare the bands 3 valid guans and vectors
	tulock.InitBuf3();
	for (int ib = 0; ib < genb12.nband3; ib++) {
		tuguan.ApplyGuanToBand3(ib);
	}
}
void TU_SMALL::BuildSmallUas() {// find sub uas size 4 in the limit of 28 (64-36)
	nsmall1 = nsmall2 = 0;
	for (uint32_t i = 0; i < genuasb12.nua; i++) {
		uint64_t R = genuasb12.tua[i] & BIT_SET_2X;
		uint32_t	ua1 = (uint32_t)R,
			ua2 = (uint32_t)(R >> 32);

		uint32_t cc = _popcnt32(ua1);
		if (cc <9) 	Addsmall1(ua1 | (cc << 27));

		cc = _popcnt32(ua2);
		if (cc <9)	Addsmall2(ua2 | (cc << 27));
	}
	if (nsmall1>128) nsmall1=128;
	if (nsmall2 > 128) nsmall2 = 128;

	// clean the count in the small tables
	for (uint32_t i = 0; i < nsmall1; i++)
		tsmall1[i] &= BIT_SET_27;
	for (uint32_t i = 0; i < nsmall2; i++)
		tsmall2[i] &= BIT_SET_27;
	if (0) {
		cout << "table small bande1 n=" << nsmall1 << endl;
		for (uint32_t i = 0; i < nsmall1; i++)
			cout << Char27out(tsmall1[i]) << endl;;
		cout << "table small bande2 n=" << nsmall2 << endl;
		for (uint32_t i = 0; i < nsmall2; i++)
			cout << Char27out(tsmall2[i]) << endl;
	}
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


void TU_GUAN::BuildCellKillersVector() {// setup for <= 256 guan
	vv.Init(nguan);
	memset(vvcells, 255, sizeof vvcells);// all bits to 1
	uint32_t cc;
	for (uint32_t i = 0; i < nguan; i++) {// set uas
		register int bloc = i >> 7, ir = i & 127;
		register uint64_t Rw = tguan[i].killer;
		Rw &= BIT_SET_2X;
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			Rw ^= bit2;// clear bit
			vvcells[cc].Clear(bloc,ir);
		}
	}
}

void TU_GUAN::ApplyGuanToBand3(int ib3) {
	STD_B3 & b = genb12.bands3[ib3];
	b.tvv6 = tulock.pvx3;
	tulock.pvx3 += b.nmybv6;
	STD_B3::GUAs &g = b.guas;
	uint32_t * td = b.fd_sols[0];// per digit status
	uint32_t uas[256]; //  active guas pattern band 3
	VECT256 vcells[27];//cells vectors for guas


	memset(uas, 0, sizeof  uas);
	memset(&b.v_active_guas, 0, sizeof b.v_active_guas); //  non active
	memset(&b.vag2, 0, sizeof b.vag2); // pairs guan mode
	memset(&b.vag3, 0, sizeof b.vag3); // triplets
	BF128 sk2_4 = g.isguasocket2 | g.isguasocket4;
	for (uint32_t ig = 0; ig < nguan; ig++) {
		register int bloc = ig >> 7, ir = ig & 127;
		GUAN gu = tguan[ig];
		if (gu.ncol == 2) {// expected pair or GUA4 GUA6
			if (sk2_4.On_c(gu.i81)) {
				b.v_active_guas.Set(bloc, ir);
				uas[ig] = g.ua_pair[gu.i81] & BIT_SET_27;
				if (g.isguasocket2.On_c(gu.i81))b.vag2.Set(bloc, ir);
			}
		}
		else if (gu.ncol == 3) {// expected triplet
			if (g.isguasocket3.On_c(gu.i81)) {
				b.v_active_guas.Set(bloc, ir);
				b.vag3.Set(bloc, ir);
				uas[ig] = g.ua_triplet[gu.i81] & BIT_SET_27;
			}
		}
		else {// expected  2 digits  4 cells 2 rows 2 boxes 4 colums
			TEMPGUAN4 &tp = tempguan4[gu.i81];
			if (tp.b3bf.On(bloc, ir)) {//ok for this band
				b.v_active_guas.Set(bloc, ir);
				uas[ig] = tp.b3pat[ib3];
			}
		}
	}
	//_____________ setup the cells vector for the given  band 3
	memset(vcells, 255, sizeof vcells);
	uint32_t cc;
	for (uint32_t i = 0; i < nguan; i++) {// set uas
		register int bloc = i >> 7, ir = i & 127;
		register uint32_t Rw = uas[i];
		Rw &= BIT_SET_2X;
		while (bitscanforward(cc, Rw)) {// look for  possible cells
			Rw ^= 1 << cc;;// clear bit
			vcells[cc].Clear(bloc, ir);
		}
	}
	//____________ and apply cells to the 6 clues table 
	uint32_t ni3 = b.nmyi3 - 1;
	for (uint32_t i3 = 0; i3 < ni3; i3++) {// use first 3 clues from index 3
		int cc;
		uint32_t bf3 = b.myi3[i3].cellsbf;
		uint32_t indd = b.myi3[i3].ideb,
			indf = b.myi3[i3 + 1].ideb;

		VECT256 v3 = b.v_active_guas;
		register uint32_t bf3w = bf3;
		while (bitscanforward(cc, bf3w)) {// setup v3
			bf3w ^= 1 << cc;// clear bit
			v3.And(vcells[cc]);
			// _________________loop on remaining 3 clues
			for (uint32_t i6 = indd; i6 < indf; i6++) {
				uint32_t bf6 = b.mybv6[i6] & (~bf3);
				VECT256 v6 = v3;
				while (bitscanforward(cc, bf6)) {// setup v3
					bf6 ^= 1 << cc;// clear bit
					v6.And(vcells[cc]);
				}
				b.tvv6[i6] = v6;
			}
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

//_______________ processing bandA bandB (5clues 6 clues)
void BANDS_AB::Go(STD_B1_2 & ba, STD_B1_2 & bb, int i, int mode) {
	if (g17b.debug17)cout << "entry BANDS_AB.go i=" <<i << endl;
	tulock.InitBuf3();// init/reinit vector buffers use.
	mybb = &bb;
	ni3 = ba.nmyi3-1;
	mybv5 = ba.mybv5;
	myi3 = ba.myi3;
	ia = i; ib = 1 - i;
	tusmall.Build_B_Initial(ib);
	p_cpt2g[7] += ba.nmybv5*bb.nmybv6/100;
	if (g17b.debug17) {
		if (_popcnt32(g17b.p17diag.bf.u32[ia]) != 5) return;
		cout << "go a b ni3=" << ni3 << " ia=" << ia << endl;
	}
	mode_ab = mode;// 1 if must be 5 clues 
	stack_filter = 6;
	// loop on index 3 
	for (i3 = 0; i3 < ni3; i3++) {
		moreuas_AB.Init();
		wi3 = myi3[i3];
		bf3 = wi3.cellsbf;
		p_cpt2g[2]++;
		if (g17b.debug17) {//skip if not ok
			if (wi3.cellsbf & (~g17b.p17diag.bf.u32[ia])) continue;
			cout << Char27out(bf3) << "right i3=" << i3<<" ia="<<ia<< endl;
		}

		tusmall.BuildInit3cluesA(bf3);
		tuguan.InitB_Guas();
		//cout << " i3=" << i3 << " otheruas n=" << tusmall.n3_others_b << " ia=" << ia << endl;

		if (g17b.debug17 > 1 ) {
			cout  << " i3=" << i3	<< " otheruas n=" << tusmall.n3_others_b << endl;
			cout << Char64out(tusmall.smallpair.bf.u64[0]) << "smallpairs" << endl;
			cout << Char64out(tusmall.smalltriplet.bf.u64[0]) << "smalltriplet" << endl;
			cout << Char64out(tusmall.vsm3.bf.u64[0]);
			cout << Char64out(tusmall.vsm3.bf.u64[1]) << " vsm3" << endl;
		}
		moreuas_AB.Init();
		indf = myi3[i3 + 1].ideb5;
		//cout  << "i3=" << i3 << " inddeb="<< wi3.ideb5 
		//	<< " indf="<< indf << " size5=" <<indf- wi3.ideb5 << endl;


		nxbanda = 0;
		nxy_filt1 = 0;// init storage of XY passing first filter
		// _________________loop on remaining 2 clues
		for ( i5 = wi3.ideb5; i5 < indf; i5++) {
			p_cpt2g[3]++;
			bf5 = mybv5[i5];
			if (g17b.debug17) {//skip if not ok
				if (bf5 & (~g17b.p17diag.bf.u32[ia])) continue;
				cout << Char27out(bf5) << "right bandA i5=" << i5 << endl;
			}
			//if(i3==27 && ia==1)cout << Char27out(bf5) << " i5=" << i5 << endl;
			ncluesa = 0;
			BitsInTable32((int *)tclues, ncluesa, bf5);
			tusmall.BuildInit5cluesA(bf5);
			//if (i3 == 27 && ia == 1)cout  << " i5back init5" << endl;
			if (g17b.debug17) {//skip if not ok
				cout << Char27out(bf5) << " 5 clues count="
					<< tusmall.smin.mincount << "\tcountplus="
					<< tusmall.smin.minplus << endl;
				if(g17b.debug17 > 1)cout << Char64out(tusmall.vsm5.bf.u64[0]);
				if (g17b.debug17 > 1)cout << Char64out(tusmall.vsm5.bf.u64[1]) << " vsm5" << endl;
				if (tusmall.smin.mincount) {
					tusmall.smin.Status();
					if(tusmall.smin.minplus==6)
					fout1 << g17b.npuz << ";" << tusmall.smin.mincount
						<< ";" << tusmall.smin.minplus << endl;
				}
			}

			if (0 &&i5 == 2998 && ia == 1) {
				cout << Char27out(bf5) << " 5 clues count="
					<< tusmall.smin.mincount << "\tcountplus="
					<< tusmall.smin.minplus << endl;
				cout << Char64out(tusmall.vsm5.bf.u64[0]);
				cout << Char64out(tusmall.vsm5.bf.u64[1]) << " vsm5" << endl;
				if (tusmall.smin.mincount) 			tusmall.smin.Status();
			}
			//____ 
			if (tusmall.smin.minplus > 6) continue; //nothing to do
			if (tusmall.smin.minplus == 6) {// expand direct
				p_cpt2g[4]++;
				GoExpandB();
				// process ans call clean if limit buffer reached
				//if (Init3_5clues()) sbb.Go();
				//g17b.a_17_found_here = 1;
				//return;
			}
			else {// build a band A 5 clues to run with all 6 band B
				p_cpt2g[5]++;
				XBANDA & myxba = xbanda[nxbanda++];
				myxba.bf5 = bf5;
				myxba.vsm5 = tusmall.vsm5;
				myxba.vv5 = tusmall.vv3;// after first 3 cells
				register uint32_t w = bf5 ^ bf3;// 2 more cells
				int cca;// this is bandA
				while (bitscanforward(cca, w)) {
					w ^= 1 << cca; //clear bit
					myxba.vv5.And(tusmall.vv_cellsA[cca]);
				}

			}
			if (g17b.aigstop)return;
		}// end lop 3 to 5 clues in band A
		if (0 &&g17b.debug17) {// diag known
			if (bands_ab.ia != 1) return;
			cout << " check of expected XY ia=1" << endl;
			XBANDA w = xbanda[0];
			int i6 = g17b.iband1;
			cout << Char27out(w.bf5) << " bf5 i6 band B="<< i6 << endl;
			uint32_t kbf6 = mybb->mybv6[i6];
			cout << Char27out(kbf6) << " bf6 \n"  << endl;
			char ws[129];
			cout << w.vsm5.String128(ws) << "A vector 128"<<endl;
			cout << bB_v128[i6].String128(ws) << "B vector 128"<<endl;
		}
		if (0 &&i3 == 27 && ia == 1) {
			cout << "exit loop 5 nxbanda=" << nxbanda << endl;
		}
		if (nxbanda)Go_Matrix_X5_Y6();// now critical code, crossing and Y
		CleanTempXY();

		if (g17b.aigstop)return;
	}
}
void BANDS_AB::Go_Matrix_X5_Y6() {// now critical code, crossing and Y
	int locdiag = 0;
	if (g17b.debug17>1) {
		cout << "entry Go_Matrix_X5_Y6() nxbanda=" << nxbanda << endl;
	}
	//if (i3 == 27 && ia == 1) locdiag = 1;
#define XCHUNK 100
#define YCHUNK 200
	int ny6 = mybb->nmybv6;
	int ideby = 0, iendy = YCHUNK;
	if (iendy > ny6)iendy = ny6;
	if (locdiag)cout << "entry Go_Matrix_X5_Y6() nxbanda=" << nxbanda 
		<<"ny6=" <<ny6 << endl;
	while (ideby < ny6) { //Y chunk
		if (locdiag)cout << " ycollect() ideb="<<ideby<<"  n="<<iendy-ideby << endl;
		XBANDA ycollect[YCHUNK];// collect y data in a fix area
		{	register uint32_t *R6 = &mybb->mybv6[ideby];
			register BF128 * Rv128 = &bB_v128[ideby];
			register VECT256  * Rv256 = &bB_v256[ideby];
			register XBANDA * Ra = ycollect;
			for (int iy = ideby; iy < iendy; iy++,
				R6++, Rv128++, Rv256++, Ra++) {
				Ra->bf5 = *R6;
				Ra->vsm5 = *Rv128;
				Ra->vv5 = *Rv256;
			}
		}
		if (locdiag)cout << "exit ycollect()" << endl;
		uint32_t idebx = 0, iendx = XCHUNK;
		if (iendx > nxbanda)iendx = nxbanda;
		while (idebx < nxbanda) {// X chunk  
			register XBANDA * Ry = ycollect,*Ryend=&Ry[iendy-ideby], 
				*Rendx = &xbanda[iendx];
			XBANDA *rxd = &xbanda[idebx] ;
			for (; Ry < Ryend;  Ry++) {
				register XBANDA * Ra = rxd;
				for (; Ra < Rendx;  Ra++) {
					if ((Ra->vsm5 & Ry->vsm5).isEmpty()) {
						if (Ra->vv5.IsEmpty(Ry->vv5)) {
							// passing first filter, store it
							GINT64 & w = tempXY[nxy_filt1++];
							w.u32[0] = Ra->bf5;
							w.u32[1] = Ry->bf5;
						}
					}
				}
			}// end of a matrix XY
			if (locdiag)cout << "nxy_filt1=" << nxy_filt1 << endl;
			if (nxy_filt1 > 1000) CleanTempXY();
			idebx = iendx; iendx += XCHUNK;
			if (iendx > nxbanda)iendx = nxbanda;
		}
		ideby = iendy; iendy += YCHUNK;
		if (iendy > ny6)iendy = ny6;
	} 
	if (g17b.debug17>1) {
		cout << "exit Go_Matrix_X5_Y6() nxy_filt1=" << nxy_filt1 << endl;
	}	
	if (locdiag)cout << "exit Go_Matrix_X5_Y6() nxy_filt1=" << nxy_filt1 << endl;

}
void TU_SMALL::Build_B_Initial(int ib) {
	uint64_t *t = genuasb12.tua;
	uint32_t n = genuasb12.nua;
	uint32_t *tb = (ib) ? tsmall2 : tsmall1;// small is band B
	nsm = (ib) ? nsmall2 : nsmall1;// small is band B
	uint32_t lastsm = tb[nsm - 1],
		ccl = _popcnt32(lastsm);
	n_others_b = 0;
	smallpair.SetAll_0(); smalltriplet.SetAll_0();
	//_____ Init sm 
	for (uint32_t ism = 0; ism < nsm; ism++) {
		register uint32_t R = tb[ism];
		sm[ism].nua = 0;
		sm[ism].pat = R;
		uint32_t cc = _popcnt32(R);
		if (cc < 4) {
			if (cc == 3) {
				smalltriplet.Set(ism);
				small_count[ism] = 3;
			}
			else {
				smallpair.Set(ism);
				small_count[ism] = 2;
			}
			// and find the minirow 
			for (int im = 0, mask = 7; im < 9; im++, mask <<= 3) {// find the mini row
				if (R & mask) { small_imini[ism] = im; break; }
			}
		}
	}
	for (uint32_t iua = 0; iua < n; iua++) {
		register uint64_t U = t[iua] & BIT_SET_2X, Ua, Ub;
		if (!ib) { Ua = U >> 32; Ub = U & BIT_SET_27; }
		else { Ub = U >> 32; Ua = U & BIT_SET_27; }
		uint32_t countb = _popcnt32((uint32_t)Ub);
		if (countb < ccl || ((countb == ccl) && ((uint32_t)Ub <= lastsm))) {
			//cout << Char27out((uint32_t)Ub) << "\t";
			//cout << Char27out((uint32_t)Ua) << " to add" << endl;;

			// find the ism, must exist except for some uas added  
			// or if a subset was found
			for (uint32_t ism = 0; ism < nsm; ism++) {
				if (tb[ism] == (uint32_t)Ub) {
					sm[ism].AddMini3((uint32_t)Ua);
					goto nextiua;
				}
			}
			//cout << "pas trouve" << endl;

		}
		{// not in the sm list store in in ua/ub mode
			Ub <<= 32;
			if (n_others_b < 500)t_others_b[n_others_b++].u64 = Ua + Ub;
			continue;
		}
	nextiua:;
	}
	for (uint32_t ism = 0; ism < nsm; ism++) {// set all killers
		sm[ism].SetKiller(ism);
		//sm[ism].Debug1();
	}
	if (g17b.debug17)cout << "end Build_B_Initial  uas others n=" << n_others_b << endl;
	//______________ Build vectors for band A  and band B
	vsm = maskLSB[nsm];
	memset(vcellskill, 255, sizeof vcellskill);// all bits to 1
	memset(vcellsB, 255, sizeof vcellsB);// all bits to 1
	uint32_t cc;
	for (uint32_t i = 0; i < nsm; i++) {// set uas
		SM w = sm[i];
		register uint32_t Rw = w.killer& BIT_SET_27;
		while (bitscanforward(cc, Rw)) {// look for  possible cells
			register uint32_t bit2 = 1 << cc;
			Rw ^= bit2;// clear bit
			vcellskill[cc].clearBit(i);
		}
		Rw = w.pat& BIT_SET_27;
		while (bitscanforward(cc, Rw)) {// look for  possible cells
			register uint32_t bit2 = 1 << cc;
			Rw ^= bit2;// clear bit
			vcellsB[cc].clearBit(i);
		}
	}
	{// _____ apply vcellsB to valid band B 6
		STD_B1_2 & myb = (ib) ? myband2 : myband1;
		register uint32_t* R6 = myb.mybv6;
		register BF128* Rv = bB_v128;
		for (int i6 = 0; i6 < myb.nmybv6; i6++, R6++, Rv++) {
			*Rv = vsm;
			register uint32_t bf = *R6;
			int ccb;
			while (bitscanforward(ccb, bf)) {
				bf ^= 1 << ccb; //clear bit
				*Rv &= vcellsB[ccb];
			}
		}
	}
}
void TU_SMALL::BuildInit3cluesA(uint32_t bfA) {// reduce tables 
	bf3 = bfA;
	vsm3 = vsm;// initial status
	{// erase killed sm[] from cells
		register uint32_t bf = bfA;
		int cca;// this is bandA
		while (bitscanforward(cca, bf)) {
			bf ^= 1 << cca; //clear bit
			vsm3 &= vcellskill[cca];
		}
	}

	//_ shrink small sub tables clear empty 
	uint32_t * pbuf = buffer3;// initial reloading point
	int tu[128], ntu = vsm3.Table128(tu);
	register uint32_t F = bfA;
	for (int itu = 0; itu < ntu; itu++) {
		int ism= tu[itu];
		SM & smo = sm[ism];
		uint32_t *p = pbuf, np = 0;
		for (uint32_t iua = 0; iua < smo.nua; iua++) {
			register uint32_t U = smo.tua[iua];
			if (!(U&F))p[np++] = U;
		}
		if (!np)vsm3.clearBit(ism);
		else {
			sm3[ism].tua = p;
			sm3[ism].nua = np;
			pbuf += np;
		}
	}

	//___ shrink t_others
	n3_others_b=0;
	for (uint32_t ith = 0; ith < n_others_b; ith++)
		if (!(F & t_others_b[ith].u32[0]))
			if(n3_others_b<256)t3_others_b[n3_others_b++] = t_others_b[ith];
	//___________ build vectors for others apply to b 6clues
	vv_others.Init(n3_others_b);
	memset(vv_cellsA, 255, sizeof vv_cellsA);// all bits to 1
	memset(vv_cellsB, 255, sizeof vv_cellsB);// all bits to 1
	uint32_t cc;
	for (uint32_t i = 0; i < n3_others_b; i++) {// set uas
		register int bloc = i >> 7, ir = i & 127;
		register uint32_t Rw = t3_others_b[i].u32[0];//band A
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			Rw ^= 1 << cc;// clear bit
			vv_cellsA[cc].Clear(i, ir);
		}
		Rw = t3_others_b[i].u32[1] & BIT_SET_2X;//band B
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			Rw ^= 1 << cc;// clear bit
			vv_cellsB[cc].Clear(i, ir);
		}
	}
	//_______ prepare the 3 cells others vector
	vv3 = vv_others;// initial status
	{// erase killed sm[] from cells
		register uint32_t bf = bfA;
		int cca;// this is bandA
		while (bitscanforward(cca, bf)) {
			bf ^= 1 << cca; //clear bit
			vv3.And(vv_cellsA[cc]);
		}
	}
	//________ apply  vv_cellsB to b 6clues
	STD_B416 & myB = *bands_ab.mybb;
	uint32_t ni3 = myB.nmyi3 - 1;
	for (uint32_t i3 = 0; i3 < ni3; i3++) {// use first 3 clues from index 3
		int cc;
		uint32_t bf3 = myB.myi3[i3].cellsbf;
		uint32_t indd= myB.myi3[i3 + 1].ideb,
			indf = myB.myi3[i3 + 1].ideb;

		VECT256 v3 = vv_others;
		register uint32_t bf3w = bf3;
		while (bitscanforward(cc, bf3w)) {// setup v3
			bf3w ^= 1 << cc;// clear bit
			v3.And(vv_cellsB[cc]);
		}
		// _________________loop on remaining 3 clues
		for (uint32_t i6 =indd; i6 < indf; i6++) {
			uint32_t bf6 = myB.mybv6[i6] & (~ bf3);
			VECT256 v6 = v3;
			while (bitscanforward(cc, bf6)) {// setup v3
				bf6 ^= 1 << cc;// clear bit
				v6.And(vv_cellsB[cc]);
			}
			bB_v256[i6]=v6;
		}
	}
}
void TU_SMALL::BuildInit5cluesA(uint32_t bf2A) {
	int locdiag = 0;
	//if (bands_ab.i5 == 40)locdiag = 1;
	if (locdiag) cout << "BuildInit5cluesA diag" << endl;
	bf5 = bf2A;
	uint32_t bfplus = bf2A ^ bf3;
	vsm5 = vsm3;
	//cout << Char64out(vsm5.bf.u64[0]) << "vsm5 1" << endl;
	register uint32_t w = bfplus;
	{// find the final 128 small vector
		int cca;// this is bandA
		while (bitscanforward(cca, w)) {
			w ^= 1 << cca; //clear bit
			vsm5 &= vcellskill[cca];
		}
	}// now vsm5 has still possible small 
	if (locdiag) cout << "loca" << endl;
	register uint32_t F = bf2A;
	{// check if a small uab not killed still exists
		BF128 vsm5w = vsm5;
		int ism;
		while ((ism = vsm5w.getFirst128()) >= 0) {
			vsm5w.clearBit(ism);
			uint32_t * tua=sm3[ism].tua, nua = sm3[ism].nua;
			//cout << "ism=" << ism << " nua=" << nua << endl;
			for (uint32_t iua = 0; iua < nua; iua++) {
				if (!(tua[iua] & F)) {// 
					goto nextism; // first is ok
				}
			}
			vsm5.clearBit(ism);
		nextism:;
		}
	}// left with still valid small uab
	//cout << Char64out(vsm5.bf.u64[0]) << "vsm5 3" << endl;
	if (locdiag) cout << "locb" << endl;


	memset(&smin, 0, sizeof smin);
	{// compute mincount B using uas B 2 3
		BF128 vsm5w = vsm5;
		vsm5w &= (smallpair|smalltriplet);// small uas of interest
		if (vsm5w.Count() < 4) return; //will be vector mode
		int tu[128], ntu = vsm5w.Table128(tu);
		for (int itu = 0; itu < ntu; itu++) {// first pairs
			int ism = tu[itu];
			if (small_count[ism] != 2) continue;;
			uint32_t pat = sm[ism].pat,
				imini= small_imini[ism],
				bit=1<<imini;
			if (smin.mini_bf2&bit) smin.mini_bf3 |= bit;
			if (smin.mini_bf1&bit) smin.mini_bf2 |= bit;
			smin.mini_bf1 |= bit;
			smin.critbf |= pat;
			smin.pairs27 |= (7<<(3*imini))^pat;

		}
		for (int itu = 0; itu < ntu; itu++) {// then triplets
			int ism = tu[itu];
			if (small_count[ism] != 3) continue;;
			uint32_t pat = sm[ism].pat,
				imini = small_imini[ism];
			smin.mini_triplet |= 1 << imini;
		}
	}
	if (locdiag) cout << "locc smin.mincount="<< smin.mincount << endl;
	if (smin.mincount < 4) return; // will be vector mode
	if (smin.mincount >6) return; // will be killed
	// this is purely to see if a "critical" in + out is possible
	// if not, the process will use X;Y 5 clues 6 clues matrix vector checking
	uint32_t temp[1000], ntemp = 0;
	{// collect all uas in a temporary table
		{//first small
			int tu[128], ntu = vsm5.Table128(tu);
			for (int itu = 0; itu < ntu; itu++) {
				int ism = tu[itu];
				register uint32_t pat = sm[ism].pat& BIT_SET_27;
				pat |= _popcnt32(pat) << 27;
				AddUA32(temp, ntemp, pat);
			}
		}
		{// then band b must be here to expand band B
			uint32_t *t = bands_ab.mybb->tua, n = bands_ab.mybb->nua;
			for (uint32_t i = 0; i < n; i++) {
				register uint32_t pat = t[i] & BIT_SET_27;
				pat |= _popcnt32(pat) << 27;
				AddUA32(temp, ntemp, pat);
			}
		}
		{// and finally still valid "others"
			register uint32_t F = bf2A;
			for (uint32_t i = 0; i < n3_others_b; i++) {
				register uint64_t U = t3_others_b[i].u64;
				if ((uint32_t)U & F) continue;
				U >>= 32 ;// now ua bandb
				register uint32_t Ub = (uint32_t)U& BIT_SET_27;
				Ub |= (_popcnt32(Ub) << 27); //insert count
				AddUA32(temp, ntemp, Ub);
			}
		}
	}
	if (0) for (uint32_t itp = 0; itp < ntemp; itp++)
		cout << Char27out(temp[itp]) << endl;
	{// split temp > 3 cells in subtables in field out field
		uint32_t *ti = bands_ab.btuaif, *to = bands_ab.btuaof,
			ni = 0, no = 0,ando= BIT_SET_27;
		register uint32_t F = smin.critbf; // in field for this band A
		for (uint32_t itp = 0; itp < ntemp; itp++) {
			register uint32_t ua = temp[itp];
			if ((ua >> 27) < 4) continue;// pairs triplets skipped
			if (ua&F) {if(ni<200) AddUA32(ti, ni, ua); }// smaller are enough 
			else {
				AddUA32(to, no, ua);// here keep all for the and
				ando &= ua;
			}
		}
		if (no) {// adust the mincount if outfield uas found
			smin.minplus++;
			if(!ando)smin.minplus++; // min 2 if no killer
		}
		//cout << Char27out(F) << "ntemp=" << ntemp << " ni=" << ni << " no=" << no << endl;
		bands_ab.nbif = ni;
		bands_ab.nbof = no;
		bands_ab.andoutf = ando;
	}
}

void BANDS_AB::BANDB::DebugIfOf() {
	cout << "In Field Bandb table nuaif= "<<nuaif << endl;
	for (uint32_t i = 0; i < nuaif; i++)
		cout << Char27out(tuaif[i]) << endl;
	cout << "Out Field Bandb table nuaof="<< nuaof << endl;
	for (uint32_t i = 0; i < nuaof; i++)
		cout << Char27out(tuaof[i]) << endl;

}

//======================== start band b expansion

void BANDS_AB::GoExpandB() {//start band b expansion
	nmiss = 6 - tusmall.smin.mincount;
	if (g17b.debug17 || diagbug) if(GoExpandBDebugInit()) return;
	//_______________prepare first  bandB
	sbb.Init(btuaif, btuaof, nbif, nbof, andoutf);
	sbb.Go();
	nmoreif = nmoreof = 0;
}
int BANDS_AB::GoExpandBDebugInit() {
	if (diagbug) {
		cout << "start B expansion nmiss=" << nmiss
			<< " nuaif=" << nbif << " nuaof=" << nbof << "\t p_cpt2g[4]=" << p_cpt2g[4] << endl;
	}
	if (g17b.debug17) {
		cout << Char27out(bf5)
			<< "we have the right band A nmiss=" << nmiss
			<< " nuaif=" << nbif << " nuaof=" << nbof << "\t p_cpt2g[4]="
			<< p_cpt2g[4] << endl;
	}
	return 0;
}

void BANDS_AB::BANDB::Go() {//start band b expansion in recursive mode
	mB = tusmall.smin;
	nmiss = 6 - mB.mincount;
	diag = 0;//GoExpandBDebugSetDiag();

	//_______________
	known_bb = rknown_bb = 0;
	ndead =  BIT_SET_27;
	active_bb = mB.critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_bb;//  active cells out field
	if (!nmiss)Go_Critical();
	else if(nmiss==1) Go_miss1();
	else Go_miss2();
}
void BANDS_AB::BANDB::GoExpandBDebugSetDiag() {
	if (0 && (!g17b.debug17) && p_cpt2g[3] == 3180) {
		cout << "bandb in diag p_cpt2g[3]=" << p_cpt2g[3] << endl;
		DebugIfOf();
		diag = 1;
	}
}

//_____________________ critical
void BANDS_AB::BANDB::CriticalAssignCell(int Ru) {// assign a cell within the critical cells
	// Ru is usually a register containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_bb |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & mB.mini_bf3) {// cell is in minirow with 3 pairs active
		active_bb &= ~Ru; //clear the cell
		mB.mini_bf3 ^= bit; // now only a pair to hit
		mB.mini_bf1 |= bit;
	}
	else {// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_bb &= (~Mask); // kill the minirow as active
		mB.mini_bf1 &= ~bit;
		mB.mini_triplet &= ~bit;
	}
}
int BANDS_AB::BANDB::IsFinalOrMultiple(uint32_t * wua) {
	diag = 0;
	//if (g17b.npuz == 2007) diag = 1;
	if(diag)cout << " entry multiple"  << endl;
	// all if uas must be hit
	register int bf = known_bb | active_bb;
	for (uint32_t iua = 0; iua < nuaif; iua++) {
		register int Ru = tuaif[iua];
		if(! (Ru & bf)) return 1;//known UA not hit multiple 
	}
	if (diag)cout << "suite3 after uas check" << endl;
	if (g17b.moreb.Check(bf)) return 1;// more hit
	//cout << "suite4 after more" << endl;
	if (!active_bb) {
		bands_ab.EnterTempxy(known_bb);
		return 1; 
	}
	if (diag)cout << "suite5 call is multiple" << endl;
	if (bf != rknown_bb) {
		rknown_bb = bf;
		if(bands_ab.IsMultiple_bb(rknown_bb, diag))return 1;
	}
	if (diag)cout << "suite6 retour is multiple" << endl;
	return 0;
}
int BANDS_AB::BANDB::BuildIF_short() {
	BANDS_AB & bab = bands_ab;
	tuaif = bab.tuaif;
	nuaif = 0;
	for (uint32_t iua = 0; iua < bab.nbif; iua++) {
		register int Ru = bab.btuaif[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(tuaif, nuaif, Ru);
	}
	for (uint32_t iua = 0; iua < bab.nmoreif; iua++) {
		register int Ru = bab.more_if[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
		if (!Ru) return 1;// dead branch
		Ru |= _popcnt32(Ru) << 27;
		AddUA32(tuaif, nuaif, Ru);
	}
	// also more if pending
	return 0;
}

int BANDS_AB::BANDB::ShrinkUas1() {
	irloop = 0;
	uint32_t * tn = &tuaif[nuaif], n = 0;
	for (uint32_t iua = 0; iua < nuaif; iua++) {
		register int Ru = tuaif[iua];
		if (Ru & known_bb) continue;// already hit, forget it
		Ru &= active_bb;
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
	tuaif = tn;
	nuaif = n;
	if (!nuaif) irloop = 0;// no need to loop again
	return 0;
}
void BANDS_AB::BANDB::Go_Critical(uint32_t * wua) {// critical situation all clues in pairs tripl:ets
	if (g17b.debug17>2 || diag) {
		cout << Char27out(active_bb) << "B active bb entry critical" << endl;
		cout << Char32out(known_bb) << "known 32 at entry" << endl;
	}
	if (mB.mini_bf2) {// assign common cell 
		for (int i = 0, bit = 1, mask = 7; i < 9; i++ , bit <<= 1, mask <<= 3) {
			if (mB.mini_bf2&bit) {
				active_bb &= (~mask); // clear the minirow
				known_bb |= mask & (~mB.pairs27);// and set the common cell a			
			}
		}
	}
	mB.mini_bf2 = 0;
	if (BuildIF_short()) return;// shrink and sort the table in field
	if (ShrinkUas1()) return;//assign possibles
	if (IsFinalOrMultiple(wua)) return;
	if (irloop)		CriticalLoop();
	else CriticalExitLoop();
}
void BANDS_AB::BANDB::CriticalLoop() {
	if (ShrinkUas1()) return;
	if (irloop)CriticalLoop();
	else CriticalExitLoop();
}
void BANDS_AB::BANDB::CriticalExitLoop() {
	int nmissb = 6 - _popcnt32(known_bb);// missing clues
	if (g17b.debug17 > 2 || diag) {
		cout << Char27out(known_bb) << "B entry critical exit loop nmissb="<< nmissb << endl;
	}
	if (nmissb < 0)return;
	if (IsFinalOrMultiple())return;
	//cout << "suite1 nuaif="<<nuaif << endl;
	if (nuaif) {		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_bb;
			for (uint32_t i = 0; i < nuaif; i++) {
				and_uas &= tuaif[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (mB.mini_bf1) {	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, mB.mini_bf1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_bb & mask;// catch the minirow
		}
		else {
			for (uint32_t i = 0; i < nuaif; i++) {
				register int ua = tuaif[i]& active_bb,
					cc = _popcnt32(ua);
				if (cc < sizeua) { wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && mB.mini_triplet) {// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, mB.mini_triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_bb & mask;// catch the minirow
			}
		}
		if (g17b.debug17 > 2 || diag)		
			cout << Char27out(wua) << " wua to use" << endl;

		while (bitscanforward(cell, wua)) {
			if (0)		
				cout << Char27out(wua) << " wua entry while nmissb=" 
				<< nmissb 		<<" nuaif="<<nuaif<< endl;
			register int bit = 1 << cell;
			wua ^= bit;// clear bit			
			active_bb ^= bit;//  now a dead cell downstream
			BANDS_AB::BANDB hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop();
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void BANDS_AB::BANDB::Critical_0_UA() {
	int nmissb = 6 - _popcnt32(known_bb);// missing clues
	if (g17b.debug17 > 2) {
		cout<<Char27out(known_bb) << "B entry critical 0_ua nmissb=" << nmissb << endl;
	}
	if (nmissb < 0)return;
	if (!nmissb) {// nothing more to assign (granted at first call in a branch)
		bands_ab.EnterTempxy(known_bb);
		return;
	}
	if (mB.mini_bf3) {// in active minirows with 3 pairs, assign 2
		while (mB.mini_bf3) {
			uint32_t mini;
			bitscanforward(mini, mB.mini_bf3);
			int shift = 3 * mini, bit = 1 << shift;
			mB.mini_bf3 ^= 1 << mini; //clear bit the mini row is always killed
			active_bb &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++) {
				int * tpi = tp[i];
				BANDS_AB::BANDB hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (mB.mini_bf1) {// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, mB.mini_bf1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_bb & mask;// catch the minirow
		active_bb &= ~mask;// clear the minirow
		mB.mini_bf1 ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			if (x&bb) {
				BANDS_AB::BANDB hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (mB.mini_triplet) {// safety control should always be
		uint32_t mini;
		bitscanforward(mini, mB.mini_triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_bb &= ~mask;// clear the minirow
		mB.mini_triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			BANDS_AB::BANDB hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
//_____________________ not critical
void BANDS_AB::BANDB::ShrinkUasOf() {
	if (known_bb) {// shrink the out field table
		uint32_t * tn = &tuaof[nuaof], n = 0;
		register uint32_t Ra = wactive0,
			Rfilt = known_bb;
		andoutf = BIT_SET_27;
		for (uint32_t iua = 0; iua < nuaof; iua++) {
			register int Ru = tuaof[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		for (uint32_t i = 0; i < bands_ab.nmoreof; i++) {// insert moreof 
			register int Ru = bands_ab.more_of[i];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			andoutf &= Ru;
			Ru |= _popcnt32(Ru) << 27;
			AddUA32(tn, n, Ru);
		}
		tuaof = tn;
		nuaof = n;
	}
}

void BANDS_AB::BANDB::Go_miss1() {//always one more ua if miss2 called minplus=6
	ShrinkUasOf();// if known from up stream   
	if (!andoutf) return;// no solution 2 digits 
	uint32_t wua = wactive0;
	wua &= andoutf;
	if (bands_ab.more_of) {
		for (uint32_t i = 0; i < bands_ab.nmoreof; i++) {// try also more(s)
			if (bands_ab.more_of[i] & known_bb)continue;
			wua &= bands_ab.more_of[i];
		}
	}

	if (g17b.debug17 > 2) 		cout << Char27out(wua) << " nmiss1 wua" << endl;
	if (wua) { // apply first UA to use or all out field cells
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			hn.Go_Critical(&wua);
		}
	}
}
void BANDS_AB::BANDB::Go_miss2() {
	uint32_t wua = wactive0 & tuaof[0];// use first ua  
	{ // apply first UA to use  
		uint32_t res;
		while (bitscanforward(res, wua)) {
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			BANDB hn = *this;	hn.nmiss--;	hn.known_bb |= bit;
			hn.Go_miss1();
		}
	}
}


//______________band B buffer of XY
void BANDS_AB::EnterTempxy(uint32_t bf) {// from direct expansion
	GINT64 & w = tempXY[nxy_filt1++];
	w.u32[0] = tusmall.bf5;
	w.u32[1] = bf;
}
void BANDS_AB::CleanTempXY() {//GINT64 tempXY[15000];
	if (g17b.debug17>1)cout << "entry CleanTempXY() ntempx=" << nxy_filt1 << endl;
	p_cpt2g[6]+= nxy_filt1;
	for (uint32_t itemp = 0; itemp < nxy_filt1; itemp++) {
		GINT64  w = tempXY[itemp];
		//cout << Char2Xout(w.u64) << endl;
		if (g17b.debug17_check) {
			if (w.u32[0] != g17b.k17x(ia)) continue;
			if (w.u32[1] != g17b.k17x(ib)) continue;
			cout << Char27out(w.u32[1]) << "we have the right band B" << endl;
		}

		bfA = w.u32[0]; bfB = w.u32[1]; 
		// check the "more" table if available
		if (moreuas_AB.Check(w.u64)) continue;
		// rebuild table of clues and build stack count for A,B
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
		tuguan.InitC_Guas();// get final gangster C vector vB
		cur_ib = 0;
		int ir3 = GetNextIb3();// return -1 if nothing, 0 if matrix
		if (ir3<0) continue;// no band 3 to process
		p_cpt2g[8]++;
		//_____ check for unique solution bands A;B
		register uint64_t myua = zh2b[0].ValidXY(tclues, nclues, 0);
		if (myua) {//___ if not do the best with the UA
			/*
			if (g17b.debug17 > 2) cout << Char2Xout(myua) << "uaret band b multiple" << endl;
			register uint64_t uab = myua >> 32, uaa = myua & BIT_SET_27, uab12 = myua;
			if (ia)uab12 = uab | (uaa << 32);
			uint64_t cc64 = _popcnt64(uab12);

			if (cc64 < 12) {// this should never be check for a bug
				cout << Char2Xout(uab12) << "ua < 12 to add ??? g17b.npuz" << g17b.npuz
					<< endl;
				cerr << "ua < 12 to add ??? g17b.npuz" << g17b.npuz << endl;
				zh2b[0].DebugValidXY(tclues, nclues, 0);
			}
			moreuas_AB.Add(myua);
			// insert it in the right place in tusmall (next 3 clues)
			tusmall.InsertMore(myua);
			if (cc64 < 20)genuasb12.AddUACheck(uab12 | ((uint64_t)cc64 << 59));// and update the UA table
			*/
			continue;
		}// end not a unique solution
		// init full check A B reordered if needed
		memcpy(tcluesb12, tclues, 11 * 4);
		if (ia) for (int i = 0; i < 11; i++)
			if (i > 4)tcluesb12[i] -= 27; else tcluesb12[i] += 27;
		if (zhou[0].PartialInitSearch17(tcluesb12, nclues))// A;B order ok to check
			return;// would be  bug

		p_cpt2g[9]++;
		if (ir3 < 6) {// apply a matrix filter
			p_cpt2g[10]++;
			MatrixB3(genb12.bands3[cur_ib]);
		}
		else {// do a band3 expansion critical miss1 miss2
			p_cpt2g[11]++;
			GoBand3();
			fout1 << g17b.npuz << endl;
			g17b.a_17_found_here = 1;// pour éviter le blocage
		}
	}
	nxy_filt1 = 0;
}
void BANDS_AB::MatrixB3( STD_B3 &  myb3) {
	int diagloc = 0;
	moreuas_b3.Init();
	moreuas_b3_small.Init();
	memcpy(&genb12.grid0[54], myb3.band0, 4 * 27);
	BF128 v1 = tuguan.vB.v[0], v2 = tuguan.vB.v[1];
	register VECT256 * Rd = myb3.tvv6, *Rf = &Rd[myb3.nmybv6];
	register uint32_t* Rb = myb3.mybv6;

	for (; Rd < Rf; Rd++, Rb++) {
		if ((v1 & Rd->v[0]).isNotEmpty())continue;
		if ((v2 & Rd->v[1]).isNotEmpty())continue;
		register uint32_t Bf = *Rb;
		// should check the stack limit
		FinalCheckB3(Bf);
	}
}
void TU_SMALL::InsertMore(uint64_t ua) {
	SM sml = sm[nsm - 1]; // last small
	uint32_t uab = (uint32_t)(ua >> 32), ccb = _popcnt32(uab)
		,patlast=sml.pat&BIT_SET_27,ccl= _popcnt32(patlast);
	if (ccb <= ccl) {// should be in a small
		for(uint32_t i=0;i<nsm;i++)
			if ((sm[i].pat &BIT_SET_27) == uab) {// right small
				sm[i].AddMini3((uint32_t)ua);// insert uaA
				return;
			}
	}
	// not a small, put it in others
	if (n_others_b < 1000)t_others_b[n_others_b++].u64 = ua;
}
void TU_GUAN::InitB_Guas() {// guas filter => already killed/forced plus sub table

	// build new subtables still active not fixed
	pguabufr = guabufr;
	// first check reduced vector to get potential active guan
	vA = vv;
	//vA.Print("vA initial") ;
	uint32_t bf3 = bands_ab.bf3;
	{// erase killed sm[] from cells
		register uint32_t bf = bf3;
		int cca;// this is bandA
		while (bitscanforward(cca, bf)) {
			bf ^= 1 << cca; //clear bit
			if (bands_ab.ia)cca += 27;
			vA.And(vvcells[cca]); 
		}
	}
	//vA.Print("vA after kill");

	// put vA in table
	int tu[256], ntu;   vA.Table(tu,ntu);
	// for each potential active guan, build reduced table
	register uint64_t F = bf3;// pattern 3 clues in band A
	if (bands_ab.ia)F <<= 32; // ua mode 2X 
	for (int itu = 0; itu < ntu; itu++) {
		int iguan = tu[itu];
		register int bloc = iguan >> 7, ir = iguan & 127;
		GUAN & gu = tguan[iguan];
		uint64_t *p = pguabufr;
		uint32_t np = 0;
		for (uint32_t iua = 0; iua < gu.nua; iua++) {
			register uint64_t U = gu.tua[iua] ;
			if (U&F) continue; // hit by banda
			AddUA64(p, np, U);
		}
		if (!np) vA.Clear(bloc, ir); // no ua kill it
		else {// lock the stored Uab
			gu.tuar = pguabufr;
			gu.nuar = np;
			pguabufr = &pguabufr[np];
		}
	}
	//vA.Print("vA final") ;
}
void TU_GUAN::InitC_Guas() {// guas vector after 6 clues B
	// first check reduced vector to get potential active reduced guan
	//cout << "entry InitC_Guas()"  << endl;
	vB = vA;
	//vB.Print("vB initial status");
	for (int icl = 3; icl < 11; icl++) {// remaining clues 
		int clue = bands_ab.tclues[icl];
		if (bands_ab.ia) {// must exchange bands 1 2
			if (icl < 5)clue += 27;
			else clue -= 27;
		}
		vB.And(vvcells[clue]);
	}
	//vB.Print("vB after killers");
	// put vB in table
	int tu[256], ntu;   vB.Table(tu,ntu);

	// for each potential active guan, check reduced table
	register uint64_t F ;// pattern in band B
	if (bands_ab.ia)F = (uint64_t)bands_ab.bfB + ((uint64_t)bands_ab.bfA << 32);
	else F = (uint64_t)bands_ab.bfA + ((uint64_t)bands_ab.bfB << 32);
	//cout <<  Char2Xout(F)<<"filter for guas final"<<endl;
	for (int itu = 0; itu < ntu; itu++) {
		int iguan = tu[itu];
		//if(iguan<64)cout << "iguan=" << iguan << endl;
		register int bloc = iguan >> 7, ir = iguan & 127;
		GUAN & gu = tguan[iguan];
		for (uint32_t iua = 0; iua < gu.nuar; iua++) {
			register uint64_t U = gu.tuar[iua] ;
			//if (iguan < 64)cout <<Char2Xout(U) << endl;
			if (!(U&F)) goto nextitu; // first forces ok
		}// if all hits,clear the bit
		vB.Clear(bloc, ir);
	nextitu:;
	}	// not cleared bits are the list of sockets to consider

	if (g17b.debug17 > 1)vB.Print("vB final status");
}
int BANDS_AB::GetNextIb3() {//check tuguan.vB => mincount per band 
	int locdiag = 0;
	if(locdiag)cout << "entry GetNextIb3() cur_ib=" <<cur_ib<< endl;
	if (locdiag)tuguan.DebugvB();
	for (; cur_ib < genb12.nband3; cur_ib++) {// look for first band
		STD_B3 & sb3 = genb12.bands3[cur_ib];
		sb3.vag_AB = sb3.v_active_guas;
		sb3.vag_AB.And(tuguan.vB);
		sb3.SetUpMincountAB();// building smin not minplus 

		sb3.smin.Status();
		sb3.smin.minplus = sb3.smin.mincount;

		if (sb3.smin.mincount > 6) continue;
		// build and check start count per stack
		stack_countf.u64 = stack_count.u64 +
			sb3.smin.Count_per_stack().u64;

		if (g17b.debug17)cout << stack_countf.u16[0] << stack_countf.u16[1] << stack_countf.u16[2]
			<< " stack count" << endl;

		if (stack_countf.u16[0] > 6 || stack_countf.u16[1] > 6 ||
			stack_countf.u16[2] > 6) continue;

		if (sb3.smin.mincount < 4) return 0; // will be matrix expansion
		// build free staks pattern
		uint32_t fstk= BIT_SET_27;
		if (stack_countf.u16[0] == 6) fstk ^= 07007007;
		if (stack_countf.u16[1] == 6) fstk ^= 070070070;
		if (stack_countf.u16[2] == 6) fstk ^= 0700700700;

		// Build in field out field and check minplus
		nuasb3_1 = nuasb3_2 = 0;
		int b3bloc = cur_ib >> 7, b3ir = cur_ib & 127;
		VECT256 vw = tuguan.vB; vw.Not(sb3.vag2); vw.Not(sb3.vag3);
		if (locdiag)vw.Print("sockets not pair or triplet");
		if (locdiag)cout << Char27out(fstk) << " out free stacks" << endl;
		int tu[256], ntu;
		vw.Table(tu, ntu);// now guan to consider as uas
		register int  Rfilt = sb3.smin.critbf;
		register uint32_t Rfstk = fstk;
		int * ua_46 = sb3.guas.ua_pair; // ua 4/6 patterns
		uint32_t ua,andout= fstk,nfree=6-sb3.smin.mincount;
		{// first sockets 2 and sockets 4
			for (int itu = 0; itu < ntu; itu++) {
				int iguan = tu[itu];
				GUAN & gu = tuguan.tguan[iguan];
				if (gu.ncol == 2)ua = ua_46[gu.i81];
				else  {
					TEMPGUAN4 & gu4 = tempguan4[gu.i81];
					if (!gu4.b3bf.On(b3bloc, b3ir))continue;
					ua =gu4.b3pat[cur_ib];
				}
				register uint32_t Ru = ua & BIT_SET_27;
				if (locdiag)cout << Char27out(Ru) << "ua4/6 or 4 cols" << endl;
				if (Ru & Rfilt) {
					Ru |= _popcnt32(Ru) << 27;
					AddUA32(uasb3_1, nuasb3_1, Ru);
				}
				else {
					Ru &= Rfstk;
					if (!Ru) goto nextband;
					andout &= Ru;
					Ru |= _popcnt32(Ru) << 27;
					AddUA32(uasb3_2, nuasb3_2, Ru);
				}
			}
		}
		if (locdiag)cout << nuasb3_1 << " " << nuasb3_2 << " if of after guas" << endl;
		if (nuasb3_2) {// check the limit 6
			if (!nfree) continue;
			if ((!andout) && nfree == 1) continue;
		}
		{
			// add now band 3 uas
			uint32_t * to = sb3.tua;
			for (uint32_t i = 0; i < sb3.nua; i++) {
				register uint32_t Ru = to[i] & BIT_SET_27;
				if (locdiag)cout << Char27out(Ru) << "uas belonging to band 3" << endl;
				Ru |= _popcnt32(Ru) << 27;
				if (Ru & Rfilt) {
					Ru |= _popcnt32(Ru) << 27;
					AddUA32(uasb3_1, nuasb3_1, Ru);
				}
				else {
					Ru &= Rfstk;
					if (!Ru) goto nextband;
					andout &= Ru;
					Ru |= _popcnt32(Ru) << 27;
					AddUA32(uasb3_2, nuasb3_2, Ru);
				}
			}

		}
		if (g17b.debug17)cout << nuasb3_1 << " " << nuasb3_2 << " if of final" << endl;
		if (nuasb3_2) {// check again the limit 6
			if (!nfree) continue;
			sb3.smin.minplus++;
			if ((!andout)) {
				if( nfree == 1) continue;
				sb3.smin.minplus++;
			}

		}
		if (locdiag)cout  << "final value sb3.smin.minplus="<< sb3.smin.minplus << endl;

		return sb3.smin.minplus;// no other possibility to reject this band, exit
	nextband:;
	}
	return -1;// return value if no band 3 to process
}
void STD_B3::SetUpMincountAB() {
	memset(&smin, 0, sizeof smin);
	//vag2.Print("band vag2");
	//vag3.Print("band vag3");
	//vag_AB.Print("band vag_AB");

	{// compute mincount pairs
		BF128 wpair = vag_AB.v[0] & vag2.v[0];// pairs always in the first 128
		int tu[128], ntu = wpair.Table128(tu);
		//cout << "expecting n=" << ntu << " pairs" << endl;
		for (int itu = 0; itu < ntu; itu++) {// first pairs
			int iguan = tu[itu], i81 = tuguan.tguan[iguan].i81;
			uint32_t pat = guas.ua_pair[i81],
				imini = guas.ua2_imini[i81],
				bit = 1 << imini;
			//cout << Char27out(pat) << " iguan=" << iguan << " i81=" << i81 
			//	<<" imini="<<imini<< endl;
			if (smin.mini_bf2&bit) smin.mini_bf3 |= bit;
			if (smin.mini_bf1&bit) smin.mini_bf2 |= bit;
			smin.mini_bf1 |= bit;
			smin.critbf |= pat;
			smin.pairs27 |= (7 << (3 * imini)) ^ pat;

		}
	}
	{// compute mincount triplets can be with index>128
		int tu[128], ntu;
		VECT256 vw = vag_AB; vw.And(vag3);		vw.Table(tu,ntu);
		for (int itu = 0; itu < ntu; itu++) {// then triplets
			int iguan = tu[itu], i81 = tuguan.tguan[iguan].i81;
			uint32_t pat = guas.ua_triplet[i81],
				imini = guas.ua3_imini[i81],
				bit = 1 << imini;
			smin.mini_triplet |= 1 << imini;
		}
		smin.SetMincount();
	}	
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




int BANDS_AB::IsMultiple_bb(int bf,int diag) {
	if (_popcnt32(bf) > 25) return 0;	
	nclues = ncluesa;// buil bandb part of tclues
	uint32_t cellb, wbf=bf;
	while( wbf){
		bitscanforward(cellb, wbf);
		wbf ^= 1 << cellb;
		tclues[nclues++] = cellb + 27;
	}
	//p_cpt2g[7]++;
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues, diag);
	if (myua) {//store the fresh ua bands ab
		/*

		myuab = (uint32_t)uab;
		int cc = _popcnt32((uint32_t)uab);
		// add uab in field or in "more outfield
		uint32_t ua = myuab | _popcnt32(myuab) << 27;
		if (uab & sbb.mB.critbf) {
			if (nmoreif < 128)AddUA32(more_if, nmoreif, ua);
		}
		else if (nmoreof < 128)AddUA32(more_of, nmoreof, ua);

		g17b.moreb.Add((uint32_t)uab);// see final or multiple

		if (cc < 2)return 1;// should never be <2
		if (cc >= 4) {// add it to the reduced table in ua ub mode
			if (ntua < 512)		tua[ntua++] = myua;
			return 1;
		}
		int i36;
		for(int imini=0,bit=1,mask=7;imini<9;imini++,bit<<=1,mask<<=3){
			if(!(uab&mask) )continue;
			if (cc == 2) {// fresh mini2
				uint32_t bit27 = (uab&mask) ^ mask;
				bitscanforward(i36, bit27);
			}
			else i36 = imini+27;
			if (!ntuasmini[i36])activemini[nactivemini++] = i36;
			if (ntuasmini[i36] < 50)
				tuasmini[i36][ntuasmini[i36]++] =(uint32_t) uaa;
		}		*/

	
	}
	return (myua>0);
}

void BANDS_AB::CriticalFinalCheck_bbx(int bf) {// no more ua is it a valid solution


	p_cpt2g[8]++;
	register int ir = IsMultiple_bb(bf);
	if (!ir) {// one valid bands A + B
		p_cpt2g[9]++;
		// must relabel clues bands a b in bands 1 2 if necessary
		memcpy(tcluesb12, tclues, sizeof tcluesb12);
		if (ia) {// exchange banda bandb 
			for (int i = 0; i < nclues; i++) {
				uint32_t c = tclues[i];
				if (c < 27)tcluesb12[i] = c + 27;
				else tcluesb12[i] = c - 27;
			}
		}
		if (zhou[0].PartialInitSearch17(tcluesb12, nclues))
			return;// would be  bug
	}
}
void BANDS_AB::GoBand3() {
	if (g17b.debug17 ) cout << "go band3 ib3=" << cur_ib << endl;
	moreuas_b3.Init();
	moreuas_b3_small.Init();
	smin = genb12.bands3[cur_ib].smin;
	int nmiss = smin.minplus - smin.mincount;
	if (g17b.debug17 ) {
		cout << " nmiss=" << nmiss << " nuas 1=" << nuasb3_1 << " nuas 2=" << nuasb3_2 << endl;
		if (g17b.debug17 ) {
			cout << "table uab3_1" << endl;
			for (uint32_t i = 0; i < nuasb3_1; i++)
				cout << Char27out(uasb3_1[i]) << " i=" << i << endl;
			cout << "table uab3_2" << endl;
			for (uint32_t i = 0; i < nuasb3_2; i++)
				cout << Char27out(uasb3_2[i] ) << " i=" << i << endl;
		}
	}
	p_cpt2g[11 + nmiss] ++;
	memcpy(&genb12.grid0[54], genb12.bands3[cur_ib].band0, 4*27);
	G17B3HANDLER hh0; hh0.Init();
	cout << "status after init" << endl;
	hh0.smin.Status();
	hh0.diagh = 1;// sbb.diag;
	if(!nmiss)hh0.Go_Critical();
	else hh0.Go_Not_Critical_missn();

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

	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 =smin.critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = smin.minplus - smin.mincount;
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
	p_cpt2g[16]++;
	active_b3 = smin.critbf;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	if (!active_b3){
		bands_ab.FinalCheckB3(known_b3);
		return; // should be seen earlier if possible
	}
	cout<< Char27out(active_b3) <<" active after 2 pairs"<<endl;smin.Status();
	int ua = IsMultiple(known_b3 | active_b3);
	if (ua) {
		if (wua) *wua &= ua;
		return;
	}
	if(BuildIfShortB3())return;
	if (ShrinkUas1()) return;// dead branch
	cout << Char27out(active_b3) << " active aftershrinkuas" << endl; smin.Status();
	if (!active_b3) {
		cout << Char27out(known_b3) << " call final check" << endl;  
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
		PrintStatus();
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
void G17B3HANDLER::PrintStatus() {
	cout << "G17B3HANDLER Band3 Status" << endl;
	smin.Status();
}


//======================================================================= not critical sequence
void G17B3HANDLER::Go_Not_Critical_missn() {
	if (g17b.debug17 > 1 || diagh) 	cout<<Char27out(known_b3) << "entry not_critical miss " << nmiss
		<<" nuasb3of="<< nuasb3of << endl;
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
	else if (nmiss == 1 ) {
		andoutf = wactive0;
		for (uint32_t iua = 0; iua < nuasb3of; iua++)
			andoutf &= uasb3of[iua];
	}
	uint32_t wua = andoutf;
	if (nmiss > 1) wua = uasb3of[0] & BIT_SET_27;
	if (!nuasb3of)wua = wactive0;
	if (g17b.debug17 > 1 || diagh)		cout << Char27out(wua) << "wua to use  " << endl;

	if(wua){ // apply first UA to use or all out field cells 
		uint32_t res;
		while (bitscanforward(res, wua)) {
			if (g17b.debug17 > 1|| diagh)		cout << Char27out(wua) << "wua used nmiss= "<<nmiss << endl;
			int bit = 1 << res; wua ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this; hn.nmiss--; hn.known_b3 |= bit;
			if (hn.nmiss) {
				if (bands_ab.stack_filter) {
					int stack = C_stack[res];
					hn.stack_count.u16[stack]++;
					if (hn.stack_count.u16[stack] >= bands_ab.stack_filter) {
						hn.wactive0 &= ~(07007007 << (3 * stack));
					}
//					if (!hn.wactive0)continue;// bug
				}
				hn.Go_Not_Critical_missn();
			}
			else hn.Go_Critical(&wua);
			
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
	GodebugInit(0);
	if (GodebugCheckUas("check uas")) return 1;
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
	cout << "n bands3      \t" << genb12.nband3  << "\tua bands1+2   \t" << genuasb12.nua;
	cout<<"\t nuaothers "<<tusmall.n_others_b	<< endl;

	cout << "\nnguas socket2  \t" << genb12.ntua2 << endl;
	cout << "nguas socket3  \t" << genb12.ntua3<<endl;
	cout << "guas number sockets  2;3;all\t"
		<<tuguan.ng2<<"\t" << tuguan.ng3 << "\t" << tuguan.nguan << endl;

	if (mode & 1)tuguan.Debug1();// guas list and killers

	if (mode & 2) {
		cout << "table uas" << endl;
		uint64_t *t = genuasb12.tua;
		uint32_t n = genuasb12.nua;
		for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << endl;

	}
	if (mode & 4) {// status of guas2 3 first band
		STD_B3::GUAs & bguas = genb12.bands3[0].guas;
		for (uint32_t i = 0; i < 128; i++) {
			GUAN wg = tuguan.tguan[i];
			if (wg.ncol > 3) continue;
			int i81 = wg.i81;
			if (wg.ncol == 2 && bguas.isguasocket2.On_c(i81)) {// here 3X 27
				cout << Char27out(bguas.ua_pair[i81]) << " sock2 i=" << i
					<< " i81=" << i81 << " "
					<< Char2Xout(wg.killer) << "killer" << endl;
				uint64_t *tua = wg.tua;
				uint32_t nua = wg.nua;
				for (uint32_t iu = 0; iu < nua; iu++)
					cout << Char2Xout(tua[iu]) << endl;
			}

		}

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
			cout << "found iband1=" << iband1 << " iband2=" << iband2 << endl;
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