//#define VALc15 13438
//#define VALc17 1032
//#define TESTXY 0
//2536435903110145  i1=2 i2=35
//#define TESTXY2 0
//22521159232789761
//536435903110145
#define LIM3Y 2000000
//#define DEBUGLEVEL 10  nb12=3461507
void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	const char * diagband = "285914637376528194914763528";
	const char * diagpuz = ".2.4........1....3........5.85......3...........76..2..6.....4......29..8....5...";
	diag = 0;
	//diag = 2; p17diag.SetAll_0();
	uint64_t diagval = 311976;
	cout << "entry m10 nb12=" << genb12.nb12 << " nbands3=" << genb12.nband3 << " p_cpt2g[0]=" << p_cpt2g[0] << endl;
	//if (p_cpt2g[0]) return;
	p_cpt2g[0] ++;
	p_cpt2g[1] +=genb12.nband3;
	if (genb12.nband3 > p_cpt2g[23])	p_cpt2g[23] = genb12.nband3;
	if (diag) {
		if (genb12.nb12 == diagval) {
			for (int i = 0; i < 54; i++)
				cout << genb12.grid0[i] + 1;
			cout << "entry m10 nb12=" << genb12.nb12 << endl;
			//if (strcmp(diagband, myband2.band)) return;
			cout << "this is the band in diag" << endl;
			if (diag == 2) {
				cout << diagpuz << " puz known" << endl;
				for (int i = 0; i < 81; i++) if (diagpuz[i] != '.')
					p17diag.Set_c(i);
			}
		}
		else return;
	}
	if (0) {
		for (int i = 0; i < 54; i++)
			cout << genb12.grid0[i] + 1;
		cout << "entry m10 "  << endl;
	}
	zh2b_g.test = 0;
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	myband2.DoExpandBand();// expand band2
	if (0) {
		myband2.DebugExpand();
		return;
	}
	cout << "myband1.n3_5=" << myband1.n3_5 << " myband2.n3_5=" << myband2.n3_5 << endl;
	if (!(myband1.n3_5| myband2.n3_5)) return; // no 656 no 566
	int nb3 = genb12.nband3;
	int ni3 = myband2.nind[2];

	//=========================== collect UAs  
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if(genuasb12.Initgen()) return;
	genb12.BuildGang9x3();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 
	if (0) GodebugInit(0);
	p_cpt2g[18] += genuasb12.nua;
	p_cpt2g[19] += genb12.ntua2;
	p_cpt2g[20] += genb12.ntua3;
	p_cpt2g[21] += genb12.nactive2;
	p_cpt2g[22] += genb12.nactive3;
	Go();// standard entry point for all 
}

void G17B::Go(){// start loops on 3_5 clues 
	if ((!sgo.vx[4]) && myband1.n3_5) {// A is band1
		bands_ab.Go(myband1, myband2, 0, 1);
	}
	if (myband2.n3_5) {
		// prepare zh2b_g with bands 2 first
		memcpy(zh2b_g.puz0, myband2.band0, sizeof myband2.band0);
		memcpy(&zh2b_g.puz0[27], myband1.band0, sizeof myband1.band0);
		zh2b_g.GetBands(myband2.gangster, myband1.gangster);// set sol/pm
		//_______________________________________________________
		bands_ab.Go(myband2, myband1, 1, 1);
	}
}

void BANDS_AB::Go(STD_B1_2 & ba, STD_B1_2 & bb, int i, int mode) {
	mybb = &bb;
	ni3 = ba.nxindex3;
	myt3_5 = ba.x_expand_3_5;
	myi3 = ba.xindex3;
	ia = i; ib = 1 - i;
	mode_ab = mode;// 1 if must be 5 clues 
	stack_filter = 6;
	// loop on index 3 
	cout << "go a b ni3="<<ni3<< endl;
	for (uint32_t i3 = 0; i3 < ni3; i3++) {
		wi3 = myi3[i3];
		Init3clues();
		// loop on remaining clues
		indf = myi3[i3 + 1].ideb;
		//cout << "go a b indf=" << indf << endl;
		for (uint32_t i5 = wi3.ideb; i5 < indf; i5++) {
			wi3_5 = myt3_5[i5];
			if (Init3_5clues()) sbb.Go();
		}
		//break;
	}
}
void BANDS_AB::Init3clues() {
	p_cpt2g[2]++;
	memcpy(tclues, wi3.tcells, sizeof wi3.tcells);
	uint64_t *t = genuasb12.tua;
	uint32_t n= genuasb12.nua;
	{ // reduce the UA table
		memset(ntuasmini, 0, sizeof ntuasmini);
		ntua = 0;
		register uint64_t F = wi3.cellsbf;
		if (ia) F <<= 32;// filter is on band 2
		//if (p_cpt2g[1] < 2) cout << Char2Xout(F) << " filter" << endl;
		for (uint32_t i = 0; i < n; i++) {
			register uint64_t U = t[i] & BIT_SET_2X, Ua, Ub;
			if (U&F)continue; // UA hit
			//if (p_cpt2g[1] < 2) cout << Char2Xout(U) << " kept i=" <<i<< endl;
			if (ia) { Ua = U >> 32; Ub = U & BIT_SET_27; }
			else { Ub = U >> 32; Ua = U & BIT_SET_27; }
			uint32_t countb = _popcnt32((uint32_t)Ub);
			if (countb > 3 ) {// store in in ua/ub mode
				Ub <<= 32;
				if(ntua<512)tua[ntua++] = Ua + Ub;
				continue;
			}
			// this is a mini row UA in band b use sub table
			for (int im = 0, mask = 7; im < 9; im++, mask <<= 3) {// find the mini row
				if (!(Ub&mask)) continue;
				if (countb == 3) im += 27;// now index for the subtable
				else {// must be count 2 and in mini row
					register uint32_t bit = (uint32_t)Ub ^ mask;
					bitscanforward(im,bit);
				}
				tuasmini[im][ntuasmini[im]++] = (uint32_t)Ua;
			}
		}
	}
	nactivemini = 0;
	for (int i = 0; i < 36; i++)if (ntuasmini[i])activemini[nactivemini++] = i;
	//if(p_cpt2g[1]<10)DebugInit();
}
void BANDS_AB::DebugInit() {
	uint32_t nred = ntua, nmax = 0, pairs = 0, triplets = 0;
	for (uint32_t i = 0; i < 36; i++) {
		register uint32_t j = ntuasmini[i];
		nred += j;
		if (nmax < j) nmax = j;
		if (j) {
			if (i < 27) pairs |= 1 << i;
			else triplets |= 1 << (j - 27);
		}
	}

	cout << " debuginit  ntua=" << ntua << " tot_red=" << nred<< " nactivemini="<< nactivemini << " max_mini=" << nmax << endl;
	cout << " trois premiers clues " << tclues[0] << " " << tclues[1] << " " << tclues[2] << endl;
	cout << Char27out(wi3.cellsbf) <<" filter ia="<<ia  << endl;
	cout << Char27out(pairs) << " pairs" << endl;
	cout << Char9out(triplets) << " triplets" << endl;
}
int BANDS_AB::Init3_5clues() {
	int locdiag = 0;
	g17b.moreb.Init();
	// build tclues band a (3 clues done)
	register uint32_t w = wi3_5.d;
	ncluesa = (w & 3) + 3;
	if (ncluesa != 5) return 0;
	p_cpt2g[3]++;
	if (ncluesa > 3) {
		for (int i = 3; i < ncluesa; i++) {
			w >>= 8;
			tclues[i] = w & 0xff;
		}
	}
	if (stack_filter) {// apply stack filter if valid
		stack_countba.u64 = 0;
		for (int i = 0; i < ncluesa; i++) {
			int cell = tclues[i], stack = C_stack[cell];
			stack_countba.u16[stack]++;
		}

	}


	filt32 = wi3.cellsbf;// build the filter filt32
	filt32 |= wi3_5.cellsbf;
	{ // find active minirows in band b
		register uint32_t F = filt32;
		sbb.pairs27 = sbb.critbf = sbb.triplet = 0;
		for (uint32_t i = 0; i < nactivemini; i++) {
			uint32_t iac = activemini[i], *t = tuasmini[iac], n = ntuasmini[iac];
			for (uint32_t iu = 0; iu < n; iu++) {
				if (!(t[iu] & F)) {// Ua not hit, must be in band b
					if (iac < 27) sbb.pairs27 |= 1 << iac;
					else sbb.triplet |= 1 << (iac - 27);
					break; //first unsolved is ok
				}
			}
		}
	}
	{ // set up mini rows status and field
		sbb.mini1 = sbb.mini2 = sbb.mini3 = 0;
		register uint32_t P = sbb.pairs27;
		for (int im = 0, mask = 7,bit=1; im < 9; im++, mask <<= 3,bit<<=1) {// check mini rows
			register uint32_t M = P & mask;
			if (!M) {
				if (sbb.triplet&bit)sbb.critbf |= mask;
			}
			else {
				sbb.triplet &= ~bit;// triplet would be redundant
				uint32_t cc = _popcnt32(M);
				if (cc > 1) {
					sbb.critbf |= mask;
					if(cc==2)sbb.mini2 |= bit;
					else sbb.mini3 |= bit;
				}
				else {
					sbb.critbf |= mask ^ M;
					sbb.mini1 |= bit;

				}
			}
		}
	}
	ncluesbandb = sbb.ncb2=6; //<<<<< this is for 656 ; 566 distribution 
	sbb.Ncrit();
	if(sbb.ncrit> ncluesbandb)return 0;
	sbb.tuaif = g17b.btuaif;
	sbb.tuaof = g17b.btuaof;

	{	// build remaining uas in and out field
		sbb.nuaif = sbb.nuaof = 0;
		sbb.andoutf = BIT_SET_27;
		register uint32_t F = filt32,IF=sbb.critbf;
		//cout << Char27out(IF) << " In Field " << endl;
		//cout << Char27out(F) << " Filter "   << endl;
		for (uint32_t i = 0; i < ntua; i++) {
			register uint64_t U = tua[i];
			//cout << Char2Xout(U) << " U "<< (U>>59) << endl;
			if ((uint32_t)U & F) continue;
			U >>= 32;// now ua bandb
			register uint32_t Ub = (uint32_t)U;
			if (Ub & IF) {
				//cout << Char27out(Ub) << " Ubif b12 " << endl;
				sbb.AddIF(Ub);
			}
			else {
				if (sbb.ncrit == ncluesbandb)return 0;
				//cout << Char27out(Ub) << " Ubof b12 "  << endl;
				sbb.AddOF(Ub);
				sbb.andoutf &= Ub;
				if ((!sbb.andoutf) && sbb.ncrit == (ncluesbandb - 1))return 0;
			}
		}
	}
	{	// add band b uas
		uint32_t *t=mybb->tua, n= mybb->nua;
		register uint32_t  IF = sbb.critbf;
		for (uint32_t i = 0; i < n; i++) {
			register uint32_t Ub = t[i];
			if (Ub & IF) {
				//cout << Char27out(Ub) << " Ubif bb " << endl;
				sbb.AddIF(Ub);
			}
			else {
				if (sbb.ncrit == ncluesbandb)return 0;
				//cout << Char27out(Ub) << " Ubof bb "  << endl;
				sbb.AddOF(Ub);
				sbb.andoutf &= Ub;
				if ((!sbb.andoutf) && sbb.ncrit == (ncluesbandb - 1))return 0;
			}
		}
	}
	if (locdiag)cout << "exit3_5 nuaif/of " << sbb.nuaif << " " << sbb.nuaof << endl;
	if (0) {
		cout << "table if" << endl;
		for (uint32_t i = 0; i < sbb.nuaif; i++)
			cout << Char27out(sbb.tuaif[i]) << endl;
		cout << "table of" << endl;
		for (uint32_t i = 0; i < sbb.nuaof; i++)
			cout << Char27out(sbb.tuaof[i]) << endl;
		sbb.Status();
	}
	{	//_______________ guas filter => already killed/forced plus sub table
		// build new subtables still active not fixed
		register uint32_t F = filt32;
		nguared_2 = 0;
		memset(ntuar2, 0, sizeof ntuar2);
		forced81_2.SetAll_0(); forced81_3.SetAll_0();
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			uint32_t * tuasw = tuar2[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t U = w.tua[iua] & BIT_SET_2X;
				register uint32_t Ua, Ub;
				if (ia) { Ua = U >> 32; Ub = U & BIT_SET_27; }
				else { Ub = U >> 32; Ua = U & BIT_SET_27; }
				if (Ua&F) continue; // hit by banda
				if (!Ub) {// empty in band 2 stop and force
					forced81_2.Set_c(i81);
					nt = 0;
					goto nexti81_2;
				}
				Ub |= (_popcnt32(Ub) << 27);
				AddUA32(tuasw, nt, Ub);
			}
			if (nt) {
				ntuar2[i81] = nt;
				guar2i81[nguared_2++] = i81;
			}
		nexti81_2:;
		}
		nguared_3 = 0;
		memset(ntuar3, 0, sizeof ntuar3);
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			uint32_t * tuasw = tuar3[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t U = w.tua[iua] & BIT_SET_2X;
				register uint32_t Ua, Ub;
				if (ia) { Ua = U >> 32; Ub = U & BIT_SET_27; }
				else { Ub = U >> 32; Ua = U & BIT_SET_27; }
				if (Ua&F) continue; // hit by banda
				if (!Ub) {// empty in band 2 stop and force
					forced81_3.Set_c(i81);
					nt = 0;
					goto nexti81_3;
				}
				Ub |= (_popcnt32(Ub) << 27);
				AddUA32(tuasw, nt, Ub);
			}
			if (nt) {
				ntuar3[i81] = nt;
				guar3i81[nguared_3++] = i81;
			}
		nexti81_3:;
		}
	}
	return 1;// can continue searching bandb solutions
}
void BANDS_AB::BANDB::Status() {
	cout << "Bandb Status ncrit=" <<ncrit<< endl;
	cout << Char27out(critbf) << " in field bf" << endl ;
	cout << Char27out(pairs27) << " pairs27" << endl ;
	cout << Char9out(mini_all) << " all minis2" << endl;
	cout << Char9out(mini1) << "     minis1" << endl;
	cout << Char9out(mini2) << "     minis2" << endl;
	cout << Char9out(mini3) << "     minis3" << endl ;
	cout << Char9out(triplet) << " mini triplets" << endl << endl;

}
//======================== start band b expansion

void BANDS_AB::BANDB::Go() {//start band b expansion
	//if(!ncrit) p_cpt2g[7] ++;
	//else 
	if (!nmiss) p_cpt2g[4] ++;
	else if (nmiss < 2)p_cpt2g[5] ++;
	else p_cpt2g[6] ++;
	//_______________
	g17b.SetUp(this);//catch ua tables location
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 = critbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	if (!nmiss)Go_Critical();
	else Go_Not_Critical_missn();
}

//_____________________ critical
void BANDS_AB::BANDB::CriticalAssignCell(int Ru) {// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & mini3) {// cell is in minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		mini3 ^= bit; // now only a pair to hit
		mini1 |= bit;
	}
	else {// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_b3 &= (~Mask); // kill the minirow as active
		mini1 &= ~bit;
		triplet &= ~bit;
	}
}
int BANDS_AB::BANDB::IsFinalMultiple() {
	// all if uas must be hit
	register int bf = known_b3 | active_b3;
	for (uint32_t iua = 0; iua < nuaif; iua++) {
		register int Ru = tuaif[iua];
		if(! (Ru & bf)) return 1;//know UA not hit multiple 
	}
	//check also nothing in more
	if (g17b.moreb.Check(bf)) return 1;
	if (!active_b3) {
		g17b.bands_ab.CriticalFinalCheck(known_b3);
		return 1;
	}
	if (bf != rknown_b3) {
		rknown_b3 = bf;
		if (g17b.bands_ab.IsMultiple(rknown_b3))
			return 1;
	}
	return 0;
}
void BANDS_AB::BANDB::Go_Critical() {// critical situation all clues in pairs tripl:ets
	//if (p_cpt2g[4] > 5) return;
	int locdiag = 0;
	if (locdiag)cout << Char27out(active_b3) << " active b3 entry critical" << endl;
	if (mini2) {// assign common cell 
		for (int i = 0, bit = 1, mask = 7; i < 9; i++ , bit <<= 1, mask <<= 7) {
			if (mini2&bit) {
				active_b3 &= (~mask); // clear the minirow
				known_b3 |= mask & (~pairs27);// and set the common cell a			
			}
		}
	}
	mini2 = 0;
	if (ShrinkUas1()) return;//assign possibles
	if (locdiag) {
		cout << Char27out(known_b3) << "known after shrink)" << endl;
		cout << Char27out(active_b3) << " active b3 after shrink" << endl;
		//return;
	}
	if (IsFinalMultiple())return;
	if (1) return;
	if (irloop)		CriticalLoop();
	else CriticalExitLoop();
}
int BANDS_AB::BANDB::ShrinkUas1() {
	irloop = 0;
	otuaif = tuaif;
	onuaif = nuaif;
	tuaif = &tuaif[nuaif];
	nuaif = 0;
	for (uint32_t iua = 0; iua < onuaif; iua++) {
		register int Ru = otuaif[iua];
		if (Ru & known_b3) continue;// already hit, forget it
		Ru &= active_b3;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else tuaif[nuaif++] = Ru;
	}
	if (!nuaif) irloop = 0;// no need to loop again
	return 0;

}
void BANDS_AB::BANDB::CriticalLoop() {
	if (ShrinkUas1()) return;
	if (irloop)CriticalLoop();
	else CriticalExitLoop();
}
void BANDS_AB::BANDB::CriticalExitLoop() {
	int nmissb = ncb2 - _popcnt32(known_b3);// missing clues
	if (nmissb < 0)return;
	if (!active_b3) {// nothing more to assign
		if (nuaif)return; // still not hit uas
		if (nmissb)return;// dead cell in a mini row 3 pairs
		g17b.bands_ab.CriticalFinalCheck(known_b3);
		return;
	}
	// check known + active with brute force
	int wknown = known_b3 | active_b3;
	if (rknown_b3 != wknown) {
		if (g17b.bands_ab.IsMultiple(wknown))	return;// not valid using all cells
		rknown_b3 = wknown;
	}
	if (nuaif) {		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_b3;
			for (uint32_t i = 0; i < nuaif; i++) {
				and_uas &= tuaif[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (mini1) {	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, mini1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3 & mask;// catch the minirow
		}
		else {
			for (uint32_t i = 0; i < nuaif; i++) {
				register int ua = tuaif[i], cc = _popcnt32(ua);
				if (cc < sizeua) { wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && triplet) {// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 & mask;// catch the minirow

			}
		}
		if (p_cpt2g[0] == sgo.vx[3])		cout << Char27out(wua) << " wua to use" << endl;

		while (bitscanforward(cell, wua)) {
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			BANDS_AB::BANDB hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop();
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void BANDS_AB::BANDB::Critical_0_UA() {
	int nmissb = ncb2 - _popcnt32(known_b3);// missing clues
	if (nmissb < 0)return;
	if (!nmissb) {// nothing more to assign (granted at first call in a branch)
		g17b.bands_ab.CriticalFinalCheck(known_b3);
		return;
	}
	if (mini3) {// in active minirows with 3 pairs, assign 2
		while (mini3) {
			uint32_t mini;
			bitscanforward(mini, mini3);
			int shift = 3 * mini, bit = 1 << shift;
			mini3 ^= 1 << mini; //clear bit the mini row is always killed
			active_b3 &= ~(7 << shift); // clear also the bitfield of active cells
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
	if (mini1) {// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, mini1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_b3 & mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		mini1 ^= 1 << mini;// and clear the minirow bit as active
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
	if (triplet) {// safety control should always be
		uint32_t mini;
		bitscanforward(mini, triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_b3 &= ~mask;// clear the minirow
		triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++) {
			int bb = bit << i;
			BANDS_AB::BANDB hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
// ____________________ sub critical
void  BANDS_AB::BANDB::Go_Subcritical() {
	active_b3 = active_sub = critbf;
	// check first if a global solution  is still possible (likely of no use here)
	if (IsFinalMultiple())return;// not valid using all cells
	uint32_t cct = _popcnt32(critbf) - ncrit;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	ndead = 0;
	Go_SubcriticalMiniRow();// find the first miss
}
void BANDS_AB::BANDB::Go_SubcriticalMiniRow() {
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1 << ndead, mask = 7 << (3 * ndead);
	for (int i = ndead; i < 9; i++, bit <<= 1, mask <<= 3) {
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & mini1) {// gua2 pair assign both
			BANDS_AB::BANDB hn = *this;
			hn.mini1 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & mini2) {// 2 gua2 pairs assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				BANDS_AB::BANDB hn = *this;
				hn.mini2 ^= bit;
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(M, mask);
			}
		}
		else if (bit & mini3) {// 3 gua2 pairs assign all
			BANDS_AB::BANDB hn = *this;
			hn.mini3 ^= bit;
			hn.SubMini(M, mask);
		}
		else if (bit & triplet) {// gua3 assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				BANDS_AB::BANDB hn = *this;
				hn.triplet ^= bit;
				hn.SubMini(M, mask);
			}
		}
		else { // second add in the mini row one residual cell take it
			BANDS_AB::BANDB hn = *this;
			hn.SubMini(M, mask);
		}
	}
}
void BANDS_AB::BANDB::SubMini(int M, int mask) {
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	if (nmiss) Go_SubcriticalMiniRow();// continue till no missing clue 
	else 		Go_Critical();	// leave sub critical   enter  critical
}
//_____________________ not critical

void BANDS_AB::BANDB::Go_Not_Critical_missn() {
	//	cout << Char27out(known_b3) << "entry not_critical miss " << nmiss << " nuaof="<< nuaof << endl;
	int wua = wactive0, ncells = 27, rawua = 0;
	{  // select ua to use
		register int Ra = wactive0, Rfilt = known_b3;
		for (uint32_t iua = 0; iua < nuaof; iua++) {
			register int Ru =tuaof[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			if (nmiss == 1) Ra = wua=Ru;
			else {
				register int cc = _popcnt32(Ru);
				if (cc < ncells) { ncells = cc; wua = Ru; rawua = tuaof[iua]; }
				if (cc > 6 && ncells < 5) break;
			}
		}
	}
	//cout << Char27out(wua) << "wua to use " << endl;
	if(wua){ // apply first UA to use or all out field cells
		if (nmiss == 1 && ncells < 27) {// time to do a global check
			uint32_t bf = critbf | known_b3 | wua;
			if (g17b.bands_ab.IsMultiple(bf)) 		return;
			rknown_b3 = bf; // to save redundant check multiple later
		}
		uint32_t res;
		int x = wua;
		while (bitscanforward(res, x)) {
			int bit = 1 << res; x ^= bit; wactive0 ^= bit;
			BANDS_AB::BANDB hn = *this; hn.nmiss--; hn.known_b3 |= bit;
			if (hn.nmiss) hn.Go_Not_Critical_missn();
			else 	hn.Go_Critical();
		}
	}
	if (ncells == 27 && ncrit)	Go_Subcritical();// finish in Subcritical if no ua
}


//_________________ 
int BANDS_AB::IsMultiple(int bf) {
	if (_popcnt32(bf) > 25) return 0;	
	nclues = ncluesa;// buil bandb part of tclues
	uint32_t cellb, wbf=bf;
	while( wbf){
		bitscanforward(cellb, wbf);
		wbf ^= 1 << cellb;
		tclues[nclues++] = cellb + 27;
	}
	int locdiag = 0;// (++p_cpt2g[25] < 20);
	if (locdiag) {
		cout<< nclues << " table clues\t";
		for (int i = 0; i < nclues; i++) cout << tclues[i] << " ";
		cout << endl;
		//return 1;
	}
	p_cpt2g[7]++;
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues, 0);
	//if (locdiag) cout <<Char2Xout(myua) << "uaret"<<endl;
	if ( myua) {//store the fresh ua bands ab
		if (locdiag) cout << Char2Xout(myua) << "uaret" << endl;
		register uint64_t uab = myua>>32,uaa= myua&BIT_SET_27,uab12=myua;
		if (ia)uab12 = uab | (uaa << 32);
		int cc = _popcnt32((uint32_t) uab),i36;
		genuasb12.AddUACheck(myua | ((uint64_t)cc << 59));// and update the UA table
		g17b.moreb.Add((uint32_t)uab);
		if(cc<2)return 1;// should never be <2
		if (cc >= 4) {// add it to the reduced table 
			if (ntua < 512)		tua[ntua++] = myua;
			return 1; 
		}
		for(int imini=0,bit=1,mask=7;imini<9;imini++,bit<<=1,mask<<=3){
			if(!(uab&mask) )continue;
			if (cc == 2) {// fresh mini2
				uint32_t bit27 = (uab&mask) ^ mask;
				bitscanforward(i36, bit27);
			}
			else i36 = imini+27;
			if (!ntuasmini[i36])activemini[nactivemini++] = i36;
			if (ntuasmini[i36] < 50)
				tuasmini[i36][ntuasmini[i36]++] =(uint32_t) uab;
		}
	}
	return (myua>0);
}
void BANDS_AB::CriticalFinalCheck(int bf) {// no more ua is it a valid solution
	p_cpt2g[8]++;
	register int ir = IsMultiple(bf);
	if (!ir) {// one valid bands A + B
		p_cpt2g[9]++;
		if (stack_filter) {// apply stack filter if valid
			stack_count = stack_countba;
			for (int i = ncluesa; i < nclues; i++) {
				int cell = tclues[i],stack= C_stack[cell];
				stack_count.u16[stack]++;
			}

		}
		final81_2 = forced81_2;// collect guas2 active
		for (uint32_t i = 0; i < nguared_2; i++) {
			int i81 = guar2i81[i];
			uint32_t * tua = tuar2[i81];
			for (uint32_t iua = 0; iua < ntuar2[i81]; iua++) {
				register uint32_t Ru = tua[iua];
				if (Ru&bf)continue;
				final81_2.Set_c(i81);
				break;// one ua not hit is enough here
			}
		}
		final81_3 = forced81_3;// collect guas3 active
		for (uint32_t i = 0; i < nguared_3; i++) {
			int i81 = guar3i81[i];
			uint32_t * tua = tuar3[i81];
			for (uint32_t iua = 0; iua < ntuar3[i81]; iua++) {
				register uint32_t Ru = tua[iua];
				if (Ru&bf)continue;
				final81_3.Set_c(i81);
				break;// one ua not hit is enough here
			}
		}
		int zinitdone = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			STD_B3::GUAs & myb = genb12.bands3[ib3].guas;
			BF128 ws2 = final81_2 & myb.isguasocket2;
			BF128 ws3 = final81_3 & myb.isguasocket3;
			// switch to mini rows patterns
			int tix[81], ntix = ws2.Table3X27(tix);
			mini_bf1 = mini_bf2 = mini_bf3 = pairsbf = pairs27 = mini_triplet = 0;
			for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
				int i81 = tix[i], imini = myb.ua2_imini[i81], bit = 1 << imini;
				if (mini_bf2&bit) mini_bf3 |= bit;
				if (mini_bf1&bit) mini_bf2 |= bit;
				mini_bf1 |= bit;
				pairsbf |= myb.ua_pair[i81];
				pairs27 |= 1 << myb.ua2_i27[i81];
			}
			ntix = ws3.Table3X27(tix);// now triplets to mini rows
			for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
				int imini = myb.ua3_imini[tix[i]], bit = 1 << imini;
				mini_triplet |= bit;
			}
			//___________________________ prepare a new band to process
			all_used_minis = mini_bf1 | mini_triplet;
			mini_triplet &= ~mini_bf1;// count only triplets with no pair
			mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
			if (mincount > ncluesb3) continue; // too many clues
			nmiss = ncluesb3 - mincount;
			mini_bf1 &= ~mini_bf2;// now pure one pair
			mini_bf2 &= ~mini_bf3;// now pure 2 pairs
			// set up pair + triplet bitfield
			if (mini_triplet) {// must add triplet minirow
				for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
					if (mini_triplet&bit)
						pairsbf |= field;
			}
			if (stack_filter) {//check stack limit
				stack_countf = stack_count;
				register uint32_t w = all_used_minis << 9| mini_bf3,
					stk=0111111,cc;
				for (int i = 0; i < 3; i++, stk <<= 1) {
					cc= stack_countf.u16[i]+ _popcnt32(w&stk);
					if (cc > 6) goto nextib3;
					stack_countf.u16[i] = (uint16_t)cc;
				}
			}
			{ //============= collect Gua46 and uas b3 for the band split them "in-field" "out-field"
				nuasb3_1 = nuasb3_2 = 0;
				register int  Rfilt = pairsbf;
				// first GUA46 usually shorter than UAs band3
				BF128  socket4 = genb12.bands3[ib3].guas.isguasocket4;// i81 3X
				socket4 &= final81_2;
				int * ua_46 = genb12.bands3[ib3].guas.ua_pair; // ua pattern
				int i81;
				while ((i81 = socket4.getFirstCell()) >= 0) {
					socket4.Clear_c(i81);// clear bit
					register uint32_t Ru = ua_46[i81] & BIT_SET_27;
					Ru |= _popcnt32(Ru) << 27;
					if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
					else if (!nmiss) goto nextib3;// critical + outfield uas
					else AddUA32(uasb3_2, nuasb3_2, Ru);
				}
				uint32_t * to = genb12.bands3[ib3].tua;
				for (uint32_t i = 0; i < genb12.bands3[ib3].nua; i++) {
					register uint32_t Ru = to[i] & BIT_SET_27;
					Ru |= _popcnt32(Ru) << 27;
					if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
					else if (!nmiss) goto nextib3;// critical + outfield uas
					else AddUA32(uasb3_2, nuasb3_2, Ru);
				}
				if ((!nmiss) && nuasb3_2) continue; // critical + outfield uas
			}
			//cout << " nmiss=" << nmiss <<" nuas 1="<< nuasb3_1 << " nuas 2=" << nuasb3_2 << endl;
			if (nmiss < 2)p_cpt2g[10 + nmiss] ++;
			else if (nmiss < 5)p_cpt2g[12] ++;
			else p_cpt2g[13] ++;
			p_cpt2g[15] ++;
			memcpy(&genb12.grid0[54], genb12.bands3[ib3].band0, 4 * 27);
			if (!zinitdone) {
				zinitdone = 1;
				if (zhou[0].PartialInitSearch17(tclues, nclues))
					return;// would be  bug
			}
			G17B3HANDLER hh0; hh0.Init(ib3);
			//if (!mincount) ExpandBand3();
			//else if (!nmiss)hh0.Go_Critical();
			//else hh0.Go_Not_Critical_missn();
		nextib3:;
		}
	
		if (0) {
			//cout << "one sol to print" << endl;
			char ws[82];
			strcpy(ws, empty_puzzle);
			for (int i = 0; i < nclues; i++) {
				int cell = tclues[i];
				ws[cell] = genb12.grid0[cell] + '1';
			}
			fout1 << ws << endl;
		}
	}

}


/*


void G17B3HANDLER::PrintStatus(int mode) {
	cout << "remaining uas table" << endl;
	if(mode&1)for (int i = 0; i < nuasb3; i++) 		cout << Char27out(uasb3[i]) << endl;
	cout << Char27out(active_b3) << " active_b3" << endl;
	//cout << Char27out(genb12.pairs27) << "genb12.pairs27" << endl;
	cout << Char9out(mini_bf1) << "mini bf1" << endl;
	cout << Char9out(mini_bf2) << "mini bf2" << endl;
	cout << Char9out(mini_bf3) << "mini bf3" << endl;
	cout << Char9out(mini_triplet) << "mini triplet" << endl;

}



*/


/*
void G17XY::Go_0(){// check stacks, check more uas, collect Guas status
	indexstep.diag_on = 0;
	p_cpt2g[3] ++;
	Init();
	if (g17b.debug17){
		if ((cellsbf&g17b.band1_17) != g17b.band1_17) return;
		if (((cellsbf >> 32) &g17b.band2_17) != g17b.band2_17) return;
		cout << Char2Xout(cellsbf) << " expected XY" << endl;
	}

	if (indexstep.diag_on) {
		indexstep.diag_on = 1;
		if ((cellsbf&g17b.band1_17) == g17b.band1_17) 
			if (((cellsbf >> 32) &g17b.band2_17) == g17b.band2_17) {
				cout << Char2Xout(cellsbf) << " expected XY" << endl;
				indexstep.diag_on = 2;
		}

	}
	if ( g17b.diag &4) {
		//if (cellsbf != g17b.p17diag.bf.u64[0]) return;
		if (cellsbf == g17b.p17diag.bf.u64[0]) {
			cout << Char2Xout(cellsbf) << " expected XY" << endl;
			indexstep.diag_on = 2;
		}
	}

	stacksxy.u64 = xxe.stacks.u64 + yye.stacks.u64;
	if (g17b.debug17) cout << "stacks " << stacksxy.u16[0] << stacksxy.u16[1] << stacksxy.u16[2] << endl;
	if (g17more.Check(cellsbf))return;		// check here more uas
	if (g17morebig.Check(cellsbf))return;
	Go_Guas_collect();	// collect guas still active
	if (g17b.debug17 > 1) cout << "call FirstCheckActiveBands()" << endl;
	if (indexstep.diag_on > 1)cout << "call FirstCheckActiveBands()" << endl;
	if (FirstCheckActiveBands()) return;
	p_cpt2g[4]++;
	if (g17b.debug17 > 1) cout << "call CheckValidBand12()" << endl;
	if (indexstep.diag_on > 1)cout << "call CheckValidBand12()" << endl;
	if (CheckValidBand12()) return;
	p_cpt2g[5]++;
	//if (TESTDEBUG)fout1 << "XY" << cellsbf << endl;
	//if (DEBUGLEVEL == 1) return;
	if (g17b.debug17) {
		cout << "entry Go_0 active nclues=" << nclues << " n=" << p_cptg[8]
			<< " stacks " << stacksxy.u16[0] << stacksxy.u16[1] << stacksxy.u16[2]
			<< "\t" << stacksxy.u16[3] << endl;
	}
	//genb12.Sockets2SetupForB12(cellsbf);// seems very expensive
	if (indexstep.diag_on > 1)cout << "call BuildActiveBands()" << endl;
	BuildActiveBands();// apply pairs triplets to bands 
	if (g17b.debug17)cout << "exit buidactivebands nb3=" << ntb3 << endl;
	if (indexstep.diag_on > 1)cout << "exit buidactivebands nb3=" << ntb3 << endl;
	if (!ntb3) return; // no active band b3lim and stack filters
	// at least one band is ok, time to check validity of the band 1+2
	p_cpt2g[8] ++;
	p_cpt2g[9] += ntb3;
	if (indexstep.diag_on > 1)cout << "go valide" << endl;
	Go_ValideB12();
}



int G17XY::CheckValidBand12(){
	// passing more filter now check validity or add a "more" ua
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues);
	//if (indexstep.diag_on > 1)
	if (myua){
		uint64_t cc = _popcnt64(myua);
		if (cc <= UALIMSIZE) {
			g17more.Add(myua);
			genuasb12.AddUACheck(myua | (cc << 59));// and update the UA table
		}
		else g17morebig.Add(myua);
		return 1; // not a solution unique for bands 1+2
	}
	return 0;
}

void G17XY::Go_ValideB12(){// UA2 and UA3 known not yet dead with min clues in band 3
	//cout <<Char2Xout(g17xy.cellsbf)<< "entry Go_ValideB12() ntb3=" << ntb3 << endl;

	int zinitdone = 0;
	for (int i3 = 0; i3 < ntb3; i3++){
		wg3 = g17tb3go[i3];
		int nmiss= 17 - wg3.countstack.u16[3];
		if(nmiss<2)p_cpt2g[10+nmiss] ++;
		else if (nmiss < 5)p_cpt2g[12] ++;
		else p_cpt2g[13] ++;
		if(wg3.minirows_triplets)p_cpt2g[14] ++;
		if (g17b.debug17>1) 	wg3.Debug();
		if (indexstep.diag_on > 1)wg3.Debug();


		if ((!nmiss) && nuasb3_2) continue; // critical + outfield uas
		p_cpt2g[15] ++;
		if (!zinitdone) {
			zinitdone = 1;
			if (zhou[0].PartialInitSearch17(tclues, nclues))
				return;// would be  bug 
		}
		memcpy(&genb12.grid0[54], genb12.bands3[wg3.ib3].band0, 4*27);
		//if (DEBUGLEVEL == 2) continue;
		//if (++p_cpt2g[19] > 10) return;
		//wg3.Debug();
		//for (int i = 54; i < 81; i++)
			//cout << zh_g2.grid0[i]+1;
		//cout << " band3 en zh_h2" << endl;


		g17hh0.Go();
//<<<<<<<<<<<<<<<<< a revoir
//		zhou[0].InitBand3PerDigit(genb12.bands3[wg3.ib3].band0);

		//if (g17b.debug17 > 1)		cout << " valide b12 go"  << endl;
	}
}

void G17XY::FoutValid17(int bf3, int ib3){
	char zs[82];
	int *g = genb12.grid0;
	strcpy(zs, empty_puzzle);
	for (int i = 0; i < nclues; i++) {
		int cell = tclues[i];
		zs[cell] = g[cell] + '1';
	}
	int bit = 1;
	g = genb12.bands3[ib3].band0;
	for (int i = 0; i < 27; i++, bit <<= 1)if (bf3&bit)
		zs[i + 54] = g[i] + '1';
	if (g17b.debug17) {
		fout1 << zs << ";" << g17b.npuz << endl;
		cout << zs << " valid found" << endl;
	}
	else fout1 << zs << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
}
*/
//================ part 2  band 3 procvessing
void ZHOU::InitBand3PerDigit(int * grid0b3) {
}
/*
void ZHOU::InitBand3PerDigit(int * grid0b3){
	memset(glb.band3digits, 0, sizeof glb.band3digits);
	register int * t = glb.band3digits;
	for (int i = 0; i < 27; i++){
		t[grid0b3[i]] |= 1 << i;
	}
	glb.digsols = zhxy27[0].glb.digsols;// catch pointer to solution per digit
}*/
int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{	
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = zh_g2.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) return 1;// can not be one solution
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
		if (p_cpt2g[0] == sgo.vx[3])cout << "exit mini aig=" << aig << endl;
		if (aig) {
			int ir = Apply17SingleOrEmptyCells();// restore zh_g
			if (p_cpt2g[0] == sgo.vx[3]) {
				ImageCandidats();
				cout << "retour apply17single " << ir << " zh_g.go_back= " << zh_g.go_back << endl;
				char ws[82];
				cout << zh_g.pairs.String3X(ws) << " paires" << endl;
				cout << zh_g2.triplets.String3X(ws) << " triplets" << endl;
			}
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
			zh_g.nsol++;
		}
		zh_g.go_back = 1;// closed anyway
		return;
	}
	Guess17(index , diag);// continue the process
}

void G17B3HANDLER::Init(int i) {
	BANDS_AB & bab = g17b.bands_ab;
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	ib3 = i;
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 = bab.pairsbf;// active cells in field
	pairs27 = bab.pairs27;
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = bab.nmiss;
	ncritical = bab.mincount;
	nb3 = nmiss + ncritical;
	mini_bf1 = bab.mini_bf1;
	mini_bf2 = bab.mini_bf2;
	mini_bf3 = bab.mini_bf3;
	mini_triplet = bab.mini_triplet;
	stack_count = bab.stack_countf;
	diagh = 0;
	if (bab.stack_filter) {// no active out in critical stacks
		for(int istack=0,stp=0111;istack<3;istack++,stp<<=1)
			if (stack_count.u16[istack] > 5) {// critical stack
				wactive0 &= ~(07007007 << (3 * istack));// clear outfield
				register int m2stack = stp & mini_bf2,shrink= TblShrinkMask[m2stack];
				if (m2stack) {// common cell(s) to assign
					register int Mask = tbitsrows[shrink] << (3 * istack);
					//adjust count and known
					known_b3 |= Mask & (~pairs27);// and set the common cell as assigned
					mini_bf2 &= ~stp; // clear the 2pairs bit(s) in stack
					active_b3 &= (~Mask);// clear the pairs in the in field bf
					pairs27 &= (~Mask);// clear the pairs in the pair bf
					ncritical -= _popcnt32(shrink);
				}
			}
	}
}

int G17B3HANDLER::IsMultiple(int bf) {
	if (_popcnt32(bf) > 25) return 0;
	BANDS_AB & bab = g17b.bands_ab;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, diagh);
	if (ir) {//consider store the fresh ua b3
		STD_B3 & myb3 = genb12.bands3[ib3];
		BF128 wua = zh_g2.cells_assigned;
		uint32_t ua = wua.bf.u32[2];
		int cc = _popcnt32(ua);
		if (cc < 4) {
			if (cc == 2) {// fresh gua2
				int i81 = myb3.GetI81_2(ua);
				if (i81 >= 0) {
					bab.final81_2.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar2[i81]) bab.guar2i81[bab.nguared_2++] = i81;;
					if (bab.ntuar2[i81] < GUAREDSIZE)
						bab.tuar2[i81][bab.ntuar2[i81]++] = wua.bf.u32[1];
					GEN_BANDES_12::SGUA2 & sg = genb12.tsgua2[i81];
					if (sg.nua >= SIZETGUA)sg.nua = SIZETGUA - 1;
					AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
			if (cc == 3) {// fresh gua2
				int i81 = myb3.GetI81_3(ua);
				if (i81 >= 0) {
					bab.final81_3.Set_c(i81);// new valid gua2 for other bands
					if (!bab.ntuar3[i81]) bab.guar3i81[bab.nguared_3++] = i81;;
					if (bab.ntuar3[i81] < GUAREDSIZE)
						bab.tuar3[i81][bab.ntuar3[i81]++] = wua.bf.u32[1];
					GEN_BANDES_12::SGUA3 & sg = genb12.tsgua3[i81];
					if (sg.nua >= SIZETGUA)sg.nua = SIZETGUA - 1;
					AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
			}
		}
	}
	return ir;
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~pairs27);// and set the common cell as assigned
			}
		}
		mini_bf2 = 0;
	}
}

void  G17B3HANDLER::Go(){
	if (nmiss){
		if (nmiss == g17b.b3lim){// no ua2_3 expand the uas table within the limits
			SPOT17_NOUAS2_3 s;
			memset(&s, 0, sizeof s);
			s.active = wactive0;
			s.stacks = wg3.countstack;
			s.stacks.u16[3] = nmiss;
			//s.newspot(g17xy.uasb3_2, g17xy.nuasb3_2, 0, 0);
		}
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
	if (bit & mini_bf3){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		mini_bf3 ^= bit; // now only a pairto hit
		mini_bf1 |= bit;
	}
	else{// either one pair or a triplet in the minirow
		active_b3 &= (~Mask); // kill the minirow as active
		mini_bf1 &= ~bit;
		mini_triplet &= ~bit;
	}
}
int G17B3HANDLER::ShrinkUas1(int * to, int no) {
	irloop = 0;
	uasb3 = &to[no];
	nuasb3 = 0;
	for (int iua = 0; iua < no; iua++) {
		register int Ru = to[iua];
		if (Ru & known_b3) continue;// already hit, forget it
		Ru &= active_b3;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else uasb3[nuasb3++] = Ru;
	}
	if (!nuasb3) irloop = 0;// no need to loop again
	return 0;

}
void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	p_cpt2g[16]++;
	active_b3 = wg3.pairs.u32[1];
	Critical2pairs();// assign 2 pairs in minirow to common cell
	if (!active_b3){
		if (diagh)cout << "final check1" << endl;
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (diagh) {
		cout << "critical status after assign 2 pairs" << endl;
		cout << Char27out(active_b3) << " active b3"<< endl;
		cout << Char27out(known_b3) << " known b3" << endl;
		wg3.Debug();
	}
	int bf = known_b3 | active_b3;
	if (bf != rknown_b3) {
		rknown_b3 = bf;
		if (IsMultiple(rknown_b3))	return;
	}

	if (g17b.debug17) {
		if ((rknown_b3 & g17b.band3_17) != g17b.band3_17) return;
		else 	cout << Char27out(rknown_b3) << " known b3 start critical" << endl;

	}
	if (diagh)cout << "back from is multiple)" << endl;
	//if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
	int wknown = known_b3 | active_b3;
	if (!active_b3) {
		if (diagh)cout << "final check2" << endl;
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (rknown_b3 != wknown)
	if (IsMultiple(wknown))	return;// not valid using all cells
	rknown_b3 = wknown;
	if (diagh || g17b.debug17)cout << "start  irloop="<<irloop << endl;
	if (irloop)		CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalLoop(int *to, int no){
	if (g17b.debug17>1)cout << "critical loop" << endl;
	if (ShrinkUas1(to, no)) return;
	if (irloop)CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalExitLoop(int *uasb3, int nuasb3){
	int nmissb = g17b.b3lim - _popcnt32(known_b3);// missing clues
	if (diagh || g17b.debug17>1 ){
		cout << Char27out(known_b3) << "known exit loop nuasb3= " << nuasb3 << endl;
		cout << Char27out(active_b3) << "active exit loop nmissb="<<nmissb << endl;
	}
	if (nmissb < 0)return;
	if (!active_b3){// nothing more to assign 
		if (nuasb3)return; // still not hit uas
		if (nmissb)return;// dead cell in a mini row 3 pairs
		CriticalFinalCheck();
		return;
	}
	// check known + active with brute force
	int wknown = known_b3 | active_b3;
	if (rknown_b3 != wknown) {
		if (IsMultiple(wknown))	return;// not valid using all cells
		rknown_b3 = wknown;
	}
	if (g17b.debug17>1) {
		cout<<Char9out(wg3.count.u16[0]) <<   "\tlook for next uab3 nuab3="<<nuasb3   << endl
			<<" priority to a smallest nuab3 if nmissb=1 else priority to pair if exists"<<endl;
	}

	if (nuasb3){		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case 
			register int and_uas = active_b3;
			for (int i = 0; i < nuasb3; i++) {
				and_uas &= uasb3[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (wg3.count.u16[0]){	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[0]);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3&mask;// catch the minirow
		}
		else{
			for (int i = 0; i < nuasb3; i++){
				register int ua = uasb3[i], cc = _popcnt32(ua);
				if (cc < sizeua){ wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum 
			}
			if (sizeua >= 2 && wg3.minirows_triplets){// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, wg3.minirows_triplets);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 &mask;// catch the minirow

			}
		}
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			if (g17b.debug17>1)cout << Char27out(wua) << " wua to assign in g17b.debug17 mode" << endl;
		}

		while (bitscanforward(cell, wua)){
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop(uasb3, nuasb3);
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void G17B3HANDLER::Critical_0_UA(){
	int nmissb = g17b.b3lim - _popcnt32(known_b3);// missing clues
	if (g17b.debug17>1 ){
		cout << Char27out(known_b3) << "known 0 ua" << endl;
		cout << Char27out(active_b3) << "active 0 ua " << endl;
	}
	if (nmissb < 0)return;
	if (!nmissb){// nothing more to assign (granted at first call in a branch)
		CriticalFinalCheck();
		return;
	}
	if (wg3.count.u16[2])	{// in active minirows with 3 pairs, assign 2
		while (wg3.count.u16[2]){
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[2]);
			int shift = 3 * mini, bit = 1 << shift;
			wg3.count.u16[2] ^= 1 << mini; //clear bit the mini row is always killed
			active_b3 &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++){
				int * tpi = tp[i];
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (wg3.count.u16[0]){// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, wg3.count.u16[0]);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		//wg3.count.u16[0] ^= 1 << mini; //clear bit the mini row is always killed
		int x = active_b3&mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		wg3.count.u16[0] ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++){
			int bb = bit << i;
			if (x&bb){
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (wg3.minirows_triplets){// safety control should always be
		uint32_t mini;
		bitscanforward(mini, wg3.minirows_triplets);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		//wg3.minirows_triplets ^= 1 << mini; //clear bit the mini row is always killed
		active_b3 &= ~mask;// clear the minirow
		wg3.minirows_triplets ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++){
			int bb = bit << i;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
void G17B3HANDLER::CriticalFinalCheck(){// no more ua is it a valid solution 
	//if (p_cpt2g[17] == 2252)
	int ncl = _popcnt32(known_b3);
	//if (g17b.debug17) cout << "final check test ncl=" << ncl  << endl;
	if (ncl != g17b.b3lim) return; // should be seen earlier if possible
	if (diagh)cout <<Char27out(known_b3 )<< "critical final check" << p_cpt2g[17] << endl;
	register int ir = IsMultiple(known_b3);// , g17b.debug17);
	if (diagh)cout << "retour" << ir << endl;
	if (ir){
		if (g17b.debug17&& nmiss != g17b.b3lim) {
			cout << "final check retour faux ir=" <<ir<< endl;
			//diagh = 1;
			//cout << "retour check debug=" << IsMultiple(known_b3) << endl;;
		}
		return;// mode debug pour voir
	}
	//g17xy.FoutValid17(known_b3, wg3.ib3);
	g17b.a_17_found_here = 1;
}
//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow(){
	if (g17b.debug17)cout << "entry Go_SubcriticalMiniRow() ndead="<<ndead << endl;
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	G17B3HANDLER hn;
	int bit = 1<<ndead, mask = 7<<(3*ndead), stack = ndead%3;
	for (int i = ndead; i < 9; i++, stack++, bit <<= 1, mask <<= 3){
		if (stack >2) stack = 0;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (g17b.debug17)cout << Char27out(M) << " mini row to process i="<<i << endl;
		if (bit & wg3.count.u16[0])// it was a gua2 pair assign both 
			hn.SubMini(this, M, mask, stack, 1, bit);
		else if (bit & wg3.count.u16[1])// it was 2 gua2 pair assign 2 out of 3 
			for (int j = 0; j < 3; j++){
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(this, M, mask, stack, 2, bit);
			}
		else if (bit & wg3.count.u16[2])// it was 3 gua2 pair assign 3 out of 3 
			hn.SubMini(this, M, mask, stack, 3, bit);
		else if (bit & wg3.minirows_triplets)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++){
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.wg3.minirows_triplets ^= bit;
				hn.SubMini(this, M, mask, stack, 4, bit);
			}
		else // second add in the mini row one residual cell take it
			hn.SubMini(this, M, mask, stack, 0, 0);
	}
}
void G17B3HANDLER::SubMini(G17B3HANDLER * o, int M, int mask, int stack, int imode, int bit){
	*this = *o;
	if (imode){
		if (imode<4)wg3.count.u16[imode-1] ^= bit;
		else wg3.minirows_triplets ^= bit;
	}
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added 
	active_b3 &= ~mask;
	active_sub ^= M;
	if (g17b.debug17) {
		if (g17b.debug17>1)cout << Char27out(known_b3) << " known sub mini  imode=" << imode << endl;
		if ((known_b3 & g17b.band3_17) != known_b3) {
			cout << "not the right one, goback" << endl;
			return;
		}
		if (g17b.debug17 > 1)cout << Char27out(active_b3) << " active_b3 nmiss=" << nmiss << endl;
	}
	// now adjust the stack count
	wg3.countstack.u16[stack]++;
	wg3.countstack.u16[3]++;
	if (wg3.countstack.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else{	// leave sub critical mode and enter the critical mode 
		if (g17b.debug17>1) wg3.Debug();
		Critical2pairs();// assign 2 pairs in minirow to common cell
		//if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
		rknown_b3 = known_b3 | active_b3;
		if (IsMultiple(rknown_b3))return;// not valid using all cells
		if (g17b.debug17>1)cout << "still valide" << endl;
		if (irloop)		CriticalLoop(uasb3, nuasb3);
		else CriticalExitLoop(uasb3, nuasb3);
	}
}
void G17B3HANDLER::Go_Subcritical(int docheck){// nmiss to select in the critical field
	p_cpt2g[16]++;
	active_b3 = active_sub = wg3.pairs.u32[1];
	// check first if a global solution  is still possible
	if(docheck)if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	int cct = _popcnt32(wg3.pairs.u32[1]) - ncritical;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	for (int ist = 0; ist < 3; ist++){// check stacks 
		if (wg3.countstack.u16[ist]>5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
//	uasb3 = g17xy.uasb3_1;
//	nuasb3 = g17xy.nuasb3_1;
	if (g17b.debug17)cout <<Char27out(active_sub )
		<< " active sub initial call to Go_SubcriticalMiniRow() nmiss="<<nmiss << endl;
	Go_SubcriticalMiniRow();// find the first miss
}
//======================================================================= not critical sequence
void G17B3HANDLER::Go_Not_Critical_missn() {
	if (0) 	cout << "entry not_critical miss " << nmiss << endl;
	
	int wua = wactive0, ncells = 27, rawua = 0;
	{  // select ua to use
		register int Ra = wactive0, Rfilt = known_b3;
		for (uint32_t iua = 0; iua < g17b.bands_ab.nuasb3_2; iua++) {
			register int Ru = g17b.bands_ab.uasb3_2[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			if (nmiss == 1) Ra = Ru;
			else {
				register int cc = _popcnt32(Ru);
				if (cc < ncells) { 
					ncells = cc; 
					wua = Ru; 
					rawua = g17b.bands_ab.uasb3_2[iua]; 
				}
				if (cc > 6 && ncells < 5) break;
			}
		}
	}
	if (0)		cout << Char27out(wua) << "wua to use " << endl;

	{ // apply first UA to use or all out field cells 
		uint32_t res;
		int x = wua;
		while (bitscanforward(res, x)) {
			int bit = 1 << res; x ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this; hn.nmiss--; hn.known_b3 |= bit;
			if (hn.nmiss) {
				if (g17b.bands_ab.stack_filter) {
					int stack = C_stack[res];
					hn.stack_count.u16[stack]++;
					if (hn.stack_count.u16[stack] >= g17b.bands_ab.stack_filter) {
						hn.wactive0 &= ~(07007007 << (3 * stack));
					}
					if (!hn.wactive0)continue;
				}
				hn.Go_Not_Critical_missn();
			}
			else 	hn.Go_Critical();
		}
	}
	if (ncells == 27)	Go_Subcritical();// finish in Subcritical if no ua
}




//==================================== handler direct 17 to test

void SPOT17_NOUAS2_3::newspot(int * oua, int onua, SPOT17_NOUAS2_3 * sp, int cell, int bit){
	if (sp){
		*this = *sp;// copy previous except if first
		// apply cell to stack count and check for the limit in stack
		known |= bit;
		int stack = C_box[cell];
		stacks.u16[stack]++;
		if (stacks.u16[stack] > 5)	active &= ~(07007007 << (3 * stack));
	}
	stacks.u16[3]--;// reduce the "missing cells" number
	// early check if 2 stacks filled 
	if (_popcnt32(active)<10 && g17b.g17hh0.IsMultiple(known | active)) return;
	tua = &oua[onua];
	nua = 0;
	int wua = active, ncells = 27;
	// if last step, must hit all remaining uas
	if (stacks.u16[3]){//================= first find a small ua
		register int Ra = active, Rfilt = known;
		for (int iua = 0; iua < onua; iua++){//and now UAs band3
			register int Ru = oua[iua];
			if (Ru & Rfilt) continue;// alreadyhit 
			Ru &= Ra;// no dead cells
			register int cc = _popcnt32(Ru);
			if (!cc)return; // dead branch
			if (cc < ncells)	{ ncells = cc; wua = Ru; }
			tua[nua++] = Ru;
		}
	}
	else{
		register int Ra = active, Rfilt = known;
		for (int iua = 0; iua < onua; iua++){//and now UAs band3
			register int Ru = oua[iua];
			if (Ru & Rfilt) continue;// alreadyhit 
			Ra &= Ru;// common cells
			if (!Ra)return; // can not hit all uas
		}
		if (g17b.g17hh0.IsMultiple(known | Ra)) return;// early check for last
		wua = Ra;
	}
	// apply the found wua
	uint32_t res;
	int x = wua;
	while (bitscanforward(res, x)){
		int bit = 1 << res;
		x ^= bit;// clear bit
		active ^= bit;// clear cell in that branch (dead cell)
		if (stacks.u16[3]){
			SPOT17_NOUAS2_3 sn;
			sn.newspot(tua, nua, this, res, bit);
		}
		else{// last step
			g17b.g17hh0.known_b3 = known | bit;
			g17b.g17hh0.CriticalFinalCheck();
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
void G17B::GodebugInit(int mode) {
	cout << "n bands3      \t" << genb12.nband3 << endl;
	cout << "ua bands1+2   \t" << genuasb12.nua << endl;
	cout << "guas socket2  \t" << genb12.ntua2 << endl;
	cout << "guas socket3  \t" << genb12.ntua3 << endl;
	cout << "active socket2\t" << genb12.nactive2 << endl;
	cout << "active socket3\t" << genb12.nactive3 << endl;
	if (mode & 1) {
		cout << "table uas" << endl;
		uint64_t *t = genuasb12.tua;
		uint32_t n = genuasb12.nua;
		for (uint32_t i = 0; i < n; i++) cout << Char2Xout(t[i]) << endl;

	}
	if (mode & 2) {
		cout << "sockets 2 table" << endl;
		int n2 = 0;
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			cout << i81 << " " << w.nua << endl;
			n2 += w.nua;
		}
		cout << "cumul=" << n2 << endl;
		cout << "sockets 3 table" << endl;
		int n3 = 0;
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			cout << i81 << " " << w.nua << endl;
			n3 += w.nua;
		}
		cout << "cumul=" << n3 << " total=" << n2 + n3 << endl;
	}

}
int G17B::GodebugFindKnown17() {// locate the known 17 

	return 0;
}
int G17B::GodebugCheckUas(const char * lib) {
	uint32_t nua = genuasb12.nua;
	uint64_t * tua= genuasb12.tua;
	for (uint32_t i = 0; i < nua; i++) {
		if (tua[i] & g17b.band12_17) continue;
		cout << lib << "check ua failed" << endl;
		cout << Char2Xout(tua[i]) << " not hit by the known 17" << endl;
		return 1;
	}
	return 0;
}

void V256_UAS::Debug(const char * lib, int mirror) {
	cout << " v256_ua for " << lib << endl;
	uint32_t * t32=v[0].bf.u32;
	for (int i = 0; i < 8; i++){
		int w = t32[i];
		if (mirror) w =~w;
		if(w)
		cout << Char32out(w) << " " << 32 * i << "-" << 32 * (i + 1) - 1 << endl;
	}
}
void V256_UAS::Fout(const char * lib) {
	fout1 << lib << v[0].bf.u64[0] << " " << v[0].bf.u64[1] << " "
		<< v[1].bf.u64[0] << " " << v[1].bf.u64[1] << endl;
}
void V256_UAS::Cout() {
	cout << Char64out(v[0].bf.u64[0]);
	cout<< endl;
}

void G17B::PrintEndPuz() {
	cout << endl;
	for (int i = 0; i < 10; i++) if (p_cpt1[i])
		cout << g17tl1[i] << "\t" << p_cpt1[i] << endl;
	cout << endl;
	for (int i = 0; i < 32; i++) if (p_cpt[i])
		cout << g17tl[i] << "\t" << p_cpt[i] << endl;
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