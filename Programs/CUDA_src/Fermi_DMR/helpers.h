void next_line(std::ifstream& scanner) {
	scanner.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void readNodes(std::string filename, FORD * &nodex, FORD * &nodey, unsigned &nnodes) {
	std::ifstream scanner(filename.append(".node").c_str());
	scanner >> nnodes;
	next_line(scanner);

	nodex = (FORD *) malloc(nnodes * sizeof(FORD));
	nodey = (FORD *) malloc(nnodes * sizeof(FORD));

	for (size_t i = 0; i < nnodes; i++) {
		size_t index;
		FORD x;
		FORD y;
		scanner >> index >> x >> y;
		next_line(scanner);
		nodex[index] = x;
		nodey[index] = y;
	}
}

void readTriangles(std::string basename, unsigned * &tnodes, unsigned &ntriangles, unsigned nnodes) {
	unsigned ntrianglesone, ntrianglestwo;
	unsigned i, index, n1, n2, n3, row;
	std::string filename;

	filename = basename;
	std::ifstream scanner(filename.append(".ele").c_str());
	scanner >> ntrianglesone;
	next_line(scanner);

	filename = basename;
	std::ifstream scannerperimeter(filename.append(".poly").c_str());
	next_line(scannerperimeter);
	scannerperimeter >> ntrianglestwo;
	next_line(scannerperimeter);

	ntriangles = ntrianglesone + ntrianglestwo;
	tnodes = (unsigned *) malloc(3 * ntriangles * sizeof(unsigned));

	cout << "ntrianglestwo : " << ntrianglestwo << endl;
	cout << "ntrianglesone : " << ntrianglesone << endl;
	for (i = 0; i < ntrianglesone; i++) {
		scanner >> index >> n1 >> n2 >> n3;
		row = 3 * index;
		tnodes[row + 0] = n1;
		tnodes[row + 1] = n2;
		tnodes[row + 2] = n3;
	}

	for (i = 0; i < ntrianglestwo; i++) {
		scannerperimeter >> index >> n1 >> n2 >> n3; //vivek n3 wasn't there before, it is must
		row = 3 * (ntrianglesone + index);
		tnodes[row + 0] = n1;
		tnodes[row + 1] = n2;
		tnodes[row + 2] = INVALIDID;

		//cout << "("<<tnodes[row+0] <<","<<tnodes[row+1] <<","<<tnodes[row+2] <<")"<<endl;
	}

}

void viv_readTriangles(std::string basename, unsigned * &tnodes, unsigned &ntriangles, unsigned nnodes, FORD * &nodex, FORD * &nodey) {
	unsigned ntrianglesone, ntrianglestwo;
	unsigned i, index, n1, n2, n3, row;
	std::string filename;

	filename = basename;
	std::ifstream scanner(filename.append(".ele").c_str());
	scanner >> ntrianglesone;
	next_line(scanner);

	filename = basename;
	std::ifstream scannerperimeter(filename.append(".poly").c_str());
	next_line(scannerperimeter);
	scannerperimeter >> ntrianglestwo;
	next_line(scannerperimeter);

	ntriangles = ntrianglesone + ntrianglestwo;
	tnodes = (unsigned *) malloc(3 * ntriangles * sizeof(unsigned));

	cout << "ntrianglestwo : " << ntrianglestwo << endl;
	cout << "ntrianglesone : " << ntrianglesone << endl;

#define LESSTHAN(NN1,NN2)  (nodex[NN1] < nodex[NN2]) || ( nodex[NN1] == nodex[NN2] && nodey[NN1] < nodey[NN2])
	for (i = 0; i < ntrianglesone; i++) {
		scanner >> index >> n1 >> n2 >> n3;
		row = 3 * index;

		tnodes[row + 0] = n1;
		tnodes[row + 1] = n2;
		tnodes[row + 2] = n3;

		if (LESSTHAN(n2,n1) || LESSTHAN(n3,n1)) {
			if (LESSTHAN(n2,n3)) {
				tnodes[row + 0] = n2;
				tnodes[row + 1] = n3;
				tnodes[row + 2] = n1;
			} else {
				tnodes[row + 0] = n3;
				tnodes[row + 1] = n1;
				tnodes[row + 2] = n2;
			}
		}
	}

	for (i = 0; i < ntrianglestwo; i++) {
		scannerperimeter >> index >> n1 >> n2 >> n3; //vivek n3 wasn't there before, it is must
		row = 3 * (ntrianglesone + index);
		tnodes[row + 0] = n1;
		tnodes[row + 1] = n2;
		tnodes[row + 2] = INVALIDID;

		if (LESSTHAN(n2,n1)) {
			tnodes[row + 0] = n2;
			tnodes[row + 1] = n1;
		}
		//cout << "("<<tnodes[row+0] <<","<<tnodes[row+1] <<","<<tnodes[row+2] <<")"<<endl;
	}

}

__device__
FORD distanceSquare(FORD onex, FORD oney, FORD twox, FORD twoy) {
	FORD dx = onex - twox;
	FORD dy = oney - twoy;
	FORD dsq = dx * dx + dy * dy;
	return dsq;
}
__device__
FORD distanceSquare(unsigned one, unsigned two, FORD *nodex, FORD *nodey) {
	return distanceSquare(nodex[one], nodey[one], nodex[two], nodey[two]);
}
__device__
FORD distance(unsigned one, unsigned two, FORD *nodex, FORD *nodey) {
	return sqrtf(distanceSquare(one, two, nodex, nodey));
}
__device__
FORD radiusSquare(FORD centerx, FORD centery, unsigned tri, FORD *nodex, FORD *nodey, unsigned *tnodes) {
	unsigned row = 3 * tri;
	unsigned first = tnodes[row + 0];
	return distanceSquare(centerx, centery, nodex[first], nodey[first]);
}
__device__
bool checkbad(unsigned id, FORD *nodex, FORD *nodey, unsigned *tnodes, DIMSTYPE *obtuse, unsigned ntriangles) {
	//if (id < ntriangles) {
	unsigned row = 3 * id;
	DIMSTYPE dims = (tnodes[row + 2] == INVALIDID ? 2 : 3);
	bool flag = false;
	//vivek changed because it is not guaranteed to find obtuse (not compatible with the javascript code)
	for (unsigned ii = 0; ii < dims; ++ii) {
		unsigned curr = tnodes[row + ii];
		unsigned aa = tnodes[row + (ii + 1) % dims];
		unsigned bb = tnodes[row + (ii + 2) % dims];
		if (curr < ntriangles && aa < ntriangles && bb < ntriangles) {
			FORD vax = nodex[aa] - nodex[curr];
			FORD vay = nodey[aa] - nodey[curr];
			FORD vbx = nodex[bb] - nodex[curr];
			FORD vby = nodey[bb] - nodey[curr];
			FORD dp = vax * vbx + vay * vby;

			if (dp < 0) {
				// id is obtuse at point ii.
				obtuse[id] = ii;
			} else {
				FORD dsqaacurr = distanceSquare(aa, curr, nodex, nodey);
				FORD dsqbbcurr = distanceSquare(bb, curr, nodex, nodey);
				FORD c = dp * rsqrtf(dsqaacurr * dsqbbcurr);
				if (c > cos(MINANGLE * (PI / 180))) {
					flag = true; //vivek
				}
			}
		}
	}
	//}
	return flag; //vivek
}
__global__
void dinit(FORD *nodex, FORD *nodey, unsigned *tnodes, bool *isbad, DIMSTYPE *obtuse, bool *isdel, unsigned nnodes, unsigned ntriangles) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < ntriangles) {
		obtuse[id] = 3;
		isbad[id] = checkbad(id, nodex, nodey, tnodes, obtuse, ntriangles);
		isdel[id] = false;
	}
}
__global__
void dverify(FORD *nodex, FORD *nodey, unsigned *tnodes, bool *isbad, bool *isdel, unsigned nnodes, unsigned ntriangles, bool *changed, unsigned *nchanged) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < ntriangles && !isdel[id] && isbad[id]) {
		*changed = true;
		++*nchanged;
	}
}
__device__
unsigned oadjacent(unsigned trione, unsigned tritwo, unsigned *tnodes, unsigned nnodes, unsigned ntriangles) {
	unsigned rowone = 3 * trione;
	unsigned rowtwo = 3 * tritwo;
	unsigned dimsone = (tnodes[rowone + 2] == INVALIDID ? 2 : 3);
	unsigned dimstwo = (tnodes[rowtwo + 2] == INVALIDID ? 2 : 3);
	unsigned ncommon = 0;
	unsigned firstmatch = 3; // not adjacent.

	for (unsigned ii = 0; ii < dimsone; ++ii) {
		for (unsigned jj = 0; jj < dimstwo; ++jj) {
			if (tnodes[rowone + ii] == tnodes[rowtwo + jj]) {
				if (++ncommon == 2) {
					return firstmatch;
				} else {
					firstmatch = ii;
				}
			}
		}
	}
	return 3; // not adjacent.
}

//vivek
__device__
unsigned adjacent(unsigned trione, unsigned tritwo, unsigned *tnodes, unsigned nnodes, unsigned ntriangles) {
	unsigned rowone = 3 * trione;
	unsigned rowtwo = 3 * tritwo;
	unsigned dimsone = (tnodes[rowone + 2] == INVALIDID ? 2 : 3);
	unsigned dimstwo = (tnodes[rowtwo + 2] == INVALIDID ? 2 : 3);

	for (unsigned ii = 0; ii < dimsone; ++ii) {
		unsigned pt1 = tnodes[rowone + ii];
		unsigned pt2 = tnodes[rowone + (ii + 1) % dimsone];
		for (unsigned jj = 0; jj < dimstwo; jj++) {
			unsigned t_pt1 = tnodes[rowtwo + jj];
			unsigned t_pt2 = tnodes[rowtwo + (jj + 1) % dimstwo];
			if ((pt1 == t_pt1 && pt2 == t_pt2) || (pt1 == t_pt2 && pt2 == t_pt1)) {
				return ii;
			}
		}
	}
	return 3;
}

__global__
void dfindneighbors(FORD *nodex, FORD *nodey, unsigned *tnodes, unsigned *neighbors, DIMSTYPE *neighboredges, unsigned nnodes, unsigned ntriangles, unsigned nblocks, unsigned starttri, unsigned endtri) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads=blockDim.x * gridDim.x;
	unsigned wpt = (endtri - starttri + nthreads - 1) / (nthreads);
	unsigned start = starttri + id * wpt;
	unsigned end = start + wpt;

	for (unsigned tt = start; tt < end && tt < ntriangles; ++tt) {
		unsigned row = 3 * tt;
		unsigned iirow = 0;
		for (unsigned ii = 0; ii < ntriangles; ++ii) {
			if (ii != tt) {
				unsigned commonedgestart = adjacent(tt, ii, tnodes, nnodes, ntriangles);
				if (commonedgestart < 3) { // common edge, adjacent.
					neighbors[row + iirow] = ii;
					neighboredges[row + iirow] = commonedgestart; // store the common edge for the first triangle, another thread will store it for the second triangle.
					++iirow;
					if(iirow>=3){
						break;
					}
				}
			}
		}
		// fill the remaining entries by invalid data.
		for (; iirow < 3; ++iirow) {
			neighbors[row + iirow] = INVALIDID;
			neighboredges[row + iirow] = 3;
		}
	}
}
__device__
unsigned getOpposite(unsigned centerelement, unsigned obtuse, unsigned *neighbors, DIMSTYPE *neighboredges, unsigned *tnodes, unsigned nnodes, unsigned ntriangles) {
	unsigned row = 3 * centerelement;
	DIMSTYPE dims = (tnodes[row + 2] == INVALIDID ? 2 : 3);
	DIMSTYPE commonedgepoint1 = (obtuse + 1) % dims;
	//unsigned commonedgepoint2 = (obtuse + 2) % dims;

	for (unsigned ii = 0; ii < 3; ++ii) { // iterate over neighbors.
		DIMSTYPE nnedgestart = neighboredges[row + ii];
		if (nnedgestart == commonedgepoint1) {
			return neighbors[row + ii];
		}
	}
	return INVALIDID;
}

__device__
void getCenter(unsigned centerelement, FORD &centerx, FORD &centery, FORD *nodex, FORD *nodey, unsigned *tnodes, unsigned nnodes, unsigned ntriangles) {
	unsigned row = 3 * centerelement;
	DIMSTYPE dims = (tnodes[row + 2] == INVALIDID ? 2 : 3);

	unsigned aa = tnodes[row + 0];
	unsigned bb = tnodes[row + 1];
	unsigned cc = tnodes[row + 2];


	if (dims == 3 && !(aa < ntriangles && bb < ntriangles && cc < ntriangles)) {
		centerx = centery = 0.0;
		printf("3 centerelement=%d ntriangles=%d cp\n", centerelement,ntriangles);
		return;
	}

	if (dims == 2 && !(aa < ntriangles && bb < ntriangles)) {
		centerx = centery = 0.0;
		printf("2 centerelement=%d ntriangles=%d cp\n", centerelement,ntriangles);
		return;
	}

	if (dims == 2) {
		centerx = (nodex[aa] + nodex[bb]) * 0.5;
		centery = (nodey[aa] + nodey[bb]) * 0.5;
		return;
	}
	FORD xxx = nodex[bb] - nodex[aa];
	FORD xxy = nodey[bb] - nodey[aa];
	FORD yyx = nodex[cc] - nodex[aa];
	FORD yyy = nodey[cc] - nodey[aa];

	FORD xxlen = distance(aa, bb, nodex, nodey);
	FORD yylen = distance(aa, cc, nodex, nodey);
	FORD cosine = (xxx * yyx + xxy * yyy) / (xxlen * yylen);
	FORD sinesq = 1.0 - cosine * cosine;
	FORD plen = yylen / xxlen;
	FORD ss = plen * cosine;
	FORD tt = plen * sinesq;
	FORD wp = (plen - cosine) / (2 * tt);
	FORD wb = 0.5 - (wp * ss);

	centerx = nodex[aa] * (1 - wb - wp) + nodex[bb] * wb + nodex[cc] * wp;
	centery = nodey[aa] * (1 - wb - wp) + nodey[bb] * wb + nodey[cc] * wp;
}
__device__
bool inCircumcircle(FORD xx, FORD yy, unsigned tri, FORD *nodex, FORD *nodey, unsigned *tnodes, unsigned nnodes, unsigned ntriangles) {
	// check if point (xx, yy) is in the circumcircle of tri.
	FORD centerx, centery;
	getCenter(tri, centerx, centery, nodex, nodey, tnodes, nnodes, ntriangles);
	FORD dd = distanceSquare(centerx, centery, xx, yy);
	return dd <= radiusSquare(centerx, centery, tri, nodex, nodey, tnodes);
}
__device__
unsigned addPoint(FORD xx, FORD yy, FORD *nodex, FORD *nodey, unsigned *pnnodes, unsigned &nnodes) {
	//unsigned newpoint = *pnnodes; ++*pnnodes;
	unsigned newpoint = atomicInc(pnnodes, MAXID); //vivek
	nodex[newpoint] = xx;
	nodey[newpoint] = yy;
	nnodes = newpoint; // update.
	return newpoint;
}
__device__
void addPoint(FORD xx, FORD yy, FORD *nodex, FORD *nodey, unsigned newpoint) {
	nodex[newpoint] = xx;
	nodey[newpoint] = yy;
}
__device__
void initNeighbors(unsigned tri, unsigned *neighbors, DIMSTYPE *neighboredges) {
	unsigned row = 3 * tri;
	for (unsigned ii = 0; ii < 3; ++ii) {
		neighbors[row + ii] = INVALIDID;
		neighboredges[row + ii] = 3;
	}
}
__device__
unsigned addTriangle(unsigned point0, unsigned point1, unsigned point2, unsigned *tnodes, unsigned *pntriangles, unsigned &ntriangles, bool *isdel, DIMSTYPE *obtuse, unsigned *neighbors, DIMSTYPE *neighboredges) {

	unsigned newtriid = atomicInc(pntriangles, MAXID);

	//printf("Adding tri:%d\n",newtriid);
	unsigned newrow = 3 * newtriid;
	tnodes[newrow + 0] = point0;
	tnodes[newrow + 1] = point1;
	tnodes[newrow + 2] = point2;
	initNeighbors(newtriid, neighbors, neighboredges);
	isdel[newtriid] = false;
	obtuse[newtriid] = 3;
	ntriangles = newtriid; // update.

	return newtriid;
}
__device__
void copyNeighbors(unsigned to, unsigned from, unsigned *neighbors, DIMSTYPE *neighboredges) {
	unsigned torow = 3 * to;
	unsigned fromrow = 3 * from;
	for (unsigned ii = 0; ii < 3; ++ii) {
		neighbors[torow + ii] = neighbors[fromrow + ii];
		neighboredges[torow + ii] = neighboredges[fromrow + ii]; // ???
	}
}
__device__
bool updateNeighbor(unsigned of, unsigned oldn, unsigned newn, unsigned *neighbors, unsigned *tnodes) {
	unsigned row = 3 * of;
	DIMSTYPE dims = (tnodes[row + 2] == INVALIDID ? 2 : 3);
	for (unsigned ii = 0; ii < dims; ++ii) {
		if (neighbors[row + ii] == oldn) {
			neighbors[row + ii] = newn;
			return true;
		}
	}
	// no need to update neighboredges, as the index won't change.
	return false;
}
void *mycudarealloc(void *oldptr, unsigned oldsize, unsigned newsize) {
	void *newptr;
	if (cudaMalloc((void **) &newptr, newsize) != cudaSuccess)
		CudaTest("allocating newptr failed");
	cudaMemcpy(newptr, oldptr, oldsize, cudaMemcpyDeviceToDevice);
	cudaFree(oldptr);
	return newptr;
}


//GPU lock-free synchronization function
__device__
void __gpu_sync(unsigned goalVal, volatile unsigned *Arrayin, volatile unsigned *Arrayout) {
	// thread ID in a block
	unsigned tid_in_blk = threadIdx.x * blockDim.y + threadIdx.y;
	unsigned nBlockNum = gridDim.x * gridDim.y;
	unsigned bid = blockIdx.x * gridDim.y + blockIdx.y;
	// only thread 0 is used for synchronization
	if (tid_in_blk == 0) {
		Arrayin[bid] = goalVal;
		__threadfence();
	}
	if (bid == 0) {
		if (tid_in_blk < nBlockNum) {
			while (Arrayin[tid_in_blk] != goalVal) {
				//Do nothing here
			}
		}
		__syncthreads();
		if (tid_in_blk < nBlockNum) {
			Arrayout[tid_in_blk] = goalVal;
			__threadfence();
		}
	}
	if (tid_in_blk == 0) {
		while (Arrayout[bid] != goalVal) {
			//Do nothing here
		}
	}
	__syncthreads();
}

__global__
void countbad(bool *isdel, bool *isbad, unsigned ntriangles, unsigned *nbad, unsigned goal, volatile unsigned *arrayin, volatile unsigned *arrayout, unsigned *blockcount) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned wpt = (ntriangles + nthreads - 1) / nthreads;
	unsigned start = id * wpt;
	unsigned end = start + wpt;
	__shared__ unsigned tcount[BLOCKSIZE];

	unsigned imemyself = threadIdx.x;

	/*
	 if(id==0)
	 printf("blockDim.x=%d\n",blockDim.x);
	 */

	tcount[imemyself] = 0;
	for (unsigned ii = start; ii < end; ++ii) {
		if (ii < ntriangles && isbad[ii] && (isdel[ii] == false)) {
			++tcount[imemyself];
		}
	}
	__syncthreads();

	for (unsigned s = blockDim.x / 2; s; s >>= 1) {
		if (imemyself < s) {
			tcount[imemyself] += tcount[imemyself + s];
		}
		__syncthreads();
	}
	__syncthreads();

	if (imemyself == 0) {
		blockcount[blockIdx.x] = tcount[0];
		__threadfence();
	}
	__gpu_sync(++goal, arrayin, arrayout);
	if (id == 0) {
		unsigned lcount = 0;
		for (unsigned ii = 0; ii < gridDim.x; ++ii) {
			lcount += blockcount[ii];
		}
		*nbad = lcount;
	}
}

__global__
void countnotdel(bool *isdel, unsigned ntriangles, unsigned *nbad, unsigned goal, volatile unsigned *arrayin, volatile unsigned *arrayout, unsigned *blockcount) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned wpt = (ntriangles + nthreads - 1) / nthreads;
	unsigned start = id * wpt;
	unsigned end = start + wpt;
	__shared__ unsigned tcount[BLOCKSIZE];

	unsigned imemyself = threadIdx.x;
	tcount[imemyself] = 0;
	for (unsigned ii = start; ii < end; ++ii) {
		if (ii < ntriangles && (isdel[ii] == false)) {
			++tcount[imemyself];
		}
	}
	__syncthreads();

	for (unsigned s = blockDim.x / 2; s; s >>= 1) {
		if (imemyself < s) {
			tcount[imemyself] += tcount[imemyself + s];
		}
		__syncthreads();
	}
	__syncthreads();

	if (imemyself == 0) {
		blockcount[blockIdx.x] = tcount[0];
		__threadfence();
	}
	__gpu_sync(++goal, arrayin, arrayout);
	if (id == 0) {
		unsigned lcount = 0;
		for (unsigned ii = 0; ii < gridDim.x; ++ii) {
			lcount += blockcount[ii];
		}
		*nbad = lcount;
	}
}

__device__
inline void CHECK_EXIST_ADD(unsigned int ary[], unsigned &len, unsigned val) {
	unsigned jj;
	for (jj = 0; jj < len; ++jj) {
		if (ary[jj] == val) {
			break;
		}
	}
	if (jj == len) {
		ary[len++] = val;
	}
}

__global__
void dynamic_cavity(unsigned iipost, unsigned *post, FORD *nodex, FORD *nodey, unsigned *tnodes, bool *isbad, DIMSTYPE *obtuse, bool *changed, unsigned ntriangles,
		unsigned *neighbors, DIMSTYPE *neighboredges, unsigned *pnnodes) {

}

__global__
void dynamicP_check_bad(unsigned npost, unsigned *post, FORD *nodex, FORD *nodey, unsigned *tnodes, bool *isbad, DIMSTYPE *obtuse, bool *changed, unsigned ntriangles) {

	unsigned nthread = blockDim.x * gridDim.x;
	unsigned ppt = (npost + nthread - 1) / nthread;
	unsigned lid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned start = ppt * lid;
	unsigned end = start + ppt;

	for (unsigned it = start; it < npost && it < end; it++) {
		unsigned tri = post[it];
		obtuse[tri] = 3;
		bool bad = checkbad(tri, nodex, nodey, tnodes, obtuse, ntriangles);
		isbad[tri] = bad;
		if (bad)
			*changed = bad;
	}

}

__global__
void dynamicP_check_bad1(unsigned iipost, unsigned *post, FORD *nodex, FORD *nodey, unsigned *tnodes, bool *isbad, DIMSTYPE *obtuse, bool *changed, unsigned ntriangles,
		unsigned *neighbors, DIMSTYPE *neighboredges, unsigned *pnnodes) {

	unsigned kk = threadIdx.x;
	unsigned nnodes = *pnnodes;

	if (kk < iipost) {
		unsigned tri = post[kk];

		unsigned iineighbor = 1;
		for (unsigned jj = 0; jj < iipost && iineighbor < 3; ++jj) {
			if (kk != jj) {
				DIMSTYPE commonedgestart = adjacent(tri, post[jj], tnodes, nnodes, ntriangles);
				if (commonedgestart < 3) {
					unsigned newrow = tri * 3;
					neighbors[newrow + iineighbor] = post[jj];
					neighboredges[newrow + iineighbor] = commonedgestart;
					++iineighbor;
				}
			}
		}

		obtuse[tri] = 3;
		bool bad = checkbad(tri, nodex, nodey, tnodes, obtuse, ntriangles);
		isbad[tri] = bad;
		if (bad)
			*changed = bad;
	}
}

#define 	PRINT_TRI(tri)				printf("	%d{%d,%d,%d} ([%d,%d,%d] , [ %d,%d,%d ])\n", \
												tri, myDmesh.tnodes[3*tri + 0], myDmesh.tnodes[3*tri + 1], \
												myDmesh.tnodes[3*tri + 2], myDmesh.neighboredges[3*tri + 0], \
												myDmesh.neighboredges[3*tri + 1], myDmesh.neighboredges[3*tri + 2], \
												myDmesh.neighbors[3*tri + 0], myDmesh.neighbors[3*tri + 1], \
												myDmesh.neighbors[3*tri + 2])

#define 	PRINT_TRI_ID(id,tri)				printf("id=%d	%d{%d,%d,%d} ([%d,%d,%d] , [ %d,%d,%d ])\n", \
												id,tri, tnodes[3*tri + 0], tnodes[3*tri + 1], \
												tnodes[3*tri + 2], neighboredges[3*tri + 0], \
												neighboredges[3*tri + 1], neighboredges[3*tri + 2], \
												neighbors[3*tri + 0], neighbors[3*tri + 1], \
												neighbors[3*tri + 2])

#define MAXITR	10
