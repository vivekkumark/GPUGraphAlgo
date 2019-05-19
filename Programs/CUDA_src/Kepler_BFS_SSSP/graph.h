/* copied from 32.cu.
 * implements graph.
 */

#include <stdio.h>
#include <cuda.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <string.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <inttypes.h>
using namespace std;
#define foru unsigned
typedef struct Graph {
	enum {
		NotAllocated, AllocatedOnHost, AllocatedOnDevice
	} memory;

	Graph();
	//~Graph();
	unsigned cudaCopy(struct Graph &copygraph);
	unsigned readFromGR(char file[]);
	unsigned getSizeafterLoadBalance(unsigned &, int,int);
	void makeLoadBalance(int);

	void printfDegreeDistribution(char *);

	void findSerialSSSPnPrint(char *, unsigned);
	void findSerialBFSnPrint(char *, unsigned);

	void findSerialSSSP(unsigned);
	void findSerialBFS(unsigned);

	void NodeSplit(int,int);
	void printGraph(char *);
	unsigned getMaxAllowedOutDegree(int);
	unsigned getMaxdeg();
	unsigned getMaxdegIndex();
	unsigned init();
	unsigned allocOnHost();
	unsigned allocOnDevice();
	void testgraph();


	unsigned nnodes, nedges;
	unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst;foru *edgessrcwt;
	unsigned *parentNode;
	unsigned *childCount, *childIndex;
	unsigned *IsParent;

	foru *sssp, max_sssp;foru *bfs, max_bfs;

	unsigned nodes_b4_ns;

} Graph;

typedef struct _SResult {
	char exename[100];
	char algo[10];
	char graphfile[200];
	unsigned bucketSize;
	unsigned nnode;
	unsigned nedge;
	unsigned ns_nnode;
	unsigned ns_nedge;
	unsigned extranode;
	unsigned maxdeg;
	unsigned splitlevel;
	unsigned MAX_EDGES_ALLOWED;
	unsigned iterations;
	unsigned runtime;
	float time;
	float ktime;
	float sktime;
} SResult;
SResult myresult;

unsigned Graph::init() {
	noutgoing = nincoming = srcsrc = psrc = edgessrcdst = NULL;
	edgessrcwt = NULL;
	parentNode = NULL;
	nnodes = nedges = 0;
	memory = NotAllocated;
	return 0;
}

#define MYINFINITY	1000000000

void Graph::printfDegreeDistribution(char *fname){
	unsigned maxd=getMaxdeg()+1;
	unsigned *dd=(unsigned *)calloc(maxd,sizeof(unsigned));
	for(unsigned i=0;i<nnodes;i++){
		dd[noutgoing[i]]++;
	}

	ofstream op;
	op.open(fname);
	for (int i = 0; i < maxd; i++){
		if(dd[i])
			op << i << "," << dd[i] << endl;
	}
	op.close();

	free(dd);
}
void Graph::findSerialSSSP(unsigned src) {
	foru *localdist = (foru *) calloc(nnodes, sizeof(foru));
	for (int i = 0; i < nnodes; i++)
		localdist[i] = MYINFINITY;

	localdist[src] = 0;

	for (int i = 0; i < nnodes; i++) {
		int flag = 1;
		for (int ver = 0; ver < nnodes; ver++) {
			foru tdist = localdist[ver];
			for (int ed = psrc[ver]; ed < (psrc[ver] + noutgoing[ver]); ed++) {
				foru wt = edgessrcwt[ed];
				if (tdist + wt < localdist[edgessrcdst[ed]]) {
					flag = 0;
					localdist[edgessrcdst[ed]] = tdist + wt;
				}
			}
		}
		if (flag)
			break;
	}

	max_sssp = 0;
	for (int ver = 0; ver < nnodes; ver++) {
		if (localdist[ver] < MYINFINITY) {
			if (localdist[ver] > max_sssp) {
				max_sssp = localdist[ver];
			}
		}
	}

	sssp = localdist;
}

void Graph::findSerialBFS(unsigned src) {

	foru *localdist = (foru *) calloc(nnodes, sizeof(foru));
	for (int i = 0; i < nnodes; i++)
		localdist[i] = MYINFINITY;

	localdist[src] = 0;

	for (int i = 0; i < nnodes; i++) {
		int flag = 1;
		for (int ver = 0; ver < nnodes; ver++) {
			foru tdist = localdist[ver];
			for (int ed = psrc[ver]; ed < (psrc[ver] + noutgoing[ver]); ed++) {
				foru wt = 1;
				if (tdist + wt < localdist[edgessrcdst[ed]]) {
					flag = 0;
					localdist[edgessrcdst[ed]] = tdist + wt;
				}
			}
		}
		if (flag)
			break;
	}

	max_bfs = 0;
	for (int ver = 0; ver < nnodes; ver++) {
		if (localdist[ver] < MYINFINITY) {
			if (localdist[ver] > max_bfs) {
				max_bfs = localdist[ver];
			}
		}
	}

	bfs = localdist;
}

void Graph::findSerialSSSPnPrint(char *fname, unsigned src) {

	foru *localdist = (foru *) calloc(nnodes, sizeof(foru));
	for (int i = 0; i < nnodes; i++)
		localdist[i] = MYINFINITY;

	localdist[src] = 0;

	for (int i = 0; i < nnodes; i++) {
		int flag = 1;void block_aware_LB(int);
		for (int ver = 0; ver < nnodes; ver++) {
			foru tdist = localdist[ver];
			for (int ed = psrc[ver]; ed < (psrc[ver] + noutgoing[ver]); ed++) {
				foru wt = edgessrcwt[ed];
				if (tdist + wt < localdist[edgessrcdst[ed]]) {
					flag = 0;
					localdist[edgessrcdst[ed]] = tdist + wt;
				}
			}
		}
		if (flag)
			break;
	}
	ofstream op;
	op.open(fname);
	op.precision(0);
	for (int i = 0; i < nnodes; i++)
		op << fixed << localdist[i] << endl;
	op.close();

	free(localdist);
}
void Graph::findSerialBFSnPrint(char *fname, unsigned src) {

	foru *localdist = (foru *) calloc(nnodes, sizeof(foru));
	for (int i = 0; i < nnodes; i++)
		localdist[i] = MYINFINITY;

	localdist[src] = 0;

	for (int i = 0; i < nnodes; i++) {
		int flag = 1;
		for (int ver = 0; ver < nnodes; ver++) {
			foru tdist = localdist[ver];
			for (int ed = psrc[ver]; ed < (psrc[ver] + noutgoing[ver]); ed++) {
				foru wt = 1;
				if (tdist + wt < localdist[edgessrcdst[ed]]) {
					flag = 0;
					localdist[edgessrcdst[ed]] = tdist + wt;
				}
			}
		}
		if (flag)
			break;
	}
	ofstream op;
	op.open(fname);
	op.precision(0);
	for (int i = 0; i < nnodes; i++)
		op << fixed << localdist[i] << endl;
	op.close();

	free(localdist);
}
unsigned Graph::allocOnHost() {
	edgessrcdst = (unsigned int *) malloc((nedges + 1) * sizeof(unsigned int)); // first entry acts as null.
	edgessrcwt = (foru *) malloc((nedges + 1) * sizeof(foru)); // first entry acts as null.
	psrc = (unsigned int *) calloc(nnodes + 1, sizeof(unsigned int)); // init to null.
	psrc[nnodes] = nedges; // last entry points to end of edges, to avoid thread divergence in drelax.
	noutgoing = (unsigned int *) calloc(nnodes, sizeof(unsigned int)); // init to 0.
	nincoming = (unsigned int *) calloc(nnodes, sizeof(unsigned int)); // init to 0.
	srcsrc = (unsigned int *) malloc(nnodes * sizeof(unsigned int));
	memory = AllocatedOnHost;
	IsParent = (unsigned *) calloc(nnodes, sizeof(unsigned));
	parentNode = (unsigned int *) malloc(nnodes * sizeof(unsigned int));
	childCount = (unsigned *) calloc(nnodes, sizeof(unsigned));
	childIndex = (unsigned *) calloc(nnodes, sizeof(unsigned));
	return 0;
}

void writeGraph(char *fname, foru *ary, unsigned size) {
	FILE *fp;
	fp = fopen(fname, "w+t");
	for (int i = 0; i < size; i++) {
		fprintf(fp, "%f\n", ary[i]);
	}
	fclose(fp);
}


typedef struct _comp {
	unsigned nout;
	unsigned src;
	unsigned n;
	unsigned cc;
	unsigned par;
} SCOMP;

int comp(const void * elem1, const void * elem2) {
	unsigned f = *((unsigned *) elem1);
	unsigned s = *((unsigned *) elem2);
	if (f > s)
		return -1;
	if (f < s)
		return 1;
	return 0;
}
int comp1(const void * elem1, const void * elem2) {
	SCOMP f = *((SCOMP *) elem1);
	SCOMP s = *((SCOMP *) elem2);
	if (f.nout > s.nout)
		return -1;
	if (f.nout < s.nout)
		return 1;
	return 0;
}
void Graph::printGraph(char *fname) {
	//writeGraph("edgewt.txt",edgessrcwt,nedges);
	//exit(0);

	SCOMP *scomp = (SCOMP *) malloc(sizeof(SCOMP) * nnodes);

	for (int i = 0; i < nnodes; i++) {
		scomp[i].nout = noutgoing[i];
		scomp[i].src = psrc[i];
		scomp[i].n = i;
		scomp[i].cc = childCount[i];
		scomp[i].par = IsParent[i];
	}

	//qsort(scomp,nnodes,sizeof(SCOMP),comp1);

	//unsigned *Lnoutgoing;
	//Lnoutgoing=(unsigned *)malloc(sizeof(unsigned)*nnodes);
	//memcpy(Lnoutgoing,noutgoing,sizeof(unsigned)*nnodes);

	FILE *fp;
	fp = fopen(fname, "w+t");
	for (int i = 0; i < nnodes; i++) {
		if (scomp[i].cc || scomp[i].par == 0) {
			fprintf(fp, "%10u %5u :", scomp[i].n, scomp[i].nout);
			unsigned index = scomp[i].src;
			for (int j = 0; j < scomp[i].nout; j++) {
				fprintf(fp, ",%u", edgessrcdst[index + j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	free(scomp);
	//free(Lnoutgoing);
}
unsigned Graph::getMaxdeg() {
	unsigned maxdeg = noutgoing[0];
	for (int j = 1; j < nnodes; j++) {
		if (maxdeg < noutgoing[j]) {
			maxdeg = noutgoing[j];
		}
	}
	myresult.maxdeg = maxdeg;
	return maxdeg;
}
unsigned Graph::getMaxdegIndex() {
	unsigned maxdeg = noutgoing[0];
	unsigned index = 0;
	for (int j = 1; j < nnodes; j++) {
		if (maxdeg < noutgoing[j]) {
			maxdeg = noutgoing[j];
			index = j;
		}
	}
	return index;
}
unsigned Graph::getMaxAllowedOutDegree(int maxHisto) {
	unsigned *histo = (unsigned *) calloc(maxHisto, sizeof(unsigned));
	int maxdeg = getMaxdeg();
	int maxdeg1 = maxdeg + 1;

	for (int i = 0; i < nnodes; i++) {
		int index = (noutgoing[i] / (float) maxdeg1) * maxHisto;
		histo[index]++;
	}

	printf("Histo\n");
	for (int i = 0; i < maxHisto; i++) {
		printf("%d , ", histo[i]);
	}
	printf("nnodes=%d\n", nnodes);

	int max = -1, maxindex = -1;
	for (int i = 0; i < maxHisto; i++) {
		if ((int) histo[i] > max) {
			max = histo[i];
			maxindex = i;
		}
	}
	maxindex++;
	int maxOutDeg = (maxdeg1 / maxHisto) * maxindex;
	return maxOutDeg;
}
unsigned Graph::getSizeafterLoadBalance(unsigned &MAX_EDGES_ALLOWED, int NodeSplitLevel,int maxHisto) {

	unsigned maxdeg = getMaxdeg();

	MAX_EDGES_ALLOWED = (maxdeg + NodeSplitLevel - 1) / NodeSplitLevel;
	if (NodeSplitLevel < 0) {
		MAX_EDGES_ALLOWED = getMaxAllowedOutDegree(maxHisto);
		printf("Approx Level=%f\n", (float) maxdeg / MAX_EDGES_ALLOWED);
	}
	printf("maxOutDegree=%d\n", MAX_EDGES_ALLOWED);

	myresult.MAX_EDGES_ALLOWED = MAX_EDGES_ALLOWED;

	unsigned i, MY_VERTICE_INDEX = 0;

	for (i = 0; i < nnodes; i++) {
		unsigned nout = noutgoing[i];
		if (nout == 0) {
			MY_VERTICE_INDEX++;
			continue;
		}
		int SPLIT;
		SPLIT = (nout + MAX_EDGES_ALLOWED - 1) / MAX_EDGES_ALLOWED;
		MY_VERTICE_INDEX += SPLIT;
	}

	return MY_VERTICE_INDEX;

}
void Graph::NodeSplit(int NodeSplitLevel,int maxHisto) {
	unsigned MAX_EDGES_ALLOWED;
	unsigned getLoadSize = getSizeafterLoadBalance(MAX_EDGES_ALLOWED, NodeSplitLevel,maxHisto);
	//MAX_EDGES_ALLOWED=80;
	unsigned int *MY_hpsrc = (unsigned int *) calloc(getLoadSize + 1, sizeof(unsigned int));
	unsigned int *MY_hsrcdist = (unsigned int *) calloc(getLoadSize, sizeof(unsigned int)); // init to null.
	unsigned int *MY_hsrcsrc = (unsigned int *) malloc(getLoadSize * sizeof(unsigned int));
	unsigned int *MY_hnoutgoing = (unsigned int *) malloc(getLoadSize * sizeof(unsigned int));
	unsigned int *MY_hnincoming = (unsigned int *) malloc(getLoadSize * sizeof(unsigned int));

	unsigned int *MY_IsParent = (unsigned *) calloc(getLoadSize, sizeof(unsigned));
	unsigned int *MY_childCount = (unsigned *) calloc(getLoadSize, sizeof(unsigned));

	unsigned int i, MY_VERTICE_INDEX = 0;
	unsigned int *MY_presum = (unsigned int *) malloc((nnodes + 1) * sizeof(unsigned int));
	MY_presum[0] = 0;

	for (i = 0; i < nnodes; i++) {
		if (noutgoing[i] == 0) {
			MY_presum[i + 1] = MY_presum[i];
			MY_hnoutgoing[MY_VERTICE_INDEX] = 0;

			if (MY_VERTICE_INDEX == 0) {
				MY_hpsrc[MY_VERTICE_INDEX] = 1;
			} else {
				MY_hpsrc[MY_VERTICE_INDEX] = MY_hpsrc[MY_VERTICE_INDEX - 1] + MY_hnoutgoing[MY_VERTICE_INDEX - 1];
			}

			MY_hsrcdist[MY_VERTICE_INDEX] = i;
			MY_hsrcsrc[MY_VERTICE_INDEX] = MY_VERTICE_INDEX;

			MY_IsParent[MY_VERTICE_INDEX] = 1;
			MY_VERTICE_INDEX++;
			continue;
		}

		int SPLIT = (noutgoing[i] + MAX_EDGES_ALLOWED - 1) / MAX_EDGES_ALLOWED;
		if (SPLIT > 1) {
			MY_presum[i + 1] = MY_presum[i] + SPLIT - 1;
			MY_childCount[MY_VERTICE_INDEX] = SPLIT - 1;
		} else
			MY_presum[i + 1] = MY_presum[i];

		MY_IsParent[MY_VERTICE_INDEX] = 1;
		int count = 0;

#ifdef INT_SPLIT
		int local_max_deg=(noutgoing[i]+SPLIT-1)/SPLIT;
		int so_for=0;
#else
		int rc = noutgoing[i];
#endif

		while (SPLIT--) {
#ifdef INT_SPLIT
			if(SPLIT==0) {
				count = noutgoing[i] - so_for;
			}
			else {
				count = local_max_deg;
				so_for += count;
			}
#else
			if (rc > MAX_EDGES_ALLOWED) {
				count = MAX_EDGES_ALLOWED;
				rc = rc - MAX_EDGES_ALLOWED;
			} else {
				count = rc;
				rc = 0;
			}
#endif
			MY_hnoutgoing[MY_VERTICE_INDEX] = count;
			if (MY_VERTICE_INDEX == 0) {
				MY_hpsrc[MY_VERTICE_INDEX] = 1;
			} else {
				MY_hpsrc[MY_VERTICE_INDEX] = MY_hpsrc[MY_VERTICE_INDEX - 1] + MY_hnoutgoing[MY_VERTICE_INDEX - 1];
			}
			MY_hsrcdist[MY_VERTICE_INDEX] = i;
			MY_hsrcsrc[MY_VERTICE_INDEX] = MY_VERTICE_INDEX;
			MY_VERTICE_INDEX++;
		}
	}

	MY_hpsrc[MY_VERTICE_INDEX] = nedges;

	for (i = 1; i < (nedges + 1); i++) {
		edgessrcdst[i] = edgessrcdst[i] + MY_presum[edgessrcdst[i]];
	}

	free(srcsrc);
	free(psrc);
	free(noutgoing);
	free(nincoming);
	free(IsParent);
	free(childCount);

	nodes_b4_ns = nnodes;
	nnodes = MY_VERTICE_INDEX;

	srcsrc = MY_hsrcsrc;
	psrc = MY_hpsrc;
	noutgoing = MY_hnoutgoing;
	nincoming = MY_hnincoming;
	IsParent = MY_IsParent;
	parentNode = MY_hsrcdist;
	childCount = MY_childCount;

}
static unsigned CudaTest(char *msg) {
	cudaError_t e;
	cudaThreadSynchronize();
	if (cudaSuccess != (e = cudaGetLastError())) {
		fprintf(stderr, "%s: %d\n", msg, e);
		fprintf(stderr, "%s\n", cudaGetErrorString(e));
		exit(-1);
		//return 1;
	}
	return 0;
}
unsigned Graph::allocOnDevice() {
	if (cudaMalloc((void **) &edgessrcdst, (nedges + 1) * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating edgessrcdst failed");
	if (cudaMalloc((void **) &edgessrcwt, (nedges + 1) * sizeof(foru)) != cudaSuccess)
		CudaTest("allocating edgessrcwt failed");
	if (cudaMalloc((void **) &psrc, (nnodes + 1) * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating psrc failed");
	if (cudaMalloc((void **) &noutgoing, nnodes * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating noutgoing failed");
	if (cudaMalloc((void **) &nincoming, nnodes * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating nincoming failed");
	if (cudaMalloc((void **) &srcsrc, nnodes * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating srcsrc failed");
	if (cudaMalloc((void **) &parentNode, nnodes * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating parentNode failed");
	if (cudaMalloc((void **) &childCount, nnodes * sizeof(unsigned int)) != cudaSuccess)
		CudaTest("allocating childCount failed");
	memory = AllocatedOnDevice;
	return 0;
}
Graph::Graph() {
	init();
}
void Graph::testgraph() {
	nnodes = 11;
	nedges = 10;

	noutgoing[0] = 3;
	noutgoing[1] = 0;
	noutgoing[2] = 5;
	noutgoing[3] = 2;

	noutgoing[4] = 0;
	noutgoing[5] = 0;
	noutgoing[6] = 0;
	noutgoing[7] = 0;
	noutgoing[8] = 0;
	noutgoing[9] = 0;
	noutgoing[10] = 0;

	psrc[0] = 1;
	psrc[1] = 4;
	psrc[2] = 4;
	psrc[3] = 9;
	psrc[4] = 11;
	psrc[5] = 11;
	psrc[6] = 11;
	psrc[7] = 11;
	psrc[8] = 11;
	psrc[9] = 11;
	psrc[10] = 11;
	psrc[11] = 11;

	edgessrcdst[0] = 0;
	edgessrcdst[1] = 1;
	edgessrcdst[2] = 2;
	edgessrcdst[3] = 3;
	edgessrcdst[4] = 4;
	edgessrcdst[5] = 5;
	edgessrcdst[6] = 6;
	edgessrcdst[7] = 7;
	edgessrcdst[8] = 8;
	edgessrcdst[9] = 9;
	edgessrcdst[10] = 10;
}
unsigned Graph::readFromGR(char file[]) {
	std::ifstream cfile;
	cfile.open(file);

	// copied from GaloisCpp/trunk/src/FileGraph.h
	int masterFD = open(file, O_RDONLY);
	if (masterFD == -1) {
		printf("FileGraph::structureFromFile: unable to open %s.\n", file);
		return 1;
	}

	struct stat buf;
	int f = fstat(masterFD, &buf);
	if (f == -1) {
		printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
		abort();
	}
	size_t masterLength = buf.st_size;

	int _MAP_BASE = MAP_PRIVATE;
//#ifdef MAP_POPULATE
//  _MAP_BASE  |= MAP_POPULATE;
//#endif

	void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
	if (m == MAP_FAILED ) {
		m = 0;
		printf("FileGraph::structureFromFile: mmap failed.\n");
		abort();
	}

	//parse file
	uint64_t* fptr = (uint64_t*) m;
	__attribute__((unused))  uint64_t version = le64toh(*fptr++);
	assert(version == 1);
	uint64_t sizeEdgeTy = le64toh(*fptr++);
	uint64_t numNodes = le64toh(*fptr++);
	uint64_t numEdges = le64toh(*fptr++);
	uint64_t *outIdx = fptr;
	fptr += numNodes;
	uint32_t *fptr32 = (uint32_t*) fptr;
	uint32_t *outs = fptr32;
	fptr32 += numEdges;
	if (numEdges % 2)
		fptr32 += 1;
	unsigned *edgeData = (unsigned *) fptr32;

	// cuda.
	nnodes = numNodes;
	nedges = numEdges;

	printf("nnodes=%d, nedges=%d.\n", nnodes, nedges);
	allocOnHost();

	for (unsigned ii = 0; ii < nnodes; ++ii) {
		// fill unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst; foru *edgessrcwt;
		srcsrc[ii] = ii;
		IsParent[ii] = 1;
		childIndex[ii] = ii;
		childCount[ii] = 0;
		if (ii > 0) {
			psrc[ii] = le64toh(outIdx[ii - 1])+ 1;
			noutgoing[ii] = le64toh(outIdx[ii]) - le64toh(outIdx[ii - 1]);
		} else {
			psrc[0] = 1;
			noutgoing[0] = le64toh(outIdx[0]);
		}
		for (unsigned jj = 0; jj < noutgoing[ii]; ++jj) {
			unsigned edgeindex = psrc[ii] + jj;
			unsigned dst = le32toh(outs[edgeindex - 1]);
			if (dst >= nnodes)
				printf("\tinvalid edge from %d to %d at index %d(%d).\n", ii, dst, jj, edgeindex);
			edgessrcdst[edgeindex] = dst;
			edgessrcwt[edgeindex] = edgeData[edgeindex - 1];

			++nincoming[dst];
			//if (ii == 194 || ii == 352) {
			//	printf("edge %d: %d->%d, wt=%d.\n", edgeindex, ii, dst, edgessrcwt[edgeindex]);
			//}
		}
	}
	cfile.close(); // probably galois doesn't close its file due to mmap.
	return 0;
}
unsigned Graph::cudaCopy(struct Graph &copygraph) {
	copygraph.nnodes = nnodes;
	copygraph.nedges = nedges;

	copygraph.allocOnDevice();

	cudaMemcpy(copygraph.edgessrcdst, edgessrcdst, (nedges + 1) * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.edgessrcwt, edgessrcwt, (nedges + 1) * sizeof(foru), cudaMemcpyHostToDevice);

	cudaMemcpy(copygraph.psrc, psrc, (nnodes + 1) * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.noutgoing, noutgoing, nnodes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.nincoming, nincoming, nnodes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.srcsrc, srcsrc, nnodes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.parentNode, parentNode, nnodes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(copygraph.childCount, childCount, nnodes * sizeof(unsigned int), cudaMemcpyHostToDevice);
	return 0;
}

