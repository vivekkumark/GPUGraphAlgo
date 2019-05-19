typedef struct Dmesh {


#ifdef DATA_DRIVEN
	//For DataDriven Approach
	unsigned *workList;
	unsigned *workList_Count;
	void Create_workList();
	__global__ void dinitWorklist(unsigned ntriangles);
#endif

	//unsigned subsetLevel;
	unsigned curralloc;
	unsigned currsizenodes;

	unsigned *pnnodes;
	FORD *nodex;
	FORD *nodey;

	unsigned *pntriangles;
	unsigned *tnodes;
	unsigned *neighbors;
	DIMSTYPE *neighboredges;

	DIMSTYPE *obtuse;
	bool *isbad;
	bool *isdel;
	unsigned *owner;

	bool *changed;
	bool *possible;

	Dmesh(unsigned nnodes, unsigned ntriangles);
	void ensureMemory(unsigned ntriangles, unsigned hnbad);

	unsigned get_ntriangles();
	bool get_changed();
	void init_changed();

	bool get_possible();
	void init_possible();
} Dmesh;

#ifdef DATA_DRIVEN
__global__ void Dmesh::dinitWorklist(unsigned ntriangles) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < ntriangles && !isdel[id] && isbad[id]) {
		unsigned myid = atomicInc(workList_Count,MAXID);
		workList[myid] = id;
	}
}
void Dmesh::Create_workList(){
	unsigned wl=0;
	cudaMemcpy(workList_Count, &wl, sizeof(unsigned), cudaMemcpyHostToDevice);
	unsigned ntriangles = get_ntriangles();

	unsigned NB = (ntriangles + _BS - 1)/_BS;
	dinitWorklist<<<NB,_BS>>>(ntriangles);
}
#endif

void Dmesh::init_changed(){
	bool hchanged = false;
	cudaMemcpy(changed, &hchanged, sizeof(bool), cudaMemcpyHostToDevice);
}
void Dmesh::init_possible(){
	bool hpossible = false;
	cudaMemcpy(possible, &hpossible, sizeof(bool), cudaMemcpyHostToDevice);
}
unsigned Dmesh::get_ntriangles() {
	unsigned hntriangles;
	cudaMemcpy(&hntriangles, pntriangles, sizeof(unsigned), cudaMemcpyDeviceToHost);
	return hntriangles;
}

bool Dmesh::get_changed() {
	bool hchanged;
	cudaMemcpy(&hchanged, changed, sizeof(bool), cudaMemcpyDeviceToHost);
	return hchanged;
}

bool Dmesh::get_possible() {
	bool hpossible;
	cudaMemcpy(&hpossible, possible, sizeof(bool), cudaMemcpyDeviceToHost);
	return hpossible;
}

Dmesh::Dmesh(unsigned nnodes, unsigned ntriangles) {

	if (cudaMalloc((void **) &possible, sizeof(bool)) != cudaSuccess)
		CudaTest("allocating possible failed");
	if (cudaMalloc((void **) &changed, sizeof(bool)) != cudaSuccess)
		CudaTest("allocating changed failed");

	if (cudaMalloc((void **) &pnnodes, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating pnnodes failed");
	if (cudaMalloc((void **) &pntriangles, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating pntriangles failed");

	cudaMemcpy(pnnodes, &nnodes, sizeof(unsigned), cudaMemcpyHostToDevice);
	cudaMemcpy(pntriangles, &ntriangles, sizeof(unsigned), cudaMemcpyHostToDevice);

	curralloc = ALLOCFACTOR * ntriangles;
	currsizenodes = ALLOCFACTOR * nnodes;

	if (cudaMalloc((void **) &nodex, ALLOCFACTOR * nnodes * sizeof(FORD)) != cudaSuccess)
		CudaTest("allocating nodex failed");
	if (cudaMalloc((void **) &nodey, ALLOCFACTOR * nnodes * sizeof(FORD)) != cudaSuccess)
		CudaTest("allocating nodey failed");

	if (cudaMalloc((void **) &tnodes, ALLOCFACTOR * 3 * ntriangles * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating tnodes failed");
	if (cudaMalloc((void **) &neighbors, ALLOCFACTOR * 3 * ntriangles * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating neighbors failed");
	if (cudaMalloc((void **) &neighboredges, ALLOCFACTOR * 3 * ntriangles * sizeof(DIMSTYPE)) != cudaSuccess)
		CudaTest("allocating neighboredges failed");

	if (cudaMalloc((void **) &isbad, ALLOCFACTOR * ntriangles * sizeof(bool)) != cudaSuccess)
		CudaTest("allocating isbad failed");
	if (cudaMalloc((void **) &obtuse, ALLOCFACTOR * ntriangles * sizeof(DIMSTYPE)) != cudaSuccess)
		CudaTest("allocating obtuse failed");
	if (cudaMalloc((void **) &isdel, ALLOCFACTOR * ntriangles * sizeof(bool)) != cudaSuccess)
		CudaTest("allocating isdel failed");
	if (cudaMalloc((void **) &owner, ALLOCFACTOR * ntriangles * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating owner failed");

#ifdef DATA_DRIVEN
	if (cudaMalloc((void **) &workList, ALLOCFACTOR * 3 * ntriangles * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating workList failed");
	if (cudaMalloc((void **) &workList_Count, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating workList failed");
#endif

}

void Dmesh::ensureMemory(unsigned ntriangles, unsigned hnbad) {
	if (ntriangles + 2 * hnbad > curralloc) { // here 2 is the no of new triangles added for each bad triangle.
		nodex = (FORD *) mycudarealloc(nodex, currsizenodes * sizeof(FORD), ALLOCMULTIPLE * currsizenodes * sizeof(FORD));
		nodey = (FORD *) mycudarealloc(nodey, currsizenodes * sizeof(FORD), ALLOCMULTIPLE * currsizenodes * sizeof(FORD));
		currsizenodes = ALLOCMULTIPLE * currsizenodes;

		tnodes = (unsigned *) mycudarealloc(tnodes, 3 * curralloc * sizeof(unsigned), ALLOCMULTIPLE * 3 * curralloc * sizeof(unsigned));
		neighbors = (unsigned *) mycudarealloc((unsigned *) neighbors, 3 * curralloc * sizeof(unsigned), ALLOCMULTIPLE * 3 * curralloc * sizeof(unsigned));
		neighboredges = (unsigned *) mycudarealloc((unsigned *) neighboredges, 3 * curralloc * sizeof(DIMSTYPE), ALLOCMULTIPLE * 3 * curralloc * sizeof(DIMSTYPE));

		isbad = (bool *) mycudarealloc(isbad, curralloc * sizeof(bool), ALLOCMULTIPLE * curralloc * sizeof(bool));
		obtuse = (DIMSTYPE *) mycudarealloc(obtuse, curralloc * sizeof(DIMSTYPE), ALLOCMULTIPLE * curralloc * sizeof(DIMSTYPE));
		isdel = (bool *) mycudarealloc(isdel, curralloc * sizeof(bool), ALLOCMULTIPLE * curralloc * sizeof(bool));
		owner = (unsigned *) mycudarealloc((void *) owner, curralloc * sizeof(unsigned), ALLOCMULTIPLE * curralloc * sizeof(unsigned));

		#ifdef DATA_DRIVEN
			workList = (unsigned *) mycudarealloc(workList, 3 * curralloc * sizeof(unsigned), ALLOCMULTIPLE * 3 * curralloc * sizeof(unsigned));
		#endif

		curralloc *= ALLOCMULTIPLE;
		printf("\t\tallocating memory to %d.\n", curralloc);
	}
}
