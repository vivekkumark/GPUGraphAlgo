#define FORD float

typedef struct Cavity {

	unsigned MAXCAVITY;
	unsigned *cavity;
	unsigned *frontier;
	unsigned *pre;
	unsigned *post;
	unsigned *connections;

	unsigned *iicavity;
	unsigned *iifrontier;
	unsigned *iipre;
	unsigned *iipost;
	unsigned *iiconnections;
	unsigned *centerelement;
	unsigned *marker;
	bool *cavityProcessed;

	unsigned *activeCavityCount;

	FORD *centerx;
	FORD *centery;

	unsigned *post_bad_work;
	unsigned *post_bad_work_count;

	void init_post_bad_work_count();
	unsigned get_post_bad_work_count();

	Cavity(unsigned numCavity);
	void set_activeCavityCount(unsigned val);
	unsigned get_activeCavityCount();

}Cavity;

unsigned Cavity::get_post_bad_work_count(){
	unsigned temp;
	cudaMemcpy(&temp, post_bad_work_count, sizeof(unsigned), cudaMemcpyDeviceToHost);
	return temp;
}

void  Cavity::init_post_bad_work_count(){
	unsigned hpost_bad_work_count=0;
	cudaMemcpy((void *)post_bad_work_count, &hpost_bad_work_count, sizeof(unsigned), cudaMemcpyHostToDevice);
}
void Cavity::set_activeCavityCount(unsigned val){
	unsigned hval=0;
	cudaMemcpy((void *)activeCavityCount, &hval, sizeof(unsigned), cudaMemcpyHostToDevice);
}

unsigned Cavity::get_activeCavityCount(){
	unsigned hactiveCavityCount;
	cudaMemcpy(&hactiveCavityCount, activeCavityCount, sizeof(unsigned), cudaMemcpyDeviceToHost);
	return hactiveCavityCount;
}

Cavity::Cavity(unsigned numCavity){

	printf("numCavity=%d\n",numCavity);
	MAXCAVITY = numCavity;
	if (cudaMalloc((void **) &iicavity, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating iicavity failed");
	if (cudaMalloc((void **) &iifrontier, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating iifrontier failed");
	if (cudaMalloc((void **) &iipre, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating iipre failed");
	if (cudaMalloc((void **) &iipost, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating iipost failed");
	if (cudaMalloc((void **) &iiconnections, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating iiconnections failed");
	if (cudaMalloc((void **) &centerelement, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating centerelement failed");
	if (cudaMalloc((void **) &centerx, numCavity * sizeof(FORD)) != cudaSuccess)
		CudaTest("allocating centerx failed");
	if (cudaMalloc((void **) &centery, numCavity * sizeof(FORD)) != cudaSuccess)
		CudaTest("allocating centery failed");
	if (cudaMalloc((void **) &marker, numCavity * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating marker failed");
	if (cudaMalloc((void **) &cavityProcessed, numCavity * sizeof(bool)) != cudaSuccess)
		CudaTest("allocating cavityProcessed failed");

	if (cudaMalloc((void **) &activeCavityCount,sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating activeCavityCount failed");

	if (cudaMalloc((void **) &cavity, numCavity * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating cavity failed");
	if (cudaMalloc((void **) &frontier, numCavity * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating frontier failed");
	if (cudaMalloc((void **) &pre, numCavity * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating pre failed");
	if (cudaMalloc((void **) &post, numCavity * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating post failed");
	if (cudaMalloc((void **) &connections, numCavity * 4 * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating connections failed");

	if (cudaMalloc((void **) &post_bad_work_count,sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating post_bad_work_count failed");

	if (cudaMalloc((void **) &post_bad_work, numCavity * CAVITY_SIZE * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating post_bad_work failed");

}
