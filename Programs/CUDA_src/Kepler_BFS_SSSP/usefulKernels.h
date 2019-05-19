__global__
void manyblk_my_drelax_oneit(foru *dist, Graph graph, bool *changed) {

	unsigned totalVertices = graph.nnodes;
	unsigned myver = (blockIdx.x * blockDim.x + threadIdx.x);
	int flag = 0;

	if (myver < totalVertices) {

		unsigned int start = graph.psrc[myver];
		unsigned int end = start + graph.noutgoing[myver];

		for (int ii = start; ii < end; ++ii) {
			unsigned int v = graph.edgessrcdst[ii];
			#ifdef SSSP
				foru wt = graph.edgessrcwt[ii];
			#else
				foru wt = 1;
			#endif

			foru alt = dist[myver] + wt;
			if (alt < dist[v]) {
				dist[v] = alt;
				flag = 1;
			}
		}

	}

	if (flag) {
		*changed = true;
	}
}

__global__
void my_drelax_oneit(float *dist, Graph graph, bool *changed) {

	unsigned totalVertices = graph.nnodes;
	unsigned bid = blockIdx.x;
	unsigned totalSM = gridDim.x;
	unsigned totalIteration = (totalVertices + totalSM * blockDim.x - 1) / (totalSM * blockDim.x);

	int flag = 0;

	for (int it = 0; it < totalIteration; it++) {
		unsigned myver = (bid * blockDim.x + threadIdx.x);
		bid = bid + totalSM;
		if (myver < totalVertices) {
			unsigned int start = graph.psrc[myver];
			unsigned int end = start + graph.noutgoing[myver];

			for (int ii = start; ii < end; ++ii) {
				unsigned int v = graph.edgessrcdst[ii];
				#ifdef SSSP
					float wt = graph.edgessrcwt[ii];
				#else
					float wt = 1;
				#endif

				float alt = dist[myver] + wt;
				if (alt < dist[v]) {
					dist[v] = alt;
					flag = 1;
				}
			}
		}
	}

	if (flag) {
		*changed = true;
	}
}

__device__
void global_sync(unsigned goalVal, volatile unsigned *Arrayin, volatile unsigned *Arrayout) {
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
void my_drelax_oneitA(float *dist, Graph graph, volatile bool *changed, volatile unsigned *arrayin, volatile unsigned *arrayout) {

	unsigned totalVertices = graph.nnodes;
	unsigned totalSM = gridDim.x;
	unsigned totalIteration = (totalVertices + totalSM * blockDim.x - 1) / (totalSM * blockDim.x);
	unsigned int globalit = 1;

	while (1) {
		volatile int flag = 0;
		unsigned bid = blockIdx.x;
		for (int it = 0; it < totalIteration; it++) {
			unsigned myver = (bid * blockDim.x + threadIdx.x);
			bid = bid + totalSM;
			if (myver < totalVertices) {
				unsigned int start = graph.psrc[myver];
				unsigned int end = start + graph.noutgoing[myver];

				for (int ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii];
#ifdef SSSP
					float wt = graph.edgessrcwt[ii];
#else
					float wt = 1;
#endif

					float alt = dist[myver] + wt;
					if (alt < dist[v]) {
						dist[v] = alt;
						flag = 1;
					}
				}
			}
		}
		if (flag) {
			*changed = true;
			__threadfence();
		}

		global_sync(globalit, arrayin, arrayout);
		globalit++;

		if (*changed == false)
			break;

		global_sync(globalit, arrayin, arrayout);
		globalit++;

		*changed = false;
		__threadfence();

		global_sync(globalit, arrayin, arrayout);
		globalit++;

	}

}

