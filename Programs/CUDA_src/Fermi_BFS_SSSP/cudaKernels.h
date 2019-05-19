__global__
void my_drelax_metric(foru *dist, Graph graph, bool *changed, unsigned *nb) {

	clock_t start_time;
	if (threadIdx.x == 0)
		start_time = clock64();

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src]) {
				//src=graph.srcsrc[src];
				unsigned int start = graph.psrc[src];
				unsigned int end = start + graph.noutgoing[src];
				// go over all the target nodes for the source node.
				unsigned int u = src;
				for (ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii]; // target node.
#ifdef SSSP
					foru wt = graph.edgessrcwt[ii];
#else
					foru wt = 1;
#endif

					foru alt = dist[u] + wt;
					if (alt < dist[v]) {
						atomicMin(&dist[v], alt);
						if (graph.childCount[v]) {
							for (unsigned k = 1; k <= graph.childCount[v]; k++) {
								atomicMin(&dist[v + k], alt);
							}
						}
						if (++nprocessed < VERTICALWORKPERTHREAD) {
							changedv[iichangedv] = v;
							iichangedv += BANKSIZE;
						}
					}

				}
			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}

	__syncthreads();
	if (threadIdx.x == 0) {
		clock_t stop_time = clock64();
		nb[blockIdx.x] = (stop_time - start_time) / 1000;
	}
}

__global__
void my_drelax(foru *dist, Graph graph, bool *changed) {

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src]) {
				//src=graph.srcsrc[src];
				unsigned int start = graph.psrc[src];
				unsigned int end = start + graph.noutgoing[src];
				// go over all the target nodes for the source node.
				unsigned int u = src;
				for (ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii]; // target node.
#ifdef SSSP
					foru wt = graph.edgessrcwt[ii];
#else
					foru wt = 1;
#endif

					foru alt = dist[u] + wt;
					if (alt < dist[v]) {
						dist[v] = alt;
						if (++nprocessed < VERTICALWORKPERTHREAD) {
							// add work to the worklist.
							changedv[iichangedv] = v;
							iichangedv += BANKSIZE;
						}
					}

				}
			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}
}

__global__
void my_updateDist(foru *dist, unsigned *myChildCount, unsigned *myChildIndex, unsigned myParentCount) {
	unsigned nn = (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = myParentCount;
	if (nn < nv) {
		unsigned index = myChildIndex[nn];
		unsigned chcount = myChildCount[nn];
		foru mydist = dist[index];
		for (unsigned i = 1; i <= chcount; i++) {
			dist[index + i] = mydist;
		}
	}
}

__global__
void no_atomic_dist_check(foru *dist, Graph graph, bool *changed) {
	unsigned nn = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nv = graph.nnodes;
	if (nn < nv) {
		unsigned cc = graph.childCount[nn];
		if (cc) {
			foru nndist = dist[nn];
			for (unsigned k = 1; k <= cc; k++) {
				if (dist[nn + k] != nndist) {
					dist[nn + k] = nndist;
					*changed = true;
				}
			}
		}
	}
}

__global__
void my_drelax_no_atomic(foru *dist, Graph graph, bool *changed) {

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src]) {
				//src=graph.srcsrc[src];
				unsigned int start = graph.psrc[src];
				unsigned int end = start + graph.noutgoing[src];
				// go over all the target nodes for the source node.
				unsigned int u = src;
				for (ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii]; // target node.
#ifdef SSSP
					foru wt = graph.edgessrcwt[ii];
#else
					foru wt = 1;
#endif

					foru alt = dist[u] + wt;
					if (alt < dist[v]) {
						dist[v] = alt;
						if (graph.childCount[v]) {
							for (unsigned k = 1; k <= graph.childCount[v]; k++) {
								dist[v + k] = alt;
							}
						}
						if (++nprocessed < VERTICALWORKPERTHREAD) {
							changedv[iichangedv] = v;
							iichangedv += BANKSIZE;
						}
					}

				}
			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}
}



__global__
void my_drelax_no_atomic_level(foru *dist, Graph graph, bool *changed,unsigned level) {

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src] && dist[src]==level) {
				//src=graph.srcsrc[src];
				unsigned int start = graph.psrc[src];
				unsigned int end = start + graph.noutgoing[src];
				// go over all the target nodes for the source node.
				unsigned int u = src;
				for (ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii]; // target node.
				#ifdef SSSP
					foru wt = graph.edgessrcwt[ii];
				#else
					foru wt = 1;
				#endif

					foru alt = dist[u] + wt;
					if (alt < dist[v]) {
						dist[v] = alt;
						if (graph.childCount[v]) {
							for (unsigned k = 1; k <= graph.childCount[v]; k++) {
								dist[v + k] = alt;
							}
						}
						if (++nprocessed < VERTICALWORKPERTHREAD) {
							changedv[iichangedv] = v;
							iichangedv += BANKSIZE;
						}
					}

				}
			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}
}

__global__
void M_my_drelax(foru *dist, Graph graph, bool *changed) {

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src]) {
				//src=graph.srcsrc[src];
				unsigned int start = graph.psrc[src];
				unsigned int end = start + graph.noutgoing[src];
				// go over all the target nodes for the source node.
				unsigned int u = src;
				for (ii = start; ii < end; ++ii) {
					unsigned int v = graph.edgessrcdst[ii]; // target node.
					#ifdef SSSP
						foru wt = graph.edgessrcwt[ii];
					#else
						foru wt = 1;
					#endif

					foru alt = dist[u] + wt;
					if (alt < dist[v]) {
						atomicMin(&dist[v], alt);
						if (graph.childCount[v]) {
							for (unsigned k = 1; k <= graph.childCount[v]; k++) {
								atomicMin(&dist[v + k], alt);
							}
						}
						if (++nprocessed < VERTICALWORKPERTHREAD) {
							changedv[iichangedv] = v;
							iichangedv += BANKSIZE;
						}
					}

				}
			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}
}
#if 0
__global__
void dynamicP_relax_coal(foru *dist, Graph graph, bool *changed,unsigned u,unsigned mystart,unsigned nout){

	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned start = mystart + id;
	unsigned myend = mystart + nout;

	bool flag =false;
	for(unsigned ii=start;ii<myend;ii+=nthreads){
		unsigned int v = graph.edgessrcdst[ii];
		#ifdef SSSP
			foru wt = graph.edgessrcwt[ii];
		#else
			foru wt = 1;
		#endif
		foru alt = dist[u] + wt;
		if (alt < dist[v]){
			#ifndef NO_ATOMIC_KERNEL
				atomicMin(&dist[v], alt);
			#else
				dist[v] = alt;
			#endif
			flag = true;
		}
	}

	if (flag) {
		*changed = true;
	}
}

#define MAX_DEPTH 2
__global__
void dynamic_Relax_Recursive(foru *dist, Graph graph, bool *changed,unsigned depth,unsigned src){

	unsigned start = graph.psrc[src];
	unsigned nout = graph.noutgoing[src];
	unsigned end = start + nout;
	unsigned nthreads = blockDim.x*gridDim.x;
	unsigned id=threadIdx.x + blockIdx.x * blockDim.x;

	if(depth>MAX_DEPTH)
		return;

	unsigned myWork = start + id;

	for(unsigned ii=myWork;ii<end;ii+=nthreads){
		#ifdef SSSP
			foru wt = graph.edgessrcwt[ii];
		#else
			foru wt = 1;
		#endif

		unsigned int v = graph.edgessrcdst[ii];
		foru alt = dist[src] + wt;
		if (alt < dist[v]) {
			#ifndef NO_ATOMIC_KERNEL
				atomicMin(&dist[v], alt);
			#else
				dist[v] = alt;
			#endif
			*changed = true;

			cudaStream_t St;
			cudaStreamCreateWithFlags(&St, cudaStreamNonBlocking);
			dynamic_Relax_Recursive<<<1, 64, 0, St>>>(dist, graph, changed, depth+1,v);
			//if(cudaSuccess != cudaGetLastError())
			//	printf("p\n");
			cudaStreamDestroy(St);
		}
	}
}


__global__
void dynamicP_relax_2(foru *dist, Graph graph, bool *changed,unsigned u,unsigned mystart,unsigned nout){

	unsigned int v = graph.edgessrcdst[mystart+threadIdx.x];
	#ifdef SSSP
		foru wt = graph.edgessrcwt[mystart+threadIdx.x];
	#else
		foru wt = 1;
	#endif
	foru alt = dist[u] + wt;
	if (alt < dist[v]){
		dist[v] = alt;
		*changed = true;
	}
}

__global__
void my_drelax_dynamicP(foru *dist, Graph graph, bool *changed,unsigned maxDegreeAllowed,unsigned *dSplitCount) {

	unsigned int workperthread = WORKPERTHREAD;
	unsigned nn = workperthread * (blockIdx.x * blockDim.x + threadIdx.x);
	unsigned nv = graph.nnodes;
	unsigned int ii;
	__shared__ int changedv[VERTICALWORKPERTHREAD * BLOCKSIZE];
	int iichangedv = threadIdx.x;
	int anotheriichangedv = iichangedv;
	unsigned int nprocessed = 0;

	// collect the work to be performed.
	for (unsigned node = 0; node < workperthread; ++node, ++nn) {
		changedv[iichangedv] = nn;
		iichangedv += BANKSIZE;
	}

	// go over the worklist and keep updating it in a BFS manner.
	while (anotheriichangedv < iichangedv) {
		nn = changedv[anotheriichangedv];
		anotheriichangedv += BANKSIZE;
		if (nn < nv) {
			unsigned src = nn; // source node.
			if (graph.psrc[src]) {
				unsigned nout = graph.noutgoing[src];
				unsigned int start = graph.psrc[src];

				bool dpfailed=true;
				if(nout>maxDegreeAllowed){
					//printf("here\n");
					//atomicInc(dSplitCount,MYINFINITY);
					cudaError_t e ;
					#define DP_STREAM
					#ifdef DP_STREAM
						cudaStream_t St;
						cudaStreamCreateWithFlags(&St, cudaStreamNonBlocking);
						dynamicP_relax_coal <<<1,32,0,St>>> (dist,graph,changed,src,start,nout);
						e = cudaGetLastError();
						cudaStreamDestroy(St);
					#else
						dynamicP_relax <<<1,32>>> (dist,graph,changed,src,start,nout);
						e = cudaGetLastError();
					#endif
					//dynamicP_relax_2 <<<1,nout>>> (dist,graph,changed,src,start,nout);
					dpfailed=(cudaSuccess != e);
					//if(dpfailed)
					//	printf("p");
				}
				if(dpfailed){
					unsigned int end = start + nout;
					// go over all the target nodes for the source node.
					unsigned int u = src;
					for (ii = start; ii < end; ++ii) {
						unsigned int v = graph.edgessrcdst[ii];
						#ifdef SSSP
							foru wt = graph.edgessrcwt[ii];
						#else
							foru wt = 1;
						#endif

						foru alt = dist[u] + wt;
						if (alt < dist[v]) {
							#ifndef NO_ATOMIC_KERNEL
								atomicMin(&dist[v], alt);
							#else
								dist[v] = alt;
							#endif
							if (++nprocessed < VERTICALWORKPERTHREAD) {
								changedv[iichangedv] = v;
								iichangedv += BANKSIZE;
							}
						}
					}
				}

			}
		}
	}
	if (nprocessed) {
		*changed = true;
	}
}
#endif
