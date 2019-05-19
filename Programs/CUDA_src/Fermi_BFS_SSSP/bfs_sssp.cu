/** Breadth-first search -*- CUDA -*-
 * @file
 * @section License
 *
 * Galois, a framework to exploit amorphous data-parallelism in irregular
 * programs.
 *
 * Copyright (C) 2013, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 *
 * @section Description
 *
 * Example breadth-first search application for demoing Galois system.
 *
 * @author Rupesh Nasre <nasre@ices.utexas.edu>
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <string.h>
#include <inttypes.h>
#include <iomanip>

#define XSTR(x) #x
#define STR(x) XSTR(x)

#define BLOCKSIZE	1024
#define BANKSIZE BLOCKSIZE

#define MYINFINITY	1000000000
#define BLOCKSIZE	1024
#define WORKPERTHREAD	1
#define VERTICALWORKPERTHREAD	1
using namespace std;

#define BFS
#define INT_SPLIT
#define EVERY_IT_CHILD_DIST_UPDATE

/*The following ENABLE_SERIAL_ALGO option is used for writing serial BFS/SSSP
 answer to BFS_S.txt/SSSP_S.txt*/
//#define ENABLE_SERIAL_ALGO


/*ENABLE_METRIC is for simulating block assignment to SM by measuring block executing time
(by the use of clock count). And shows the effectiveness of node splitting.
For this atomicMin operation is used for relaxation*/
//#define ENABLE_METRIC
//#define PRINT_BSM

/*ENABLE_IT_DIST is for measuring convergence of BFS/SSSP.
After each iteration program writes dist to file.
Then finally it calculates the convergence.*/
//#define ENABLE_IT_DIST

/*ENABLE_RTIME_DIST is for convergence rate. Here at first searial algo is executed first,
then after each iteration program calculates the convergence without writing result into a file.*/
//#define ENABLE_RTIME_DIST

//#define NO_ATOMIC_KERNEL
//#define DYNAMICP

//#define ENABLE_PRINT_DEGREE_DISTRIBUTION

#ifdef SSSP
	#define BFS_SSSP(index) edgessrcwt[index]
#else
	#define BFS_SSSP(index) 1
#endif

#include "graph.h"
extern SResult myresult;
#include "cudaKernels.h"
#include "metrics.h"
#include "convergence.h"
#include "myutils.h"

unsigned int NVERTICES;

int main(int argc, char *argv[])
{

	#ifdef SSSP
		strcpy(myresult.algo,"SSSP");
		printf("SSSP\n");fflush(stdout);
	#else
		strcpy(myresult.algo,"BFS");
		printf("BFS\n");fflush(stdout);
	#endif

	foru *dist;
	foru *hdist;
	foru floatzero = 0.0;
	unsigned int NBLOCKS, FACTOR = 128;

	Graph hgraph, graph;
	bool *changed, hchanged;
	int iteration = 0;

	clock_t starttime, endtime;
	int runtime;
	cudaEvent_t start, stop;
	float time;

	std::string str;

	cudaDeviceProp deviceProp;
	cudaDeviceReset();

	int dc;
	cudaGetDeviceCount(&dc);
	printf("dc=%d\n",dc);

	//cudaSetDevice(1);
	//int gc;cudaGetDevice(&gc);
	//printf("gc=%d\n",gc);

	//cudaError_t e=cudaDeviceSetLimit (cudaLimitDevRuntimePendingLaunchCount, 5000);
	//printf("Error:%d\n",e);

	cudaGetDeviceProperties(&deviceProp, 0);
	NBLOCKS = deviceProp.multiProcessorCount;
	printf("Device Name:%s\n",deviceProp.name);

	//cudaFuncSetCacheConfig(my_drelax, cudaFuncCachePreferShared);
	if (argc < 2) {
		printf("Usage: %s <file>\n", argv[0]);
		exit(1);
	}

	strcpy(myresult.exename,argv[0]);
	char gnm[100];
	getGraphname(gnm,argv[1]);
	cout<<gnm<<endl;
	strcpy(myresult.graphfile,gnm);

	cout << argv[1] <<endl;
	unsigned EnableNodeSplit=0;
	int NodeSplitLevel=1;
	unsigned w2f=0;
	char *fname;

	unsigned bucketSize = 10;
	if(argc>2){
		EnableNodeSplit=1;
		NodeSplitLevel = atoi(argv[2]);
		if (NodeSplitLevel==0)
			EnableNodeSplit=0;
		if(NodeSplitLevel==-1 && argc>3){
			bucketSize = atoi(argv[3]);
		}
	}

	if(NodeSplitLevel<0)
		myresult.bucketSize = bucketSize;
	else
		myresult.bucketSize = 0;

	if(NodeSplitLevel!=-1 && argc>3){
		w2f=1;
		fname=argv[3];
	}
	cudaGetLastError();

	hgraph.readFromGR(argv[1]);
	myresult.nnode=hgraph.nnodes;
	myresult.nedge=hgraph.nedges;
	unsigned oldnode=hgraph.nnodes;
	unsigned MaxDeg =hgraph.getMaxdeg();
	printf("MaxDeg=%d\n",MaxDeg);

	//hgraph.printGraph("BeforeSplit.txt");

#ifdef ENABLE_PRINT_DEGREE_DISTRIBUTION
	char tn[100];
	sprintf(tn,"DD_%s.txt",myresult.graphfile);
	hgraph.printfDegreeDistribution(tn);
#endif

	int set_parent_src=0;
	#ifdef ENABLE_RTIME_DIST
		#ifdef SSSP
			hgraph.findSerialSSSP(set_parent_src);
		#else
			hgraph.findSerialBFS(set_parent_src);
		#endif
	#endif

	#ifdef ENABLE_SERIAL_ALGO
		cout<<"Serial algo is enabled:Computing and writing into a file"<<endl;
		#ifdef SSSP
		cout<<"Writing SSSP_S.txt file"<<endl;
			hgraph.findSerialSSSPnPrint("SSSP_S.txt",set_parent_src);
		#else
			cout<<"Writing BFS_S.txt file"<<endl;
			hgraph.findSerialBFSnPrint("BFS_S.txt",set_parent_src);
		#endif
	#endif


	#ifdef DYNAMICP
		unsigned maxAllowed;
		if (NodeSplitLevel < 0) {
			maxAllowed =hgraph.getMaxAllowedOutDegree(bucketSize);
			printf("Approx Level=%f\n", (float) MaxDeg / maxAllowed);
		}
		else{
			unsigned temp=NodeSplitLevel;
			if(NodeSplitLevel ==0){
				temp = 1;
			}
			maxAllowed = (MaxDeg + temp - 1)/temp;
		}
	#else
		if(EnableNodeSplit) {
			hgraph.NodeSplit(NodeSplitLevel,bucketSize);
		}
		printf("After node split:\nnnodes=%d,nedges=%d (Extranodes=%d)\n",hgraph.nnodes,hgraph.nedges,hgraph.nnodes-oldnode);
	#endif

	//hgraph.printGraph("AfterSplit.txt");
	//exit(0);

	#ifdef ENABLE_RTIME_DIST
		printConv_rtime_init(NodeSplitLevel,myresult.maxdeg,gnm,myresult.algo);
	#endif

	myresult.ns_nnode=hgraph.nnodes;
	myresult.ns_nedge=hgraph.nedges;
	myresult.extranode=hgraph.nnodes-oldnode;


	#if defined(EVERY_IT_CHILD_DIST_UPDATE) && !defined(DYNAMIC)
		unsigned myParentCount=0;
		for(unsigned i=0;i<hgraph.nnodes;i++){
			if(hgraph.childCount[i]){
				myParentCount++;
			}
		}
		unsigned *hmyChildCount,*hmyChildIndex;
		unsigned *myChildCount,*myChildIndex;

		hmyChildCount=(unsigned *)malloc(sizeof(unsigned)*myParentCount);
		hmyChildIndex=(unsigned *)malloc(sizeof(unsigned)*myParentCount);
		unsigned index=0;
		for(unsigned i=0;i<hgraph.nnodes;i++){
			if(hgraph.childCount[i]){
				hmyChildCount[index]=hgraph.childCount[i];
				hmyChildIndex[index]=i;
				index++;
			}
		}
		if (cudaMalloc((void **)&myChildCount,myParentCount * sizeof(unsigned)) != cudaSuccess) CudaTest("allocating myChildCount failed");
		cudaMemcpy(myChildCount, hmyChildCount, myParentCount * sizeof(unsigned), cudaMemcpyHostToDevice);

		if (cudaMalloc((void **)&myChildIndex,myParentCount * sizeof(unsigned)) != cudaSuccess) CudaTest("allocating myChildIndex failed");
		cudaMemcpy(myChildIndex, hmyChildIndex, myParentCount * sizeof(unsigned), cudaMemcpyHostToDevice);

		unsigned DIST_FACTOR = (myParentCount + BLOCKSIZE * NBLOCKS - 1) / (BLOCKSIZE * NBLOCKS);
	#endif


	unsigned distLen=hgraph.nnodes;
	NVERTICES=hgraph.nnodes;

	hdist=(foru *)malloc(sizeof(foru)*distLen);

	/*
	for(unsigned i=0;i<distLen;i++){
		hdist[i]=MYINFINITY;
	}

	for(unsigned i=0;i<=hgraph.childCount[0];i++){
		hdist[i]=floatzero;
	}
*/

	//Setting distant 0 for src node and INF to other nodes
	int temp_par_count=-1;
	for(unsigned i=0;i<distLen;i++){
		if(hgraph.IsParent[i]){
			temp_par_count++;
		}
		if(temp_par_count==set_parent_src){
			hdist[i]=floatzero;
			//cout<<"i="<<i << " hgraph.childCount[i]="<<hgraph.childCount[i]<<endl;
			for(int k=1;k<=hgraph.childCount[i];k++){
				hdist[i+k]=floatzero;
			}
			i += hgraph.childCount[i];
		}
		else{
			hdist[i]=MYINFINITY;
		}
	}

	if (cudaMalloc((void **)&dist,distLen * sizeof(foru)) != cudaSuccess) CudaTest("allocating dist failed");
	cudaMemcpy(dist, hdist, distLen * sizeof(foru), cudaMemcpyHostToDevice);

	hgraph.cudaCopy(graph);

	#if 0
		volatile unsigned *arrayin, *arrayout;
		if (cudaMalloc((void **)&arrayin, NBLOCKS*sizeof(volatile unsigned)) != cudaSuccess) CudaTest("allocating arrayin failed");
		if (cudaMalloc((void **)&arrayout, NBLOCKS*sizeof(volatile unsigned)) != cudaSuccess) CudaTest("allocating arrayout failed");
		cudaMemset((void *)arrayin, 0, NBLOCKS*sizeof(volatile unsigned));
		cudaMemset((void *)arrayout, 0, NBLOCKS*sizeof(volatile unsigned));
	#endif

	FACTOR = (NVERTICES + BLOCKSIZE * NBLOCKS - 1) / (BLOCKSIZE * NBLOCKS);
	unsigned nblocks = NBLOCKS*FACTOR;

	#ifdef ENABLE_METRIC
		unsigned MAX_IT=32;
		unsigned *blk_counter=(unsigned *)malloc(MAX_IT * nblocks *  sizeof(unsigned));
		unsigned *hnb=(unsigned *)calloc(nblocks,sizeof(unsigned));
		unsigned *nb;
		if (cudaMalloc((void **)&nb,nblocks * sizeof(unsigned)) != cudaSuccess) CudaTest("allocating nb failed");
	#endif


	if (cudaMalloc((void **)&changed, sizeof(bool)) != cudaSuccess) CudaTest("allocating changed failed");
	printf("solving.\n");
	starttime = clock();
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

/*-------------------------------------------------------------------------------------------------------------------
 ***************************************** The Loop ****************************************************************
-------------------------------------------------------------------------------------------------------------------*/

#ifdef DYNAMICP
	unsigned *dSplitCount,SplitCount=0;
	if (cudaMalloc((void **)&dSplitCount, sizeof(unsigned)) != cudaSuccess) CudaTest("allocating arrayin failed");
	cudaMemcpy(dSplitCount, &SplitCount,  sizeof(unsigned), cudaMemcpyHostToDevice);
#endif

	#if defined(NO_ATOMIC_KERNEL) && !defined(EVERY_IT_CHILD_DIST_UPDATE)
		unsigned no_atomic_dist_check_count=0;
	#endif

	float k_time=0,sk_time=0,t_time;

//#ifdef DYNAMICP
//	dynamic_Relax_Recursive<<<1, 32>>>(dist, graph, changed, 0,set_parent_src);
//#endif

	do {
		#ifdef ENABLE_METRIC
			cudaMemcpy(nb, hnb, nblocks * sizeof(unsigned), cudaMemcpyHostToDevice);
		#endif

		++iteration;
		hchanged = false;
		cudaMemcpy(changed, &hchanged, sizeof(bool), cudaMemcpyHostToDevice);

		#if defined(NO_ATOMIC_KERNEL) && defined(EVERY_IT_CHILD_DIST_UPDATE) && !defined(DYNAMIC)
			if(DIST_FACTOR){
				unsigned dnblocks = NBLOCKS*DIST_FACTOR;
				cudaEvent_t k_start, k_stop;
				cudaEventCreate(&k_start);
				cudaEventCreate(&k_stop);
				cudaEventRecord(k_start, 0);

				my_updateDist <<<dnblocks, BLOCKSIZE>>> (dist,myChildCount,myChildIndex,myParentCount);

				cudaEventRecord(k_stop, 0);
				cudaEventSynchronize(k_stop);
				cudaEventElapsedTime(&t_time, k_start, k_stop);
				sk_time += t_time;
			}
		#endif

		cudaEvent_t k_start, k_stop;
		cudaEventCreate(&k_start);
		cudaEventCreate(&k_stop);
		cudaEventRecord(k_start, 0);

		#ifdef ENABLE_METRIC
			my_drelax_metric<<<nblocks/WORKPERTHREAD, BLOCKSIZE>>> (dist,graph,changed,nb);
		#else
			#ifdef DYNAMICP
				my_drelax_dynamicP <<<nblocks/WORKPERTHREAD, BLOCKSIZE>>> (dist,graph,changed,maxAllowed,dSplitCount);
			#elif defined(NO_ATOMIC_KERNEL)
				my_drelax_no_atomic <<<nblocks/WORKPERTHREAD, BLOCKSIZE>>> (dist,graph,changed);
			#else //Atomic Relax
				M_my_drelax <<<nblocks/WORKPERTHREAD, BLOCKSIZE>>> (dist,graph,changed);
			#endif
		#endif

		cudaEventRecord(k_stop, 0);
		cudaEventSynchronize(k_stop);
		cudaEventElapsedTime(&t_time, k_start, k_stop);
		k_time += t_time;

		#ifdef ENABLE_IT_DIST
			#ifdef ENABLE_RTIME_DIST
				cudaMemcpy(hdist, dist, distLen * sizeof(foru), cudaMemcpyDeviceToHost);
				printConv_rtime(hdist,hgraph,iteration);
			#else
				char distfname[100];
				sprintf(distfname,"Dist__L_%d_it_%d_%s_%s.txt",NodeSplitLevel,iteration-1,gnm,myresult.algo);
				cudaMemcpy(hdist, dist, distLen * sizeof(foru), cudaMemcpyDeviceToHost);
				Ary2File(distfname,hdist,hgraph.IsParent,distLen);
			#endif
		#endif

		CudaTest("solving failed");
		cudaMemcpy(&hchanged, changed, sizeof(bool), cudaMemcpyDeviceToHost);

		#ifdef ENABLE_METRIC
			if(iteration>MAX_IT){
				MAX_IT = MAX_IT *2;
				blk_counter = (unsigned *)realloc(blk_counter,MAX_IT * nblocks *  sizeof(unsigned));
			}
			//copy counter values
			cudaMemcpy(blk_counter + nblocks * (iteration-1), nb, nblocks * sizeof(unsigned), cudaMemcpyDeviceToHost);
		#endif


		#if !defined(DYNAMICP) && defined(NO_ATOMIC_KERNEL) && !defined(EVERY_IT_CHILD_DIST_UPDATE)
			if (hchanged == false) {

				cudaEvent_t k_start, k_stop;
				cudaEventCreate(&k_start);
				cudaEventCreate(&k_stop);
				cudaEventRecord(k_start, 0);

				no_atomic_dist_check<<<nblocks / WORKPERTHREAD, BLOCKSIZE>>>(dist, graph, changed);

				cudaEventRecord(k_stop, 0);
				cudaEventSynchronize(k_stop);
				cudaEventElapsedTime(&t_time, k_start, k_stop);
				sk_time += t_time;

				cudaMemcpy(&hchanged, changed, sizeof(bool), cudaMemcpyDeviceToHost);

				if (hchanged) {
					no_atomic_dist_check_count += 1;
				}
			}
		#endif

	} while (hchanged);
	cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
	endtime = clock();
	myresult.ktime=k_time;
	myresult.sktime=sk_time;

	#ifdef NO_ATOMIC_KERNEL
		#ifndef EVERY_IT_CHILD_DIST_UPDATE
			cout << "no_atomic_dist_check_count="<<no_atomic_dist_check_count<<endl;
		#endif
	#endif

	#ifdef DYNAMICP
		cudaMemcpy(&SplitCount, dSplitCount,  sizeof(unsigned), cudaMemcpyDeviceToHost);
		printf("SplitCount=%d(%d)\n",SplitCount,SplitCount/iteration);
	#endif
	/*-------------------------------------------------------------------------------------------------------------------
	 ***************************************** The End ****************************************************************
	-------------------------------------------------------------------------------------------------------------------*/

	printf("iterations = %d.\n", iteration);
	cudaMemcpy(hdist, dist, distLen * sizeof(foru), cudaMemcpyDeviceToHost);
	if (w2f){
		cout<<"Writing Result To file : "<<fname<<endl;
		fflush(stdout);
		Ary2File(fname,hdist,hgraph.IsParent,distLen);
	}
	verifysolution(hgraph.edgessrcdst, hgraph.edgessrcwt, hgraph.noutgoing,hgraph.psrc, hdist, NVERTICES, hgraph.nedges, hgraph.srcsrc);
	runtime = (int) (1000.0f * (endtime - starttime) / CLOCKS_PER_SEC);
	printf("%d ms.\n", runtime);
	printf("Time:%.2f ktime=%.2f sktime=%.2f\n",time,k_time,sk_time);
	fflush(stdout);

	myresult.iterations=iteration;
	myresult.runtime=runtime;
	myresult.time = time;
	myresult.splitlevel=NodeSplitLevel;
	// cleanup left to the OS.

	writeResult();

	#ifdef ENABLE_IT_DIST
		#ifdef ENABLE_RTIME_DIST
			foru max_dist;
			#ifdef SSSP
				max_dist = hgraph.max_sssp;
			#else
				max_dist = hgraph.max_bfs;
			#endif
			printConv_rtime_final(iteration,max_dist);
		#else
			char pattern[100];
			cout<<"Level:"<<NodeSplitLevel<<endl;
			sprintf(pattern,"Dist__L_%d_it_%%d_%s_%s.txt",NodeSplitLevel,gnm,myresult.algo);
			printConv(iteration-1,pattern);
		#endif
	#endif

	#ifdef ENABLE_METRIC
		char namepad[100];
		unsigned long tc=0;
		unsigned *aindex;
		for(int it=0;it<iteration;it++){
			sprintf(namepad,"L_%d_it_%d_%s_%s",NodeSplitLevel,it,gnm,myresult.algo);
			tc = tc + combineClocks(blk_counter+it*nblocks,nblocks,aindex,namepad);
		}
		printf("tc=%llu\n",tc);
		printf("Est time=%4.3f ms\n",tc/1147.0);
	#endif

	return 0;
}
