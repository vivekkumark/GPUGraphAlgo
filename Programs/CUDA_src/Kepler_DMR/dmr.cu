/** Delaunay refinement -*- CUDA -*-
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
 * Refinement of an initial, unrefined Delaunay mesh to eliminate triangles
 * with angles < 30 degrees, using a variation of Chew's algorithm.
 *
 * @author Rupesh Nasre <nasre@ices.utexas.edu>
 */

#include "common.h"

#define VTKI 0
#define VTKO 0

#define MINANGLE	25
#define PI		3.14159265358979323846	// from C99 standard.
#define FORD		float
#define DIMSTYPE	unsigned

#define INVALIDID	1234567890
#define MAXID		INVALIDID
#define TESTNBLOCKSFACTOR	1

#define ALLOCMULTIPLE	2	// alloc in multiples of this.
unsigned ALLOCFACTOR = 6; // initial alloc factor.

#include "helpers.h"
#include "myresult.h"
#include "NeighbournVTK.h"
#include "cavity.h"
#include "dmesh.h"
#include "m_refine.h"
#include "refine.h"

int main(int argc, char *argv[]) {
	unsigned NBLOCKS;
	unsigned int ntriangles, nnodes;
	bool hchanged;
	FORD *hnodex, *hnodey;
	unsigned *htnodes;
	int iteration = 0;
	unsigned intzero = 0;
	volatile unsigned *go;
	volatile unsigned *arrayin, *arrayout;
	unsigned *bcount, *nbad;
	unsigned hnbad;
	bool hpossible;

	clock_t starttime, endtime;
	int runtime;

	std::string str;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	unsigned _SM = deviceProp.multiProcessorCount;
	unsigned _MTPB = deviceProp.maxThreadsPerBlock;

	printf("SM Count=%d BLOCKSIZE=%d _MTPB=%d\n", _SM, BLOCKSIZE, _MTPB);

	//cudaFuncSetCacheConfig(drefine, cudaFuncCachePreferL1);
	if (argc != 2) {
		printf("Usage: %s <basefilename>\n", argv[0]);
		exit(1);
	}

	getname(myresult.fname, argv[0]);
	getname(myresult.ipname, argv[1]);

	cudaGetLastError();

	std::cout << "reading graphs...\n";
	readNodes(argv[1], hnodex, hnodey, nnodes);
	std::cout << "\t" << nnodes << " nodes\n";
	readTriangles(argv[1], htnodes, ntriangles, nnodes);
	std::cout << "\t" << ntriangles << " triangles.\n";

	vivek me;
#if VTKI
	cout << "writing vtk data file" << endl;
	me.print2file_unstructuredgrid("dmr.vtk", hnodex, hnodey, nnodes, htnodes, ntriangles);
#endif

	Dmesh myDmesh(nnodes, ntriangles);

	cudaMemcpy(myDmesh.nodex, hnodex, nnodes * sizeof(FORD), cudaMemcpyHostToDevice);
	cudaMemcpy(myDmesh.nodey, hnodey, nnodes * sizeof(FORD), cudaMemcpyHostToDevice);
	cudaMemcpy(myDmesh.tnodes, htnodes, 3 * ntriangles * sizeof(unsigned), cudaMemcpyHostToDevice);

	unsigned ntriperit = _SM * _BS;
	unsigned ntriit = (ntriangles + ntriperit - 1) / ntriperit;

#if 1
	if (me.findneigh_file_check(argv[1], myDmesh.neighboredges, myDmesh.neighbors, ntriangles) == false) {
		for (unsigned ii = 0; ii < ntriit; ++ii) {
			printf("finding neighbors: %3d%% complete.\r", (int) (ii * ntriperit * 100.0 / ntriangles));
			fflush(stdout);
			dfindneighbors<<<_SM, _BS>>>(myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.neighbors, myDmesh.neighboredges, nnodes, ntriangles, _SM, ii * ntriperit,
					(ii + 1) * ntriperit);
			CudaTest("find neighbors failed");
		}
		me.write_neigh_to_file(argv[1], myDmesh.neighboredges, myDmesh.neighbors, ntriangles);
	}
#else
	me.findneighbor_serial(argv[1], neighboredges, neighbors, htnodes, ntriangles);
#endif

	printf("\n");
	printf("init.\n");

	NBLOCKS = (ntriangles + _BS - 1) / _BS;
	dinit<<<NBLOCKS, _BS>>>(myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.isbad, myDmesh.obtuse, myDmesh.isdel, nnodes, ntriangles);
	CudaTest("initialization failed");

	if (cudaMalloc((void **) &go, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating go failed");
	cudaMemcpy((void *) go, &intzero, sizeof(unsigned), cudaMemcpyHostToDevice);

	unsigned nblockfactor = TESTNBLOCKSFACTOR; //(ntriangles < 1000000 ? 7 : (ntriangles < 10000000 ? 31 : 61));	// for 250k.2, use 7, for r1M use 31, for r5M use 61.
	unsigned nblocks = _SM * nblockfactor;
	unsigned blocksize = BLOCKSIZE;
	//bool hlchanged;

	if (cudaMalloc((void **) &arrayin, nblocks * sizeof(volatile unsigned)) != cudaSuccess)
		CudaTest("allocating arrayin failed");
	if (cudaMalloc((void **) &arrayout, nblocks * sizeof(volatile unsigned)) != cudaSuccess)
		CudaTest("allocating arrayout failed");
	if (cudaMalloc((void **) &bcount, nblocks * sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating blockcount failed");

	if (cudaMalloc((void **) &nbad, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating nbad failed");

	printf("solving.\n");
	starttime = clock();

	unsigned R_NBLOCKS = nblocks;
	unsigned R_BLOCKSIZE = blocksize;

	float k_time = 0, t_time;
	unsigned hndel;

	countbad<<<nblocks, blocksize>>>(myDmesh.isdel, myDmesh.isbad, ntriangles, nbad, 1000 + iteration, arrayin, arrayout, bcount);
	CudaTest("countbad failed");
	cudaMemcpy(&hnbad, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

	countnotdel<<<nblocks, blocksize>>>(myDmesh.isdel, ntriangles, nbad, 1001 + iteration, arrayin, arrayout, bcount);
	CudaTest("countdel failed");
	cudaMemcpy(&hndel, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

	cout << "\n\n hndel=" << hndel << " nbad=" << hnbad << endl;

	myresult.i_ntriangles = ntriangles;
	myresult.i_bad = hnbad;

	Cavity mycav(nblocks * blocksize * 50);
//	cudaError_t e=cudaDeviceSetLimit (cudaLimitDevRuntimePendingLaunchCount, nblocks * blocksize * 50);
//	printf("Error:%d\n",e);


	//myDmesh.subsetLevel = 0;
	do {

		//myDmesh.subsetLevel++;
		++iteration;
		unsigned orintriangles = ntriangles;
		myDmesh.init_changed();
		myDmesh.ensureMemory(ntriangles, hnbad);

#ifdef MKER

		if (R_NBLOCKS == 1 && R_BLOCKSIZE == 1) {
			ntriperit = ntriangles;
		} else {
			ntriperit = R_NBLOCKS * R_BLOCKSIZE * 50;
		}
		//printf("ntriperit=%d R_NBLOCKS=%d R_BLOCKSIZE=%d \n", ntriperit, R_NBLOCKS, R_BLOCKSIZE);
		ntriit = (ntriangles + ntriperit - 1) / ntriperit;

		for (unsigned ii = 0; ii < ntriit; ++ii) {
			mycav.set_activeCavityCount(0);
			//printf("-------------------------------------------------------------ii=%d/%d It=%d\n", ii, ntriit, iteration);
			fflush(stdout);

			cudaEvent_t k_start, k_stop;
			cudaEventCreate(&k_start);
			cudaEventCreate(&k_stop);
			cudaEventRecord(k_start, 0);

	#ifdef ENABLE_ATOMIC_MARKING
				ntriangles = myDmesh.get_ntriangles();
				cudaMemset(myDmesh.owner,0xFF,sizeof(unsigned)*ntriangles);
	#endif
			m_refine_cavityMark<<<R_NBLOCKS, R_BLOCKSIZE>>>(ii * ntriperit, (ii + 1) * ntriperit, myDmesh, mycav);
			if (CudaTest("solving failed m_refine_cavityMark"))
				exit(0);

			unsigned hactiveCavityCount = mycav.get_activeCavityCount();
			unsigned N_TH = 128;
			unsigned N_BLK = (hactiveCavityCount + N_TH - 1) / N_TH;
			unsigned Ccount = hactiveCavityCount;

			if (hactiveCavityCount) {
				#ifndef ENABLE_ATOMIC_MARKING
					m_refine_cavityReMark<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
					if (CudaTest("solving failed m_refine_cavityReMark"))
						exit(0);
				#endif

				#ifdef ENABLE_POST_WORK
					mycav.init_post_bad_work_count();
				#endif
				m_refine_cavityReTriangulate<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
				if (CudaTest("solving failed m_refine_cavityReTriangulate"))
					exit(0);
				#ifdef ENABLE_POST_WORK
					unsigned post_wc=mycav.get_post_bad_work_count();
					if(post_wc){
						unsigned pN_TH = 128;
						unsigned pN_BLK = (post_wc + N_TH - 1) / (N_TH);
						post_work_kernel<<<pN_BLK,pN_TH>>>(mycav,myDmesh);
						if (CudaTest("solving failed post_work_kernel"))
							exit(0);
					}
				#endif


				hactiveCavityCount = mycav.get_activeCavityCount();

				#if ENABLE_SECONDARY_MARKING
					myDmesh.init_possible();
					m_refine_checkCavity<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
					if (CudaTest("solving failed m_refine_checkCavity"))
						exit(0);
					hpossible = myDmesh.get_possible();

					while (hpossible) {
						m_refine_cavityMark2_atomic<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
						if (CudaTest("solving failed m_refine_cavityMark2_atomic"))
							exit(0);

						#ifdef ENABLE_POST_WORK
							mycav.init_post_bad_work_count();
						#endif
						m_refine_cavityReTriangulate<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
						if (CudaTest("solving failed m_refine_cavityReTriangulate"))
							exit(0);
						#ifdef ENABLE_POST_WORK
							unsigned post_wc=mycav.get_post_bad_work_count();
							if(post_wc){
								unsigned pN_TH = 128;
								unsigned pN_BLK = (post_wc + N_TH - 1) / (N_TH);
								post_work_kernel<<<pN_BLK,pN_TH>>>(mycav,myDmesh);
								if (CudaTest("solving failed post_work_kernel"))
									exit(0);
							}
						#endif

						myDmesh.init_possible();
						m_refine_checkCavity<<<N_BLK, N_TH>>>(myDmesh, mycav, Ccount);
						if (CudaTest("solving failed m_refine_checkCavity"))
							exit(0);
						hpossible = myDmesh.get_possible();

						hactiveCavityCount = mycav.get_activeCavityCount();
					}
				#endif
				cudaEventRecord(k_stop, 0);
				cudaEventSynchronize(k_stop);
				cudaEventElapsedTime(&t_time, k_start, k_stop);
				k_time += t_time;
			}
		}
#else
		ntriperit = ntriangles;
		ntriit = (ntriangles + ntriperit - 1) / ntriperit;
		for (unsigned ii = 0; ii < ntriit; ++ii) {

			cudaEvent_t k_start, k_stop;
			cudaEventCreate(&k_start);
			cudaEventCreate(&k_stop);
			cudaEventRecord(k_start, 0);

			drefine<<<R_NBLOCKS, R_BLOCKSIZE>>>(ii * ntriperit, (ii + 1) * ntriperit, myDmesh, mycav,iteration * 100, arrayin, arrayout);
			if (CudaTest("solving failed drefine"))
				exit(0);

			cudaEventRecord(k_stop, 0);
			cudaEventSynchronize(k_stop);
			cudaEventElapsedTime(&t_time, k_start, k_stop);
			k_time += t_time;
		}
#endif

		ntriangles = myDmesh.get_ntriangles();
		countbad<<<nblocks, blocksize>>>(myDmesh.isdel, myDmesh.isbad, ntriangles, nbad, 1000 + iteration, arrayin, arrayout, bcount);
		CudaTest("countbad failed");
		cudaMemcpy(&hnbad, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

		countnotdel<<<nblocks, blocksize>>>(myDmesh.isdel, ntriangles, nbad, 1001 + iteration, arrayin, arrayout, bcount);
		CudaTest("countdel failed");
		cudaMemcpy(&hndel, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

		printf("\n\n hndel=%d hnbad=%d orintriangles=%d ntriangles=%d", hndel, hnbad, orintriangles, ntriangles);

		hchanged = myDmesh.get_changed();
		if (hchanged && orintriangles == ntriangles) {
			if (R_NBLOCKS == 1 && R_BLOCKSIZE == 1) {
				printf("-------------------Termination-------------------\n");
				break;
			}
			R_NBLOCKS = 1;
			R_BLOCKSIZE = 1;
			cout << "------------------------------(1,1)\n" << endl;
		} else {
			R_NBLOCKS = nblocks;
			R_BLOCKSIZE = blocksize;
#if 0
			cout << "writing vtk output data file" << endl;
			char fname[100];
			static int ff=1;
			sprintf(fname,"odmr_%d.vtk",ff++);
			me.print2file_unstructuredgrid_isdel(fname, nodex, nodey,
					pnnodes, tnodes, pntriangles, isdel);
#endif
		}
		printf("\n");
	} while (hnbad); //(hchanged);
	endtime = clock();

	myresult.time = k_time;

	countbad<<<nblocks, blocksize>>>(myDmesh.isdel, myDmesh.isbad, ntriangles, nbad, 2000 + iteration, arrayin, arrayout, bcount);
	CudaTest("countbad failed");
	cudaMemcpy(&hnbad, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

	countnotdel<<<nblocks, blocksize>>>(myDmesh.isdel, ntriangles, nbad, 2001 + iteration, arrayin, arrayout, bcount);
	CudaTest("countdel failed");
	cudaMemcpy(&hndel, nbad, sizeof(unsigned), cudaMemcpyDeviceToHost);

	myresult.o_ntriangles = hndel;
	myresult.o_bad = hnbad;
	myresult.iteration = iteration;

	printf("NBAD=%d,NnonDEl=%d\n", hnbad, hndel);

#if VTKO
	cout << "writing vtk output data file" << endl;
	me.print2file_unstructuredgrid_isdel("odmr.vtk", nodex, nodey, pnnodes, tnodes, pntriangles, isdel);
#endif

	printf("verifying...\n");

	unsigned hnchanged = 0, *nchanged;
	if (cudaMalloc((void **) &nchanged, sizeof(unsigned)) != cudaSuccess)
		CudaTest("allocating nchanged failed");

	cudaMemcpy(nchanged, &hnchanged, sizeof(unsigned), cudaMemcpyHostToDevice);
	hchanged = false;
	cudaMemcpy(myDmesh.changed, &hchanged, sizeof(bool), cudaMemcpyHostToDevice);

	NBLOCKS = (ntriangles + _BS - 1) / _BS;
	dverify<<<NBLOCKS, _BS>>>(myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.isbad, myDmesh.isdel, nnodes, ntriangles, myDmesh.changed, nchanged);
	CudaTest("verification failed");

	cudaMemcpy(&hchanged, myDmesh.changed, sizeof(bool), cudaMemcpyDeviceToHost);
	cudaMemcpy(&hnchanged, nchanged, sizeof(unsigned), cudaMemcpyDeviceToHost);

	if (hchanged) {
		printf("verification failed: bad triangles exist: %d.\n", hnchanged);
	} else {
		printf("verification succeeded: 0 bad triangles exist.\n");
	}

	printf("iterations = %d.\n", iteration);
	runtime = (int) (1000.0f * (endtime - starttime) / CLOCKS_PER_SEC);
	printf("%d ms.\n", runtime);

	myresult.runtime = runtime;
	// cleanup left to the OS.
	writeResult();
	return 0;
}
