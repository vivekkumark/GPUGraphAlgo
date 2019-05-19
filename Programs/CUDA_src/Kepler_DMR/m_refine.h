__global__
void post_work_kernel(Cavity myCavity,Dmesh myDmesh){

	unsigned wc = *myCavity.post_bad_work_count;
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned wpt = (wc + nthreads - 1) / nthreads;

	unsigned start = id * wpt;
	unsigned end = start + wpt;
	for (unsigned myid = start; myid < end && myid < wc; ++myid){
		unsigned tri = myCavity.post_bad_work[myid];
		unsigned ntriangles = *myDmesh.pntriangles;
		bool bad = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
		myDmesh.isbad[tri] = bad;
		if (bad)
			*myDmesh.changed = bad;
	}

}





__global__
void m_dynamicP_check_bad(unsigned iipost,unsigned *post, Dmesh myDmesh) {

	unsigned wpt=(iipost +blockDim.x -1)/ blockDim.x;

	unsigned start = threadIdx.x * wpt;
	unsigned end = start + wpt;

	for (unsigned ii = start; ii < end && ii < iipost; ii++) {
		unsigned tri = post[ii];
		unsigned ntriangles = *myDmesh.pntriangles;
		myDmesh.obtuse[tri] = 3;
		bool bad = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
		myDmesh.isbad[tri] = bad;
		if (bad)
			*myDmesh.changed = bad;
	}
}

__global__
void m_dynamicP_check_bad_updateNeigh(unsigned iipost, unsigned *post, unsigned myid, Dmesh myDmesh) {
	unsigned kk = threadIdx.x;
	if (kk < iipost) {
		unsigned nnodes = *myDmesh.pnnodes;
		unsigned ntriangles = *myDmesh.pntriangles;

		unsigned tri = post[kk];
		unsigned iineighbor = 1;
		for (unsigned jj = 0; jj < iipost && iineighbor < 3; ++jj) {
			if (kk != jj) {
				DIMSTYPE commonedgestart = adjacent(tri, post[jj], myDmesh.tnodes, nnodes, ntriangles);
				if (commonedgestart < 3) {
					unsigned newrow = tri * 3;
					myDmesh.neighbors[newrow + iineighbor] = post[jj];
					myDmesh.neighboredges[newrow + iineighbor] = commonedgestart;
					++iineighbor;
				}
			}
		}
		myDmesh.obtuse[tri] = 3;
		bool bad = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
		myDmesh.isbad[tri] = bad;
		if (bad)
			*myDmesh.changed = bad;
	}

}

__device__
void RebuildCavity(Cavity &myCavity, unsigned myid, Dmesh &myDmesh) {
	unsigned tt = myCavity.marker[myid];
	bool exitflag = false;
	unsigned centerelement = 0;
	DIMSTYPE ceobtuse = 3;
	FORD centerx = 0.0, centery = 0.0;

	unsigned iicavity, iifrontier, iipre, iipost, iiconnections = 0;
	unsigned nnodes = *myDmesh.pnnodes;
	unsigned ntriangles = *myDmesh.pntriangles;

	unsigned *frontier = myCavity.frontier + myid * CAVITY_SIZE;
	unsigned *pre = myCavity.pre + myid * CAVITY_SIZE;
	unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;
	unsigned *connections = myCavity.connections + myid * 4 * CAVITY_SIZE;

	iicavity = iifrontier = iipre = iipost = iiconnections = 0;
	// cavity.initialize(tt);
	centerelement = tt;
	ceobtuse = myDmesh.obtuse[centerelement];

	unsigned itr = 0;
	while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) { // while it is obtuse.
		centerelement = getOpposite(centerelement, ceobtuse, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges, myDmesh.tnodes, nnodes, ntriangles);
		if (centerelement < ntriangles) {
			ceobtuse = myDmesh.obtuse[centerelement];
		}
	}
	if (centerelement >= ntriangles || myDmesh.isdel[centerelement]) {
		centerelement = tt;
		ceobtuse = myDmesh.obtuse[centerelement];
	}
	getCenter(centerelement, centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, *myDmesh.pntriangles);

	pre[iipre++] = centerelement;
	frontier[iifrontier++] = centerelement;
	//DEBUGCHECK(iipre);

	// cavity.build();
	while (iifrontier > 0 && exitflag == false) {
		unsigned curr = frontier[--iifrontier];
		unsigned row = 3 * curr;
		DIMSTYPE dims = (myDmesh.tnodes[row + 2] == INVALIDID ? 2 : 3);
		unsigned lc = (myDmesh.tnodes[row + 2] == INVALIDID ? 1 : 3);
		for (unsigned ii = 0; ii < lc; ++ii) {
			//expand(curr, neighbors[row + ii]);
			unsigned next = myDmesh.neighbors[row + ii];
			if (myDmesh.isdel[next]) {
				printf("ce=%d,myid=%d,curr=%d,next=%d,delAh\n", centerelement, myid, curr, next);
				continue;
			}
			unsigned nextrow = 3 * next;
			unsigned nextdims = (myDmesh.tnodes[nextrow + 2] == INVALIDID ? 2 : 3);

			if (!(dims == 2 && nextdims == 2 && next != centerelement)
					&& inCircumcircle(centerx, centery, next, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, ntriangles)) {
				// isMember says next is part of the cavity, and we're not the second
				// segment encroaching on this cavity
				if (nextdims == 2 && dims != 2) {
					// is segment, and we are encroaching.
					iicavity = iifrontier = iipre = iipost = iiconnections = 0;
					centerelement = next;
					ceobtuse = myDmesh.obtuse[centerelement];
					itr = 0;
					while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) {
						centerelement = getOpposite(centerelement, ceobtuse, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges, myDmesh.tnodes, nnodes,
								ntriangles);
						if (centerelement < ntriangles) {
							ceobtuse = myDmesh.obtuse[centerelement];
						}
					}
					if (centerelement >= ntriangles || myDmesh.isdel[centerelement]) {
						centerelement = next;
						ceobtuse = myDmesh.obtuse[centerelement];
					}
					getCenter(centerelement, centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, *myDmesh.pntriangles);
					pre[iipre++] = centerelement;
					frontier[iifrontier++] = centerelement;
					//DEBUGCHECK(iipre);
				} else {
					unsigned jj;
					for (jj = 0; jj < iipre; ++jj) {
						if (pre[jj] == next) {
							break;
						}
					}
					if (jj == iipre) {
						pre[iipre++] = next;
						frontier[iifrontier++] = next;

						if (iipre >= CAVITY_SIZE) {
							exitflag = true;
							break;
						}
					}
					//DEBUGCHECK(iipre);
				}
			} else {
				// not a member
				// add the common edge between curr and next to connections if doesn't already exist.
				DIMSTYPE cestart = myDmesh.neighboredges[row + ii]; // see definition of next above.
				if (cestart >= 3) {
					continue;
				}
				unsigned connpt1 = myDmesh.tnodes[row + cestart];
				unsigned connpt2 = myDmesh.tnodes[row + (cestart + 1) % dims];

				unsigned jj;
				for (jj = 0; jj < iiconnections; jj += 4) {
					if (connections[jj] == connpt1 && connections[jj + 1] == connpt2) {
						break;
					}
				}
				if (jj == iiconnections) {
					connections[iiconnections++] = connpt1;
					connections[iiconnections++] = connpt2;
					connections[iiconnections++] = curr;
					connections[iiconnections++] = next;

					if (iiconnections >= CAVITY_SIZE * 4) {
						exitflag = true;
						break;
					}
					//DEBUGCHECK4(iiconnections);
				}
			}
		}
	}

	iicavity = 0;

	for (unsigned ii = 0; ii < iiconnections && exitflag == false; ii += 4) {
		CHECK_EXIST_ADD(cavity, iicavity, connections[ii + 2]);
		if (iicavity >= CAVITY_SIZE) {
			exitflag = true;
			break;
		}

		CHECK_EXIST_ADD(cavity, iicavity, connections[ii + 3]);
		if (iicavity >= CAVITY_SIZE) {
			exitflag = true;
			break;
		}
	}

	for (unsigned ii = 0; ii < iipre && exitflag == false; ++ii) {
		CHECK_EXIST_ADD(cavity, iicavity, pre[ii]);
		if (iicavity >= CAVITY_SIZE) {
			exitflag = true;
			break;
		}
	}

	for (unsigned ii = 0; ii < iicavity && exitflag == false; ++ii) {
		unsigned cavtri = cavity[ii];
		myDmesh.owner[cavtri] = MAXID;
	}

	myCavity.iicavity[myid] = iicavity;
	myCavity.iifrontier[myid] = iifrontier;
	myCavity.iipre[myid] = iipre;
	myCavity.iipost[myid] = iipost;
	myCavity.iiconnections[myid] = iiconnections;
	myCavity.centerelement[myid] = centerelement;
	myCavity.centerx[myid] = centerx;
	myCavity.centery[myid] = centery;
	myCavity.marker[myid] = tt;
	myCavity.cavityProcessed[myid] = false;

	if (exitflag) {
		printf("OverFlow:exitflag\n");
	}
}

__global__
void m_refine_cavityMark(unsigned starttri, unsigned endtri, Dmesh myDmesh, Cavity myCavity) {

	bool exitflag = false;
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;

	unsigned wpt = (endtri - starttri + nthreads - 1) / nthreads;
	unsigned start = starttri + id * wpt;
	unsigned end = start + wpt;

	unsigned nnodes = *myDmesh.pnnodes;
	unsigned ntriangles = *myDmesh.pntriangles;

	unsigned centerelement = 0;
	DIMSTYPE ceobtuse = 3;
	FORD centerx = 0.0, centery = 0.0;

	unsigned myid;
	unsigned iicavity, iifrontier, iipre, iipost, iiconnections = 0;

	for (unsigned tt = start; tt < end; ++tt) {

		if (tt < ntriangles && !myDmesh.isdel[tt] && myDmesh.isbad[tt]) {

			myid = atomicInc((unsigned *) myCavity.activeCavityCount, MAXID);

			//printf("myid=%d\n",myid);
			if(myid>=myCavity.MAXCAVITY){
				printf("Enough Over size\n");
				break;
			}

			unsigned *frontier = myCavity.frontier + myid * CAVITY_SIZE;
			unsigned *pre = myCavity.pre + myid * CAVITY_SIZE;
			unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;
			unsigned *connections = myCavity.connections + myid * 4 * CAVITY_SIZE;

			iicavity = iifrontier = iipre = iipost = iiconnections = 0;
			// cavity.initialize(tt);
			centerelement = tt;
			ceobtuse = myDmesh.obtuse[centerelement];

			unsigned itr = 0;
			while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) { // while it is obtuse.
				centerelement = getOpposite(centerelement, ceobtuse, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges, myDmesh.tnodes, nnodes, ntriangles);
				if (centerelement < ntriangles) {
					ceobtuse = myDmesh.obtuse[centerelement];
				}
			}
			if (centerelement >= ntriangles || myDmesh.isdel[centerelement]) {
				centerelement = tt;
				ceobtuse = myDmesh.obtuse[centerelement];
			}
			getCenter(centerelement, centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, *myDmesh.pntriangles);

			pre[iipre++] = centerelement;
			frontier[iifrontier++] = centerelement;
			//DEBUGCHECK(iipre);

			// cavity.build();
			while (iifrontier > 0 && exitflag == false) {
				unsigned curr = frontier[--iifrontier];
				unsigned row = 3 * curr;
				DIMSTYPE dims = (myDmesh.tnodes[row + 2] == INVALIDID ? 2 : 3);
				unsigned lc = (myDmesh.tnodes[row + 2] == INVALIDID ? 1 : 3);
				for (unsigned ii = 0; ii < lc; ++ii) {
					//expand(curr, neighbors[row + ii]);
					unsigned next = myDmesh.neighbors[row + ii];
					if (myDmesh.isdel[next]) {
						printf("ce=%d,id=%d,curr=%d,next=%d,delAh\n", centerelement, id, curr, next);
						continue;
					}
					unsigned nextrow = 3 * next;
					unsigned nextdims = (myDmesh.tnodes[nextrow + 2] == INVALIDID ? 2 : 3);

					if (!(dims == 2 && nextdims == 2 && next != centerelement)
							&& inCircumcircle(centerx, centery, next, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, ntriangles)) {
						// isMember says next is part of the cavity, and we're not the second
						// segment encroaching on this cavity
						if (nextdims == 2 && dims != 2) {
							// is segment, and we are encroaching.
							iicavity = iifrontier = iipre = iipost = iiconnections = 0;
							centerelement = next;
							ceobtuse = myDmesh.obtuse[centerelement];
							itr = 0;
							while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) {
								centerelement = getOpposite(centerelement, ceobtuse, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges, myDmesh.tnodes, nnodes,
										ntriangles);
								if (centerelement < ntriangles) {
									ceobtuse = myDmesh.obtuse[centerelement];
								}
							}
							if (centerelement >= ntriangles || myDmesh.isdel[centerelement]) {
								centerelement = next;
								ceobtuse = myDmesh.obtuse[centerelement];
							}
							getCenter(centerelement, centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, *myDmesh.pntriangles);
							pre[iipre++] = centerelement;
							frontier[iifrontier++] = centerelement;
							//DEBUGCHECK(iipre);
						} else {
							unsigned jj;
							for (jj = 0; jj < iipre; ++jj) {
								if (pre[jj] == next) {
									break;
								}
							}
							if (jj == iipre) {
								pre[iipre++] = next;
								frontier[iifrontier++] = next;

								if (iipre >= CAVITY_SIZE) {
									exitflag = true;
									break;
								}
							}
							//DEBUGCHECK(iipre);
						}
					} else {
						// not a member
						// add the common edge between curr and next to connections if doesn't already exist.
						DIMSTYPE cestart = myDmesh.neighboredges[row + ii]; // see definition of next above.
						if (cestart >= 3) {
							continue;
						}
						unsigned connpt1 = myDmesh.tnodes[row + cestart];
						unsigned connpt2 = myDmesh.tnodes[row + (cestart + 1) % dims];

						unsigned jj;
						for (jj = 0; jj < iiconnections; jj += 4) {
							if (connections[jj] == connpt1 && connections[jj + 1] == connpt2) {
								break;
							}
						}
						if (jj == iiconnections) {
							connections[iiconnections++] = connpt1;
							connections[iiconnections++] = connpt2;
							connections[iiconnections++] = curr;
							connections[iiconnections++] = next;

							if (iiconnections >= CAVITY_SIZE * 4) {
								exitflag = true;
								break;
							}
							//DEBUGCHECK4(iiconnections);
						}
					}
				}
			}

			iicavity = 0;

			for (unsigned ii = 0; ii < iiconnections && exitflag == false; ii += 4) {
				CHECK_EXIST_ADD(cavity, iicavity, connections[ii + 2]);
				if (iicavity >= CAVITY_SIZE) {
					exitflag = true;
					break;
				}

				CHECK_EXIST_ADD(cavity, iicavity, connections[ii + 3]);
				if (iicavity >= CAVITY_SIZE) {
					exitflag = true;
					break;
				}
			}

			for (unsigned ii = 0; ii < iipre && exitflag == false; ++ii) {
				CHECK_EXIST_ADD(cavity, iicavity, pre[ii]);
				if (iicavity >= CAVITY_SIZE) {
					exitflag = true;
					break;
				}
			}

			// mark the triangles in the cavity.
			for (unsigned ii = 0; ii < iicavity && exitflag == false; ++ii) {
				unsigned cavtri = cavity[ii];
#ifdef ENABLE_ATOMIC_MARKING
				atomicMin(&myDmesh.owner[cavtri] , tt);
#else
				myDmesh.owner[cavtri] = tt;
#endif
			}

			if (iicavity == 0)
				printf("-------------------------------------------iicavity zero\n");

			myCavity.iicavity[myid] = iicavity;
			myCavity.iifrontier[myid] = iifrontier;
			myCavity.iipre[myid] = iipre;
			myCavity.iipost[myid] = iipost;
			myCavity.iiconnections[myid] = iiconnections;
			myCavity.centerelement[myid] = centerelement;
			myCavity.centerx[myid] = centerx;
			myCavity.centery[myid] = centery;
			myCavity.marker[myid] = tt;
			myCavity.cavityProcessed[myid] = false;

			if (nthreads == 1)
				break;

			if (exitflag) {
				printf("OverFlow:exitflag\n");
			}
		}
	}
}

__global__
void m_refine_cavityReMark(Dmesh myDmesh, Cavity myCavity, unsigned Ccount) {

	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned ncav = Ccount;
	unsigned wpt = (ncav + nthreads - 1) / nthreads;
	unsigned start = id * wpt;
	unsigned end = start + wpt;

	for (unsigned myid = start; myid < end && myid < ncav; ++myid) {
		// go over your triangles and see if they contain your id.
		if (myCavity.cavityProcessed[myid]) {
			continue;
		}
		unsigned iicavity = myCavity.iicavity[myid];
		unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;
		unsigned marker = myCavity.marker[myid];
		for (unsigned ii = 0; ii < iicavity; ++ii) {
			unsigned cavtri = cavity[ii];
			if (myDmesh.owner[cavtri] < marker) { // cavity overlap and the other thread has priority!
				break;
			} else if (myDmesh.owner[cavtri] > marker) { // cavity overlap but you have the priority.
				myDmesh.owner[cavtri] = marker; // mark it yours: due to this write, we require another checking phase.
			}
		}
	}
}

__global__
void m_refine_cavityReTriangulate(Dmesh myDmesh, Cavity myCavity, unsigned Ccount) {

	unsigned ncav = Ccount;
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned wpt = (ncav + nthreads - 1) / nthreads;

	//printf("wpt=%d\n",wpt);

	unsigned start = id * wpt;
	unsigned end = start + wpt;
	unsigned nnodes = *myDmesh.pnnodes;
	unsigned ntriangles = *myDmesh.pntriangles;

	unsigned row = 0;
	DIMSTYPE dims = 3;
	FORD centerx = 0.0, centery = 0.0;

	bool lchanged = false;

	unsigned iicavity, iipre, iipost, iiconnections;

	for (unsigned myid = start; myid < end && myid < ncav; ++myid) {
		if (myCavity.cavityProcessed[myid]) {
			continue;
		}
		iicavity = myCavity.iicavity[myid];
		iipre = myCavity.iipre[myid];
		iipost = myCavity.iipost[myid];
		iiconnections = myCavity.iiconnections[myid];
		unsigned marker = myCavity.marker[myid];

		unsigned *pre = myCavity.pre + myid * CAVITY_SIZE;
		unsigned *post = myCavity.post + myid * CAVITY_SIZE;
		unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;
		unsigned *connections = myCavity.connections + myid * 4 * CAVITY_SIZE;

		bool backoff = false;
		for (unsigned ii = 0; ii < iicavity; ++ii) {
			unsigned cavtri = cavity[ii];
			if (myDmesh.owner[cavtri] != marker) { // cavity overlap.
				backoff = true;
				break;
			}
		}
		if (backoff) {
			lchanged = true;
		} else {

			//atomicInc((unsigned *)go, MAXID);
			//atomicDec((unsigned *)myCavity.activeCavityCount, MAXID);

			unsigned centerelement = myCavity.centerelement[myid];
			myCavity.cavityProcessed[myid] = true;
			centerx = myCavity.centerx[myid];
			centery = myCavity.centery[myid];

//			if(!iicavity)
//				printf("myid=%d,iiconnections=%d,iicavity=%d,ce=%d\n",myid,iiconnections,iicavity,centerelement);

			// cavity.update(): create the new cavity based on the data of the old cavity.
			row = 3 * centerelement;
			dims = (myDmesh.tnodes[row + 2] == INVALIDID ? 2 : 3);

			unsigned newpoint = addPoint(centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.pnnodes, nnodes);
			if (dims == 2) { // we built around a segment.
				printf("dims=2\n");
				// create two segments (as triangles).
				unsigned newtriid1 = addTriangle(newpoint, myDmesh.tnodes[row + 0], INVALIDID, myDmesh.tnodes, myDmesh.pntriangles, ntriangles, myDmesh.isdel, myDmesh.obtuse,
						(unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges);
				unsigned newtriid2 = addTriangle(newpoint, myDmesh.tnodes[row + 1], INVALIDID, myDmesh.tnodes, myDmesh.pntriangles, ntriangles, myDmesh.isdel, myDmesh.obtuse,
						(unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges);
				// update triangles' neighbors: neighbors of the new triangles (segments) are the same as those of the previous segment?
				copyNeighbors(newtriid1, centerelement, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges);
				copyNeighbors(newtriid2, centerelement, (unsigned *) myDmesh.neighbors, (unsigned *) myDmesh.neighboredges);

				post[iipost++] = newtriid1;
				post[iipost++] = newtriid2;
				//DEBUGCHECK(iipost);
			}

			bool postflag = false;

			for (unsigned ii = 0; ii < iiconnections; ii += 4) {
				unsigned connpt1 = connections[ii + 0];
				unsigned connpt2 = connections[ii + 1];
				unsigned connsrc = connections[ii + 2];
				unsigned conndst = connections[ii + 3];

				unsigned newtri = atomicInc(myDmesh.pntriangles, MAXID);

				#ifdef ENABLE_POST_WORK
					unsigned pid=atomicInc((unsigned *)myCavity.post_bad_work_count,MAXID);
					myCavity.post_bad_work[pid]=newtri;
				#endif


				unsigned newrow = 3 * newtri;
				myDmesh.tnodes[newrow + 0] = newpoint;
				myDmesh.tnodes[newrow + 1] = connpt1;
				myDmesh.tnodes[newrow + 2] = connpt2;

				myDmesh.isdel[newtri] = false;
				myDmesh.obtuse[newtri] = 3;
				ntriangles = newtri; // update.

				//if pre contains conndst
				unsigned jj;
				for (jj = 0; jj < iipre; ++jj) {
					if (pre[jj] == conndst) {
						break;
					}
				}
				//unsigned newconn = (jj == iipre ? conndst : connsrc);
				if (jj != iipre)
					printf("poda\n");

				// newtri and newconn are triangles, and their common edge is (connpt1, connpt2).
				// thus they are adjacent; their neighbors need to be updated.
				//unsigned newrow = 3 * newtri;
				myDmesh.neighbors[newrow] = conndst;
				myDmesh.neighboredges[newrow] = 1; // since connpt1 is point1 (newpoint is point0).

				if (updateNeighbor(conndst, connsrc, newtri, (unsigned *) myDmesh.neighbors, myDmesh.tnodes) == false) {
					printf("ce=%d,id=%d,connsrc=%d,conndst=%d,newtri=%d\n", centerelement, id, connsrc, conndst, newtri);
					PRINT_TRI(connsrc);
					PRINT_TRI(conndst);
					PRINT_TRI(newtri);
					postflag = true;
				}

				if (iipost < CAVITY_SIZE) {
					post[iipost++] = newtri;
				} else {
					printf("problem\n");
				}
			}

			if (iipost >= CAVITY_SIZE) {
				printf("problem\n");
			}

			// remove triangles from pre.
			for (unsigned ii = 0; ii < iipre; ++ii) {
				unsigned tri = pre[ii];
				myDmesh.isdel[tri] = true;
			}


			for (unsigned kk = 0; kk < iipost; ++kk) {
				unsigned iineighbor = 1;
				for (unsigned jj = 0; jj < iipost && iineighbor < 3; ++jj) {
					if (kk != jj) {
						DIMSTYPE commonedgestart = adjacent(post[kk], post[jj], myDmesh.tnodes, nnodes, ntriangles);
						if (commonedgestart < 3) {
							unsigned newrow = post[kk] * 3;
							myDmesh.neighbors[newrow + iineighbor] = post[jj];
							myDmesh.neighboredges[newrow + iineighbor] = commonedgestart;
							++iineighbor;
						}
					}
				}
			}

			#ifdef ENABLE_DYNAMICP
			//printf("iipost=%d\n",iipost);
				cudaStream_t kernel_stream;
				cudaStreamCreateWithFlags(&kernel_stream, cudaStreamNonBlocking);
				//m_dynamicP_check_bad_updateNeigh<<<1,iipost>>>(iipost, post, myid, myDmesh);
				m_dynamicP_check_bad<<<1,iipost>>>(iipost,post,myDmesh);
				cudaError_t e = cudaGetLastError();
				if (cudaSuccess != e) {
					//printf("Child Kernel Problem:code=%d\n",e);
					// add triangles from post, mark the bad triangles.
					// triangles are already added using addTriangle(), simply mark the bad ones.

					//If dynamic parallelism fails then current thread will take care of
					//checking bad triangles
					for (unsigned ii = 0; ii < iipost; ++ii) {
						unsigned tri = post[ii];
						#ifdef ENABLE_DYNAMICP_DMR //update parent cavity
							cavity[iicavity++] = tri;
						#endif
						if (postflag) {
							printf("id=%d,post tri=%d\n", id, tri);
						}
						myDmesh.obtuse[tri] = 3;
						myDmesh.isbad[tri] = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
						lchanged |= myDmesh.isbad[tri];
					}
				}
				cudaStreamDestroy(kernel_stream);
			#else
			// add triangles from post, mark the bad triangles.
			// triangles are already added using addTriangle(), simply mark the bad ones.
				//If post work is enabled then seperate kernel launch will take care of
				//checking bad triangles
				#ifndef ENABLE_POST_WORK
				for (unsigned ii = 0; ii < iipost; ++ii) {
					unsigned tri = post[ii];
					if (postflag) {
						printf("id=%d,post tri=%d\n", id, tri);
					}
					myDmesh.obtuse[tri] = 3;
					myDmesh.isbad[tri] = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
					lchanged |= myDmesh.isbad[tri];
				}
				#endif
			#endif
		}
	}

	if (lchanged)
		*myDmesh.changed = true;
}

__global__
void m_refine_checkCavity(Dmesh myDmesh, Cavity myCavity, unsigned Ccount) {
	unsigned myid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned ncav = Ccount;
	if (myid < ncav) {
		if (myCavity.cavityProcessed[myid] == false) {

			unsigned iicavity = myCavity.iicavity[myid];
			unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;

			if (iicavity == 0) {
				myCavity.cavityProcessed[myid] = true;
				//atomicDec((unsigned *) myCavity.activeCavityCount, MAXID);
				return;
			}

			if (myDmesh.isdel[myCavity.marker[myid]]) {
				//atomicDec((unsigned *)myCavity.activeCavityCount, MAXID);
				myCavity.cavityProcessed[myid] = true;
			} else {
				unsigned ii;
				for (ii = 0; ii < iicavity; ii++) {
					unsigned cavtri = cavity[ii];
					myDmesh.owner[cavtri] = MAXID;
					if (myDmesh.isdel[cavtri]) {
						break;
					}
				}
				if (ii < iicavity) {
					RebuildCavity(myCavity, myid, myDmesh);
				}
				*myDmesh.possible = true;
			}
		}
	}
}

__global__
void m_refine_cavityMark2_atomic(Dmesh myDmesh, Cavity myCavity, unsigned Ccount) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned ncav = Ccount;
	unsigned wpt = (ncav + nthreads - 1) / nthreads;
	unsigned start = id * wpt;
	unsigned end = start + wpt;

	for (unsigned myid = start; myid < end && myid < ncav; ++myid) {
		// go over your triangles and see if they contain your id.
		if (myCavity.cavityProcessed[myid]) {
			continue;
		}
		unsigned iicavity = myCavity.iicavity[myid];
		unsigned *cavity = myCavity.cavity + myid * CAVITY_SIZE;
		unsigned marker = myCavity.marker[myid];
		for (unsigned ii = 0; ii < iicavity; ++ii) {
			unsigned cavtri = cavity[ii];
			atomicMin(&myDmesh.owner[cavtri], marker);
		}
	}
}
