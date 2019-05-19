__global__
void drefine(unsigned starttri, unsigned endtri, Dmesh myDmesh, Cavity myCavity,unsigned goal, volatile unsigned *arrayin, volatile unsigned *arrayout
		) {

	bool exitflag = false;

	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned nthreads = blockDim.x * gridDim.x;
	unsigned wpt = (endtri - starttri + nthreads - 1) / nthreads; //1;
	unsigned start = starttri + id * wpt;
	unsigned end = start + wpt;
	unsigned nnodes = *myDmesh.pnnodes;
	unsigned ntriangles = *myDmesh.pntriangles;

	unsigned centerelement = 0, row = 0;
	DIMSTYPE ceobtuse = 3, dims = 3;
	FORD centerx = 0.0, centery = 0.0;

	bool lchanged = false;

	unsigned *frontier = myCavity.frontier + id * CAVITY_SIZE;
	unsigned *pre = myCavity.pre + id * CAVITY_SIZE;
	unsigned *post = myCavity.post + id * CAVITY_SIZE;
	unsigned *cavity = myCavity.cavity + id * CAVITY_SIZE;
	unsigned *connections = myCavity.connections + id * 4 * CAVITY_SIZE;

	unsigned iifrontier,iipre,iipost,iicavity,iiconnections;

	for (unsigned tt = start; tt < end; ++tt) {

		__gpu_sync(++goal, arrayin, arrayout);
		__threadfence();
		__gpu_sync(++goal, arrayin, arrayout);

		if (tt < ntriangles && !myDmesh.isdel[tt] && myDmesh.isbad[tt]) {
			iicavity = iifrontier = iipre = iipost = iiconnections = 0;
			// cavity.initialize(tt);
			centerelement = tt;
			ceobtuse = myDmesh.obtuse[centerelement];

			unsigned itr = 0;
			while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) { // while it is obtuse.
				centerelement = getOpposite(centerelement, ceobtuse, myDmesh.neighbors, myDmesh.neighboredges, myDmesh.tnodes, nnodes, ntriangles);
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

					/*
					 if (next == INVALIDID) {
					 continue;
					 }*/

					if (myDmesh.isdel[next]) {
						printf("ce=%d,id=%d,curr=%d,next=%d,delAh\n", centerelement, id, curr, next);
						continue;
					}
					unsigned nextrow = 3 * next;
					unsigned nextdims = (myDmesh.tnodes[nextrow + 2] == INVALIDID ? 2 : 3);

					if (!(dims == 2 && nextdims == 2 && next != centerelement) && inCircumcircle(centerx, centery, next, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, nnodes, ntriangles)) {
						// isMember says next is part of the cavity, and we're not the second
						// segment encroaching on this cavity
						if (nextdims == 2 && dims != 2) {
							// is segment, and we are encroaching.
							iicavity = iifrontier = iipre = iipost = iiconnections = 0;
							centerelement = next;
							ceobtuse = myDmesh.obtuse[centerelement];
							itr = 0;
							while (ceobtuse < 3 && centerelement < ntriangles && ++itr < MAXITR) {
								centerelement = getOpposite(centerelement, ceobtuse, myDmesh.neighbors, myDmesh.neighboredges, myDmesh.tnodes, nnodes, ntriangles);
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
				myDmesh.owner[cavtri] = id;
			}
		}

		__gpu_sync(++goal, arrayin, arrayout);
		__threadfence();
		__gpu_sync(++goal, arrayin, arrayout);

		bool backoff = false;
		if (tt < ntriangles && !myDmesh.isdel[tt] && myDmesh.isbad[tt]) {
			// go over your triangles and see if they contain your id.
			for (unsigned ii = 0; ii < iicavity; ++ii) {
				unsigned cavtri = cavity[ii];
				if (myDmesh.owner[cavtri] < id) { // cavity overlap and the other thread has priority!
					backoff = true;
					break;
				} else if (myDmesh.owner[cavtri] > id) { // cavity overlap but you have the priority.
					myDmesh.owner[cavtri] = id; // mark it yours: due to this write, we require another checking phase.
				}
			}
		}

#ifdef ENABLE_POST_WORK
		*myCavity.post_bad_work_count=0;
#endif
		__gpu_sync(++goal, arrayin, arrayout);
		__threadfence();
		__gpu_sync(++goal, arrayin, arrayout);

		if (tt < ntriangles && !myDmesh.isdel[tt] && myDmesh.isbad[tt]) {
			// once again go over your triangles and see if they contain your id.
			if (!backoff) {
				for (unsigned ii = 0; ii < iicavity; ++ii) {
					unsigned cavtri = cavity[ii];
					if (myDmesh.owner[cavtri] != id) { // cavity overlap.
						backoff = true;
						break;
					}
				}
			}

			if (exitflag) {
				lchanged = true;
				exitflag = false;
				printf("id=%d,tt=%d,overflow\n", id, tt);
			}

			if (backoff) {
				lchanged = true;
			} else {
				// cavity.update(): create the new cavity based on the data of the old cavity.
				row = 3 * centerelement;
				dims = (myDmesh.tnodes[row + 2] == INVALIDID ? 2 : 3);

				unsigned newpoint = addPoint(centerx, centery, myDmesh.nodex, myDmesh.nodey, myDmesh.pnnodes, nnodes);
				if (dims == 2) { // we built around a segment.
					printf("dims=2\n");
					// create two segments (as triangles).
					unsigned newtriid1 = addTriangle(newpoint, myDmesh.tnodes[row + 0], INVALIDID, myDmesh.tnodes, myDmesh.pntriangles, ntriangles, myDmesh.isdel, myDmesh.obtuse, myDmesh.neighbors,
							myDmesh.neighboredges);
					unsigned newtriid2 = addTriangle(newpoint, myDmesh.tnodes[row + 1], INVALIDID, myDmesh.tnodes, myDmesh.pntriangles, ntriangles, myDmesh.isdel, myDmesh.obtuse, myDmesh.neighbors,
							myDmesh.neighboredges);
					// update triangles' neighbors: neighbors of the new triangles (segments) are the same as those of the previous segment?
					copyNeighbors(newtriid1, centerelement, myDmesh.neighbors, myDmesh.neighboredges);
					copyNeighbors(newtriid2, centerelement, myDmesh.neighbors, myDmesh.neighboredges);

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

					if (updateNeighbor(conndst, connsrc, newtri, myDmesh.neighbors, myDmesh.tnodes) == false) {
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


				// remove triangles from pre.
				for (unsigned ii = 0; ii < iipre; ++ii) {
					unsigned tri = pre[ii];
					myDmesh.isdel[tri] = true;
				}
				// add triangles from post, mark the bad triangles.
				// triangles are already added using addTriangle(), simply mark the bad ones.
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
			}
		}

#ifdef ENABLE_POST_WORK
		__gpu_sync(++goal, arrayin, arrayout);
		__threadfence();
		__gpu_sync(++goal, arrayin, arrayout);

		unsigned sum=*myCavity.post_bad_work_count;
		unsigned ppt=(sum+nthreads-1)/nthreads;

		for(int ii=0;ii<ppt && (ppt*id + ii) < sum;ii++) {
			unsigned tri = myCavity.post_bad_work[ppt*id+ii];
			myDmesh.isbad[tri] = checkbad(tri, myDmesh.nodex, myDmesh.nodey, myDmesh.tnodes, myDmesh.obtuse, ntriangles);
			lchanged |= myDmesh.isbad[tri];
		}
#endif

	}

	if (lchanged)
		*myDmesh.changed = true;

}
