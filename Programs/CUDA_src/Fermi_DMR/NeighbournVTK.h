#define FORD float
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <cuda.h>

using namespace std;
#define INVALIDID	1234567890


class vivek {
public:
	bool FileExist(char *fname){
		struct stat st;
		if(stat(fname,&st) == 0)
		    return true;
		else
			return false;
	}

	bool findneigh_file_check(char *bfname,unsigned *&neighboredges,unsigned *&neighbors,unsigned ntriangles){
		char fname[255];
		strcpy(fname,bfname);
		strcat(fname,".neigh");

		//unlink(fname);
		if(FileExist(fname)){
			cout<<"reading neighbor file:"<<fname<<endl;
			unsigned *h_neighbors=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));
			unsigned *h_neighboredges=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));

			ifstream fp;
			fp.open(fname);
			for(unsigned i=0;i<ntriangles;i++){
				int row=i*3;
				fp >> h_neighboredges[row+0] ;
				fp >> h_neighboredges[row+1] ;
				fp >> h_neighboredges[row+2] ;

				fp >> h_neighbors[row+0] ;
				fp >> h_neighbors[row+1] ;
				fp >> h_neighbors[row+2] ;

			}
			fp.close();

			cudaMemcpy(neighbors, h_neighbors, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);
			cudaMemcpy(neighboredges, h_neighboredges, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);

			free(h_neighbors);
			free(h_neighboredges);

			return true;
		}
		cout<<"Neigh file not found"<<endl;
		return false;
	}

	void write_neigh_to_file(char *bfname,unsigned *&neighboredges,unsigned *&neighbors,unsigned ntriangles){
		char fname[255];
		strcpy(fname,bfname);
		strcat(fname,".neigh");

		cout<<"Writing neighbor file:"<<fname<<endl;
		unsigned *h_neighbors=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));
		unsigned *h_neighboredges=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));

		cudaMemcpy(h_neighbors, neighbors, ntriangles*3*sizeof(unsigned) , cudaMemcpyDeviceToHost);
		cudaMemcpy(h_neighboredges, neighboredges, ntriangles*3*sizeof(unsigned) , cudaMemcpyDeviceToHost);


		ofstream fp;
		fp.open(fname);
		for(unsigned i=0;i<ntriangles;i++){
			int row=i*3;
			fp << h_neighboredges[row+0] << " ";
			fp << h_neighboredges[row+1] << " ";
			fp << h_neighboredges[row+2] << " ";

			fp << h_neighbors[row+0] << " ";
			fp << h_neighbors[row+1] << " ";
			fp << h_neighbors[row+2] << " ";

			fp << endl;
		}
		fp.close();

		free(h_neighbors);
		free(h_neighboredges);
	}


	void findneighbor_serial(char *bfname,unsigned *&neighboredges,unsigned *&neighbors,unsigned *&htnodes,unsigned ntriangles){
		char fname[255];
		strcpy(fname,bfname);
		strcat(fname,".neigh");

		//unlink(fname);

		if(FileExist(fname)){
			cout<<"reading neighbor file:"<<fname<<endl;
			unsigned *h_neighbors=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));
			unsigned *h_neighboredges=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));

			ifstream fp;
			fp.open(fname);
			for(unsigned i=0;i<ntriangles;i++){
				int row=i*3;
				fp >> h_neighboredges[row+0] ;
				fp >> h_neighboredges[row+1] ;
				fp >> h_neighboredges[row+2] ;

				fp >> h_neighbors[row+0] ;
				fp >> h_neighbors[row+1] ;
				fp >> h_neighbors[row+2] ;

			}
			fp.close();

			cudaMemcpy((void *)neighbors, h_neighbors, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);
			cudaMemcpy((void *)neighboredges, h_neighboredges, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);
		}
		else{
			cout<<"serial finding neighbor and writing into file:"<<fname<<endl;
			unsigned *h_neighbors=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));
			unsigned *h_neighboredges=(unsigned *)malloc(ntriangles*3*sizeof(unsigned));

			cout<<ntriangles<<endl;

			unsigned n_find=ntriangles;

			for(unsigned trione=0;trione<n_find;trione++){
				printf("finding neighbors: %8d ,  %3d%% complete.\r",trione , (int)(trione*100.0 / ntriangles));
				int iirow=0;
				unsigned rowone = 3 * trione;
				for(unsigned tritwo=0;tritwo<ntriangles;tritwo++){
					if(trione != tritwo){

						unsigned rowtwo = 3 * tritwo;

						unsigned dimsone = (htnodes[rowone + 2] == INVALIDID ? 2 : 3);
						unsigned dimstwo = (htnodes[rowtwo + 2] == INVALIDID ? 2 : 3);

						bool flag=true;
						int commonedge=3;

						for (unsigned ii = 0; ii < dimsone && flag ; ++ii){
							unsigned pt1=htnodes[rowone+ii];
							unsigned pt2=htnodes[rowone+(ii+1)%dimsone];

							for(unsigned jj=0;jj<dimstwo && flag;jj++){
								unsigned t_pt1=htnodes[rowtwo+jj];
								unsigned t_pt2=htnodes[rowtwo+(jj+1)%dimstwo];

								if((pt1==t_pt1 && pt2==t_pt2) || (pt1==t_pt2 && pt2==t_pt1)){
									commonedge=ii;
									flag=false;
									//printf("\npt1=%d,pt2=%d  t_pt1=%d,t_pt2=%d",pt1,pt2,t_pt1,t_pt2);
								}
							}

						}

						if (flag==false && iirow<3){
							h_neighbors[rowone + iirow]=tritwo;
							h_neighboredges[rowone + iirow]=commonedge;
							iirow++;
						}

						if(iirow>=3){
							break;
						}
					}

				}
				for (; iirow < 3; ++iirow) {
					h_neighbors[rowone + iirow] = INVALIDID;
					h_neighboredges[rowone + iirow] = 3;
				}
			}

			ofstream fp;
			fp.open(fname);
			for(unsigned i=0;i<n_find;i++){
				int row=i*3;
				fp << h_neighboredges[row+0] << " ";
				fp << h_neighboredges[row+1] << " ";
				fp << h_neighboredges[row+2] << " ";

				fp << h_neighbors[row+0] << " ";
				fp << h_neighbors[row+1] << " ";
				fp << h_neighbors[row+2] << " ";

				fp << endl;
			}
			fp.close();

			cudaMemcpy((void *)neighbors, h_neighbors, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);
			cudaMemcpy((void *)neighboredges, h_neighboredges, ntriangles*3*sizeof(unsigned) , cudaMemcpyHostToDevice);

			free(h_neighboredges);
			free(h_neighbors);

		}

	}

	void print2file_unstructuredgrid(char *ofname,FORD *nodex,FORD *nodey,unsigned npt,unsigned *triangle,unsigned ntri){

		unsigned int i;
		ofstream fp;
		fp.open(ofname);
		fp << "# vtk DataFile Version 4.2" << endl;
		fp << "title:Unstructuredgrid data" << endl<<endl;
		fp << "ASCII" << endl;
		fp << "DATASET UNSTRUCTURED_GRID" << endl;

		fp << "POINTS "<< npt << " float"<<endl;
		for(i=0;i<npt;i++){
			fp << nodex[i] << " ";
			fp << nodey[i] << " ";
			fp << 0 << endl;
		}


		unsigned tt=ntri;



		int count_tri=0,count_poly=0;

		for(i=0;i<tt;i++){
			int row=i*3;
			if(triangle[row+2]==INVALIDID){
				count_poly++;
				//cout << "("<<triangle[row+0] <<","<<triangle[row+1] <<","<<triangle[row+2] <<")"<<endl;
			}
			else{
				count_tri++;
			}
		}


		cout<<"count_poly : " <<count_poly <<endl;
		cout<<"count_tri : " <<count_tri <<endl;

		fp << "CELLS " << tt << " " << count_poly * 3 + count_tri *4  << endl;
		for(i=0;i<tt;i++){
			int row=i*3;
			if(triangle[row+2]==INVALIDID){
				fp << 2 << " ";
				fp << triangle[row+0] << " ";
				fp << triangle[row+1] << endl;
			}
			else{
				fp << 3 << " ";
				fp << triangle[row+0] << " ";
				fp << triangle[row+1] << " ";
				fp << triangle[row+2] << endl;
			}
		}

		fp << "CELL_TYPES " << tt <<endl;
		for(i=0;i<tt;i++){
			int row=i*3;
			if(triangle[row+2]==INVALIDID){
				fp << 3 <<endl;
			}else{
				fp << 5 <<endl;
			}
		}

		fp.close();
	}


	void print2file_unstructuredgrid_isdel(char *ofname,FORD *dnodex,FORD *dnodey,unsigned *dnpt,unsigned *dtriangle,unsigned *dntri,bool *disdel){
			unsigned npt;
			cudaMemcpy(&npt, dnpt, sizeof(unsigned), cudaMemcpyDeviceToHost);

			unsigned ntri;
			cudaMemcpy(&ntri, dntri, sizeof(unsigned), cudaMemcpyDeviceToHost);

			FORD *nodex=(FORD *)malloc(npt * sizeof(FORD));
			cudaMemcpy(nodex, dnodex, sizeof(FORD) * npt, cudaMemcpyDeviceToHost);

			FORD *nodey=(FORD *)malloc(npt * sizeof(FORD));
			cudaMemcpy(nodey, dnodey, sizeof(FORD) * npt, cudaMemcpyDeviceToHost);

			unsigned *triangle=(unsigned *)malloc(ntri * sizeof(unsigned) *3);
			cudaMemcpy(triangle, dtriangle, sizeof(unsigned) * 3 * ntri, cudaMemcpyDeviceToHost);

			bool *isdel=(bool *)malloc(ntri * sizeof(bool));
			cudaMemcpy(isdel, disdel, sizeof(bool) *  ntri, cudaMemcpyDeviceToHost);

			unsigned int i;
			ofstream fp;
			fp.open(ofname);

			fp << "# vtk DataFile Version 4.2" << endl;
			fp << "title:Unstructuredgrid data" << endl<<endl;
			fp << "ASCII" << endl;
			fp << "DATASET UNSTRUCTURED_GRID" << endl;

			fp << "POINTS "<< npt << " float"<<endl;

			//#define isnan(a) (a != a)

			for(i=0;i<npt;i++){
				if(isnan(nodex[i])){
					fp << 0 << " ";
				}else{
					fp << nodex[i] << " ";
				}

				if(isnan(nodey[i])){
					fp << 0 << " ";
				}else{
					fp << nodey[i] << " ";
				}

				fp << 0 << endl;
			}


			unsigned tt=ntri;


			int count_tri=0,count_poly=0;

			for(i=0;i<tt;i++){
				if(isdel[i]==false){
					int row=i*3;
					if(triangle[row+2]==INVALIDID){
						count_poly++;
					}
					else{
						count_tri++;
					}
				}
			}

			cout << "tri : " <<ntri<<endl;
			cout<<"count_poly : " <<count_poly <<endl;
			cout<<"count_tri : " <<count_tri <<endl;


			fp << "CELLS " << count_poly + count_tri << " " << count_poly * 3 + count_tri *4  << endl;

			for(i=0;i<tt;i++){
				if(isdel[i]==false){
					int row=i*3;
					if(triangle[row+2]==INVALIDID){
						fp << 2 << " ";
						fp << triangle[row+0] << " ";
						fp << triangle[row+1] << endl;
					}
					else{
						fp << 3 << " ";
						fp << triangle[row+0] << " ";
						fp << triangle[row+1] << " ";
						fp << triangle[row+2] << endl;
					}
				}

			}

			fp << "CELL_TYPES " << count_poly + count_tri <<endl;

			for(i=0;i<tt;i++){
				if(isdel[i]==false){
					int row=i*3;
					if(triangle[row+2]==INVALIDID){
						fp << 3 <<endl;
					}else{
						fp << 5 <<endl;
					}
				}
			}

			fp.close();
		}
};

