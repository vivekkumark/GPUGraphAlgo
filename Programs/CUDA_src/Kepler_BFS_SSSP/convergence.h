void printConv(int total_iteration,char *pattern){

//	char pattern[100]="Dist__L_1_it_%d_rmat20.gr_BFS.txt";
	cout<<"-------------Convergence Rate--------------"<<endl;
	char src_fname[100],it_fname[100];
	sprintf(src_fname,pattern,total_iteration);

	ofstream op;
	op.open("conv.vi",ios::app);
	op<<endl<<endl<<endl;
	op<<"----------------------------------------------------------------------------"<<endl;
	op << src_fname<<endl;
	op<<"----------------------------------------------------------------------------"<<endl;


	foru max_dist=0;
	ifstream ip_src,ip_it;
	ip_src.open(src_fname);
	foru fs,fi;

	unsigned nnum=0;
	while(!ip_src.eof()){
		ip_src >> fs;
		if(fs<MYINFINITY && fs>max_dist)
			max_dist=fs;
	}
	ip_src.close();

	cout<<"Max DIST:"<<max_dist<<endl;

	for(int it=0;it<total_iteration;it++){
		sprintf(it_fname,pattern,it);

		ip_src.open(src_fname);
		ip_it.open(it_fname);

		double norm2=0;
		nnum=0;
		unsigned diff_count=0;
		while(!ip_src.eof()){
			ip_src >> fs;
			ip_it  >> fi;
			diff_count+=(fs!=fi);
			norm2 +=(fs-fi) * (fs-fi);
			nnum++;
		}
		ip_src.close();
		ip_it.close();

		norm2=sqrt(norm2);

		cout<<"Iteration:"<< setw(3)<<left <<it<<" = " << norm2 <<"("<<diff_count<<")"<<endl;
		op<<"Iteration:"<< setw(3)<<left <<it<<" = " << norm2 <<"("<<diff_count<<")"<<endl;

		if(diff_count< 100){
			ip_src.open(src_fname);
			ip_it.open(it_fname);

			nnum=0;
			while(!ip_src.eof()){
				ip_src >> fs;
				ip_it  >> fi;
				if(fs != fi){
					cout <<"	vertex: "<< setw(9)<<right << nnum<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
					op <<"	vertex: "<< setw(9)<<right << nnum<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
				}
				nnum++;
			}

			ip_src.close();
			ip_it.close();
		}

		unlink(it_fname);
	}
	unlink(src_fname);
	op<<"----------------------------------------------------------------------------"<<endl;
	op.close();
}

void printConv_rtime_init(unsigned level,unsigned maxdeg,char *gnm,char *algo){
	ofstream op;
	op.open("conv.vi",ios::app);
	op << endl << endl << endl;
	op << "----------------------------------------------------------------------------------------";
	op << endl;
	op << "Dist_Level_" << level <<"_MaxDeg_"<< maxdeg <<"_" << gnm <<"_"<<algo;
	op << endl;
	op.close();
}
void printConv_rtime_final(unsigned totalIt,foru max_dist){
	ofstream op;
	op.open("conv.vi",ios::app);
	op << "Total Iteration : " << totalIt <<"     Max Dist : " << max_dist;
	op << endl;
	op << "----------------------------------------------------------------------------------------";
	op << endl;
	op.close();
}
void printConv_rtime(foru *hdist,Graph &hgraph,unsigned it){

	foru *final_dist;

	#ifdef SSSP
		final_dist=hgraph.sssp;
	#else
		final_dist=hgraph.bfs;
	#endif


	ofstream op;
	op.open("conv.vi",ios::app);

	unsigned ver=0;
	foru fi,fs;
	double norm2=0;
	unsigned diff_count=0;

	for(int k=0;k<hgraph.nnodes;k++){
		if(hgraph.IsParent[k]){
			fi = hdist[k];
			fs = final_dist[ver];
			norm2 += (fs-fi) * (fs-fi);
			diff_count+=(fs!=fi);
			ver++;
		}
	}
	norm2=sqrt(norm2);

	cout<<"Iteration:"<< setw(3)<<left <<it<<" = " << norm2 <<"("<<diff_count<<")"<<endl;
	op<<"Iteration:"<< setw(3)<<left <<it<<" = " << norm2 <<"("<<diff_count<<")"<<endl;


	if(diff_count< 100){
		ver=0;
		for(int k=0;k<hgraph.nnodes;k++){
			if(hgraph.IsParent[k]){
				fi = hdist[k];
				fs = final_dist[ver];
				if(fs != fi){
					cout <<" vertex: "<< setw(9)<<right << ver<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
					op <<" vertex: "<< setw(9)<<right << ver<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
					/*
					if(ver == hgraph.parentNode[k]){
						cout <<" vertex: "<< setw(9)<<right << ver<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
						op <<" vertex: "<< setw(9)<<right << ver<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
					}else{
						cout <<" vertex: "<< setw(9)<<right << ver<<"("<< hgraph.parentNode[k] <<")"<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
						op <<" vertex: "<< setw(9)<<right << ver<<"("<< hgraph.parentNode[k] <<")"<<"  "<<setw(9)<<right<<fi<<"->"<<setw(9)<<right<<fs<<endl;
					}
					*/
				}
				ver++;
			}
		}
	}

	op.close();
}


