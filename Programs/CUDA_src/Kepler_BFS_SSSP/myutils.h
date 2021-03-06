void verifysolution(unsigned int *edgessrcdst, foru *edgessrcwt, unsigned int *noutgoing, unsigned int *psrc, foru *dist, unsigned int nv, unsigned int nedges, unsigned *srcsrc) {
	unsigned int ii, nn;
	unsigned int nerr = 0;
	for (nn = 0; nn < nv; ++nn) {
		unsigned int start = psrc[nn];
		unsigned int nsrcedges = noutgoing[nn];
		for (ii = 0; ii < nsrcedges; ++ii) {
			unsigned int u = nn;
			unsigned int v = edgessrcdst[start + ii];
			//float wt = 1;
			foru wt = BFS_SSSP(start + ii);
			if (wt > 0 && dist[u] + wt < dist[v]) {
				++nerr;
			}
		}
	}
	printf("verifying: no of errors = %d.\n", nerr);
}

template <class myType1>
inline void AryPrint(char name[],myType1 &var,uint32_t size){
	cout<<name<<endl;
	uint32_t i;
	for(i=0;i<size;i++)
		cout<<var[i]<<" ";
	cout<<endl;
}

template <class myType1>
inline void AryPrintFile(char name[],myType1 &var,uint32_t size){
	ofstream op;
	op.open(name);
	uint32_t i;
	for(i=0;i<size;i++)
		op<<var[i]<<endl;
	op.close();
	cout<<name<<endl;
}

inline void AryFile(char *fname,float *&ary,unsigned &indexaryLen){
	uint32_t i=0;
	FILE  *opf;
	opf=fopen(fname,"w+t");
	while(i < indexaryLen){
		fprintf(opf ,"%f\n", ary[i]);
		i++;
	}
	fclose(opf);
}

inline void Ary2File(char *fname,foru *&ary,unsigned *&indexary,unsigned &indexaryLen){
	uint32_t i=0;
	FILE  *opf;
	opf=fopen(fname,"w+t");
	while(i < indexaryLen){
		if(indexary[i]){
			fprintf(opf ,"%u\n", (unsigned)ary[i]);
		}
		i++;
	}
	fclose(opf);
}

void mystringRev(char *str){
	int len=strlen(str);
	char *s1=str;
	char *s2=str+len-1;
	while(s2>s1){
		char temp=*s1;
		*s1=*s2;
		*s2=temp;
		s1++;
		s2--;
	}
}
void getGraphname(char *ans,char *src){
	char val[256];
	strcpy(val,src);
	mystringRev(val);
	char *s2=strchr(val,'/');
	strncpy(ans,val,s2-val);
	ans[s2-val]=0;
	mystringRev(ans);
}

void writeResult(){

	FILE *fp;
	fp=fopen("result.vi","a+t");
	fprintf(fp,"%s,",myresult.exename);
	fprintf(fp,"%s,",myresult.algo);
	fprintf(fp,"%s,",myresult.graphfile);

	fprintf(fp,"%u,",myresult.bucketSize);

	fprintf(fp,"%u,",myresult.nnode);
	fprintf(fp,"%u,",myresult.nedge);
	fprintf(fp,"%u,",myresult.ns_nnode);
	fprintf(fp,"%u,",myresult.ns_nedge);

	fprintf(fp,"%u,",myresult.extranode);
	fprintf(fp,"%u,",myresult.maxdeg);
	fprintf(fp,"%u,",myresult.splitlevel);
	fprintf(fp,"%u,",myresult.MAX_EDGES_ALLOWED);
	fprintf(fp,"%u,",myresult.iterations);
	fprintf(fp,"%u,",myresult.runtime);
	fprintf(fp,"%f,",myresult.time);
	fprintf(fp,"%f,",myresult.ktime);
	fprintf(fp,"%f",myresult.sktime);

	fprintf(fp,"\n");
	fclose(fp);
}
