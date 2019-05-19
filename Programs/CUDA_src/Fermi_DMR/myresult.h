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
void getname(char *ans,char *src){
	char val[256];
	strcpy(val,src);
	mystringRev(val);
	char *s2=strchr(val,'/');
	strncpy(ans,val,s2-val);
	ans[s2-val]=0;
	mystringRev(ans);
}


typedef struct Myresult{
	char fname[100];
	char ipname[100];
	unsigned i_ntriangles;
	unsigned i_bad;
	unsigned o_ntriangles;
	unsigned o_bad;
	unsigned iteration;
	float time;
	float runtime;
}Myresult;

Myresult myresult;

void writeResult(){

	FILE *fp;
	fp=fopen("result.vi","a+t");

	fprintf(fp,"%s,",myresult.fname);
	fprintf(fp,"%s,",myresult.ipname);
	fprintf(fp,"%u,",myresult.i_ntriangles);
	fprintf(fp,"%u,",myresult.i_bad);
	fprintf(fp,"%u,",myresult.o_ntriangles);
	fprintf(fp,"%u,",myresult.o_bad);
	fprintf(fp,"%u,",myresult.iteration);
	fprintf(fp,"%f,",myresult.runtime);
	fprintf(fp,"%f,",myresult.time);

	fprintf(fp,"\n");
	fclose(fp);
}
