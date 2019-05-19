unsigned getminIndex(unsigned long *ary,unsigned len){
	unsigned minindex=0;
	unsigned long minval=ary[0];
	for(int k=1;k<len;k++){
		if(ary[k]<minval){
			minval=ary[k];
			minindex=k;
		}
	}
	return minindex;
}

unsigned getmaxIndex(unsigned long *ary,unsigned len){
	unsigned maxindex=0;
	unsigned long minval=ary[0];
	for(int k=1;k<len;k++){
		if(ary[k]>minval){
			minval=ary[k];
			maxindex=k;
		}
	}
	return maxindex;
}

void indexawarePrintindex(unsigned *ary,unsigned *aryindex,unsigned len,unsigned tindex,char *namepad){
	ofstream op;
	char name[100];
	sprintf(name,"Sm_index_%s.txt",namepad);
	op.open(name);


	unsigned *sm_count=(unsigned *)calloc(tindex,sizeof(unsigned));
	for(int j=0;j<len;j++)
		sm_count[aryindex[j]]++;

	for(int sm=0;sm<tindex;sm++){
		op<<sm_count[sm]<<",";
		int count=0;
		for(int i=0;i<len;i++){
			if(aryindex[i]==sm){
				op<<","<<i;
				count++;
			}
		}
		op<<endl;
	}
	op.close();
}


void indexawarePrintVal(unsigned *ary,unsigned *aryindex,unsigned len,unsigned tindex,char *namepad){
	ofstream op;
	char name[100];
	sprintf(name,"Sm_val_%s.txt",namepad);
	op.open(name);

	for(int sm=0;sm<tindex;sm++){
		for(int i=0;i<len;i++){
			if(aryindex[i]==sm){
				op<<ary[i]<<",";
			}
		}
		op<<endl;
	}
	op.close();
}

unsigned long combineClocks(unsigned *ary,unsigned len,unsigned *&aindex,char *namepad){
	unsigned SM=14;
	aindex=(unsigned *)malloc(len * sizeof(unsigned));
	unsigned long*mincount=(unsigned long*)calloc(SM,sizeof(unsigned long));

	for(int k=0;k<len;k++){
		unsigned index=getminIndex(mincount,SM);
		aindex[k]=index;
		mincount[index]+=ary[k];
	}

#ifdef PRINT_BSM
	indexawarePrintVal(ary,aindex,len,SM,namepad);
	indexawarePrintindex(ary,aindex,len,SM,namepad);
#endif
	return mincount[getmaxIndex(mincount,SM)];
}
