#include <stdio.h>

int testfile(char *file) {

	FILE *f;
	int flag;

	if((f=fopen(file,"r"))==NULL) {
		flag = 0;
	} else  {
		flag = 1;
		fclose(f);
	}
	return flag;
}
