int reval(char *seq, int start, int end) {

    int i,j;
    char temp;
    
    for (i = start, j = end; i < j; ++i, --j){
	temp = seq[i];
	seq[i] = seq[j];
	seq[j] = temp;
    }
    return 1;
}
