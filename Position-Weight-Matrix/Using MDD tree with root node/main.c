#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#include"functions.h"

int main(int argc, char *argv[])
{

    /* consensus: string to the store the consensus sequence.
       n 	    : Number of nucleotides in the consensus sequences (it is the column).
       pos_i    : column in the chi-squared table.
       pos_j    : row in the chi-squared table.
       row	    : row holds values of one row in a chi-squared table.
       max	    : holds the max of total of rows of the chi-squared table.
       max_index    : index which is used to fetch the conserved nucleotide.
       N	    : Total number of sequences in ci.
       N_X	    : Total number of sequences that have X in a particular position (say j).
       f_X	    : Relative frequency of nucloetide X in the sequences in ci_inv.
       pseudo_n : pseudocount for the numerator.
       pseudo_d : pseudocount for the denominator.
       ci_row   : Number of sequences in ci set.
       ci_inv_row: Number of sequences in ci_inv set.
       no_sequences: number of sequences in the input file which has training data.
    */
    char *consensus = NULL, *line = NULL, *k_mer = NULL;
    long long int *prev = NULL;
    long long int *stack = NULL, top = -1;

    float *row = NULL, N, *N_X = NULL, *f_X = NULL, *total_check = NULL;
    long long int *ci_row = NULL,*ci_inv_row = NULL;
    long long int no_sequences = 3000,n=20;
    int df_ind = -1;
    extern long long int pseudo_n, pseudo_d;
    extern int no_groups;

    if((!strcmp(argv[1],"-h")) || argc < 6)
    {
	printf("Usage: <executable> [window_size] <groups> <train_file_path> <test_file_path> <no_train_seq> <train_seq_len>\n");
	exit(0);
    }

    char lines[n], filename[250];

    long long int r = 0,ri = 0, ws = 0;
    long file_len = 0;

    long long int i,j,k;
    extern long long int limit;

    no_sequences = atoi(argv[argc-2]);
    n = atoi(argv[argc-1]);
    no_groups = strlen(argv[argc-5]);
    pseudo_n = 1;
    pseudo_d = no_groups;
    extern char *groups;
    groups = (char *)malloc(sizeof(char)*no_groups);
    strcpy(groups,argv[argc-5]);

    /* root: It is a root of a tree which has all the sequences from a file.
       cur : It represents current which holds the tree while calculating chi-squared table.
    */
    NODE root = NULL;

    /* open input file for reading the training data. */
    //FILE *fp = fopen("input.txt","r");
    FILE *fp = fopen(argv[argc-4],"r");
    FILE *fptr = NULL, *fout = NULL;

    for(i=0; !feof(fp); i++) 
	fscanf(fp,"%s",lines);

    no_sequences = i-1;
    limit = 0.15*no_sequences;
    fclose(fp);
    fp = fopen(argv[argc-4],"r");

    /* Get the root ready to store the training data. Read from the file and store it in root. */
    root = getnode();
    root->seq = (char**)malloc(sizeof(char *)*no_sequences);
    root->seq_no = (long long int*)malloc(sizeof(long long int)*no_sequences);
    if(root->seq == NULL)
        printf("Allocation failed\n");
    for(i=0; i<no_sequences; i++)
    {
        root->seq[i] = (char*)malloc(sizeof(char)*(n+1));
        if(root->seq[i] == NULL)
            printf("Allocation failed\n");
        fscanf(fp,"%s",root->seq[i]);
	root->seq_no[i] = i;
    }

    root->left = NULL;
    root->right = NULL;

    /* Now, i has the actual number of sequences which is stored in appropriate variable no_sequences. */
    
    root->no_ele = no_sequences;
    //printf("no_sequences=%lld\n",no_sequences);

    n = strlen(root->seq[0]);
    consensus = (char*)malloc(sizeof(char)*(n+1));

    prev = (long long int *)malloc(sizeof(long long int)*n);
    stack = (long long int *)malloc(sizeof(long long int)*100);

    /* Allocate memory for row, N_X and f_X. */
    row = (float*)malloc(sizeof(float)*n);
    N_X = (float*)malloc(sizeof(float)*no_groups);
    f_X = (float*)malloc(sizeof(float)*no_groups);
    total_check = (float*)malloc(sizeof(float)*n);

    if(n == 9)
	df_ind = 0;
    else
    {
	if(n == 14)
		df_ind = 1;
    }

    for(i=0; i<n; i++)
    {
        prev[i] = -1;
        total_check[i] = 0.0;
    }

    /* Allocate memories to ci_row and ci_inv_row. */
    ci_row = (long long int*)malloc(sizeof(long long int));
    ci_inv_row = (long long int*)malloc(sizeof(long long int));

    //formation(root, n, row, N_X, f_X, consensus, pseudo_n, pseudo_d, no_sequences, N, ci_row, ci_inv_row, r, ri, prev, total_check, stack, top, 0, df_ind);

    root->position = -1;
    create_pwm(root, no_groups);
    display(root);

    fptr = fopen(argv[argc-3],"r");
    //sprintf(filename,"output_scores");
    //fout = fopen(filename,"a");
    fout = fopen(argv[argc-8],"a");

    if(argc == 7)
    {
	char char_read;
	while(!feof(fptr))
	{
	    fscanf(fptr,"%c",&char_read);
	    if(isalpha(char_read))
	        file_len++;
	}
	fclose(fptr);
	fptr = fopen(argv[argc-3],"r");
	line = (char*)malloc(sizeof(char)*(file_len+1));
	ws = atoi(argv[argc-6]);
        k_mer = (char *)malloc(sizeof(char)*(ws+1));
        
	i=0;
	while(!feof(fptr))
	{
	    fscanf(fptr,"%c",&char_read);
	    if(isalpha(char_read))
		line[i++] = toupper(char_read);
	}
	line[i] = '\0';
	for(j=0;j<(strlen(line)-ws +1);j++){
	    for(i=0;i<ws;i++)
                k_mer[i] = line[j+i];
	    k_mer[i] = '\0';
            fprintf(fout,"%s\t%f\t%d\n",k_mer,traverse(root,k_mer),atoi(argv[argc-7]));
	}
    }
    else
    {
	line = (char*)malloc(sizeof(char)*(n+1));
	while(!feof(fptr))
        {
            if(fscanf(fptr,"%s",line) != 1) break;
	    i=0;
	    while(i<n)
	    {
		line[i] = toupper(line[i]);
		i++;
	    }
	    fprintf(fout,"%s\t%f\t%d\n",line,traverse(root,line),atoi(argv[argc-7]));
	}
    }

    fclose(fp);
    fclose(fptr);
    fclose(fout);

    free_ptr(row, N_X, f_X, consensus, ci_row, ci_inv_row, prev, total_check, stack, line, k_mer);
    freenode(root);

    return 0;
}
