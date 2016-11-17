#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"functions.h"

long long int check(float *arr, long long int row)
{
    while(--row>0 && fabs(arr[row]-arr[0]) > 0.00001);
    return row!=0;
}

/* square() is a function which squares the given number and returns the result. */
float square(float x)
{
    return (x*x);
}

/*
   compute_chi() is a function which will compute the chi value as per the formula
   chi = Summation(square(observed-expected)/expected)
*/
float compute_chi(float N,float *N_X,float *f_X)
{
    float res = 0;
    long long int i;
    for(i=0; i<no_groups; i++)
    {
        res += (square((N*f_X[i])-N_X[i]))/(N*f_X[i]);
    }
    return res;
}

/* find_consensus() function will determine the consensus sequence given a set of sequences. */
void find_consensus(char *con, char **sequences, long long int row, long long int column)
{
    /* counter : Keeps the count of each nucleotide ACGT in a column in that order.
       max     : Keeps the maximum count of the nucleotide in a particular column.
       index   : stores the index of the consensus nucleotide.
    */
    long long int counter[no_groups];
    long long int i,j,k,max = 0,index = 0, l;
    if(sequences == NULL)
    {
        con = NULL;
        return;
    }

    /* Count the ACGT occurance in each column. */
    for(i=0; i<column; i++)
    {
	for(l=0;l<no_groups;l++)
        	counter[l] = 0;
        for(j=0; j<row; j++)
        {
	    for(l=0;l<no_groups;l++)
            	if(sequences[j][i] == groups[l])
                	counter[l]++;
        }
        /* Find the consensus nucleotide based on the count. */
        max = 0;
        index = 0;
        for(k=0; k<no_groups; k++)
        {
            if(max < counter[k])
            {
                max = counter[k];
                index = k;
            }
        }
	con[i] = groups[index];
    }
    con[i] = '\0';
}

NODE compute_weights(NODE root, int col, int no_groups)
{
	extern long long int pseudo_n, pseudo_d; 
	float exp, denominator;
	float group_count[no_groups];
	long long int i, j, temp, k, l;
	char consensus = '\0';
        if(no_groups != 0)
        	exp = 1.0/no_groups;
	
	root->pwm = (float**)malloc(sizeof(float*) * no_groups);
	for(i=0;i<no_groups;i++)
	{
		root->pwm[i] = (float *) malloc(sizeof(float) * col);
	}
	
	for(j=0;j<col;j++)
	{
		for(k=0;k<no_groups;k++)
			group_count[k] = 0.0;
		
		temp = j;
		for(i=0;i<root->no_ele;i++)
		{
			if(root->position != -1 && root->left != NULL && root->right != NULL)
			{
				j = root->position;
			}
			consensus = root->seq[i][j];
			for(l=0;l<no_groups;l++)
				if(consensus == groups[l])
					group_count[l]++;
		}
		denominator = root->no_ele+pseudo_d;

		j = temp;
		for(l=0;l<no_groups;l++)
			root->pwm[l][j] = (log(((group_count[l]+pseudo_n)/(denominator))/exp))/log(2);
	}
	return root;
}

void create_pwm(NODE root, int no_groups)
{
	int col = 0;
	if(root == NULL)
		return;
	if(root->position != -1 && root->left != NULL && root->right != NULL)
		col = 1;
	else
		col = strlen(root->seq[0]);
	compute_weights(root, col, no_groups);
	create_pwm(root->left, no_groups);
	create_pwm(root->right, no_groups);
}
