#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"functions.h"

/* freenode() will free the memory of a node. */
void freenode(NODE set)
{
    long long int i;

    if(set == NULL) return;
    freenode(set->left);
    freenode(set->right);
    for(i=0; i<set->no_ele; i++)
    {
        free(set->seq[i]);
	set->seq[i] = NULL;
    }
    for(i=0;i<4;i++)
    {
    	free(set->pwm[i]);
	set->pwm[i] = NULL;
    }
    free(set->pwm);
    set->pwm = NULL;
    free(set->seq);
    set->seq = NULL;
    free(set);
    set = NULL;
}

void free_children(NODE cur)
{
	long long int i;
	if(cur->left != NULL)	
	{
		if(cur->left->seq != NULL)
        	{
            		for(i=0; i<cur->left->no_ele; i++)
            		{
                		if(cur->left->seq[i] != NULL)
                    			free(cur->left->seq[i]);
            		}
            		if(cur->left->seq != NULL)
                		free(cur->left->seq);
        	}
		if(cur->left->seq_no != NULL)
		{
			free(cur->left->seq_no);
			cur->left->seq_no = NULL;
		}
        	if(cur->left != NULL)
		{
            		free(cur->left);
			cur->left = NULL;
		}
	}
	if(cur->right != NULL)	
	{
		if(cur->right->seq != NULL)
        	{
            		for(i=0; i<cur->right->no_ele; i++)
            		{
                		if(cur->right->seq[i] != NULL)
                    			free(cur->right->seq[i]);
            		}
            		if(cur->right->seq != NULL)
                		free(cur->right->seq);
        	}
		if(cur->right->seq_no != NULL)
		{
			free(cur->right->seq_no);
			cur->right->seq_no = NULL;
		}
        	if(cur->right != NULL)
		{
            		free(cur->right);
			cur->right = NULL;
		}
	}
}


/* free_ptr() frees all the pointers. */
void free_ptr(float *row, float *N_X, float *f_X, char *consensus, long long int *ci_row, long long int *ci_inv_row, long long int *prev, float *total_check, long long int *stack, char *line, char *k_mer)
{
    if(N_X != NULL)
    {
        free(N_X);
        N_X = NULL;
    }
    if(f_X != NULL)
    {
        free(f_X);
        f_X = NULL;
    }
    if(consensus != NULL)
    {
        free(consensus);
        consensus = NULL;
    }
    if(ci_row != NULL)
    {
        free(ci_row);
        ci_row = NULL;
    }
    if(ci_inv_row != NULL)
    {
        free(ci_inv_row);
        ci_inv_row = NULL;
    }
    if(prev != NULL)
    {
        free(prev);
        prev = NULL;
    }
    if(total_check != NULL)
    {
        free(total_check);
        total_check = NULL;
    }
    if(line != NULL)
    {
        free(line);
        line = NULL;
    }
    if(stack != NULL)
    {
        free(stack);
        stack = NULL;
    }
    if(row != NULL)
    {
        free(row);
        row = NULL;
    }
    if(k_mer != NULL)
    {
	free(k_mer);
        k_mer = NULL;
    }
}

