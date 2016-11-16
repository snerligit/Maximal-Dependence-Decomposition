#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"functions.h"

/* getnode() is a function which allocates memory for a new node. */
NODE getnode()
{
    NODE temp;
    temp=(NODE)malloc(sizeof(struct node));
    if(temp==NULL)
    {
        printf("Memory allocation failed");
        return NULL;
    }
    return temp;
}

void display(NODE tree)
{
    int i, j, k, col;
    FILE *fptr = NULL;
    if(tree == NULL)
        return;
    fptr = fopen("PWM.txt","a");
    if(tree->nucleotide != '\0' && tree->left != NULL && tree->right != NULL)
    {
        fprintf(fptr, "\n%lld|%c|%lld\nPWM\n-------\n",tree->position+1,tree->nucleotide,tree->no_ele);
	for(i=0;i<no_groups;i++)
	{
		fprintf(fptr, "%f\n",tree->pwm[i][0]);
	}
    }
    else
    {
	if(tree->nucleotide != '\0')
	{
		fprintf(fptr, "Node without one child\n");
		fprintf(fptr, "\n%lld|%c|%lld\nPWM\n-------\n",tree->position+1,tree->nucleotide,tree->no_ele);
	}
	else
		fprintf(fptr, "Leaf|%lld\n",tree->no_ele);
	col = strlen(tree->seq[0]);
	for(i=0;i<no_groups;i++)
	{
		for(j=0;j<col;j++)
			fprintf(fptr, "%f\t",tree->pwm[i][j]);
		fprintf(fptr,"\n");
	}
	/*printf("\nSequences\n");
	for(k=0;k<tree->no_ele;k++)
	{
		printf("%s %lld\n", tree->seq[k], tree->seq_no[k]);
	}*/
    }
    fclose(fptr);
    display(tree->left);
    display(tree->right);
}

/* Give only one set to split that is only one node. Based on the nucleotide at that position,
   sequences are split into two groups which are the children of set.

   Parameters:		set: It is a root for which we need to find children.
			nucleotide and pos: based on nucleotide and position, splitting happens.
			row and col: rows and columns of set(i.e, no. of sequences and no. of nucleotides in each sequence).
			ci_row and ci_inv_row: attributes which hold the count of sequences in ci and ci_inv.
*/
NODE split(NODE set, char nucleotide, long long int pos, long long int row, long long int col, long long int *ci_row, long long int *ci_inv_row, char max_nucltd, long long int max_pos, int is_root)
{
    long long int i,j,k;
    /* lchild and rchild are the children of set that are created by split. */
    NODE lchild = NULL,rchild = NULL;
    long long int lflag = 0,rflag = 0;

    if(set==NULL)return set;

    lchild = getnode();
    rchild = getnode();

    lchild->seq = (char**)malloc(sizeof(char*)*row);
    rchild->seq = (char**)malloc(sizeof(char*)*row);

    lchild->seq_no = (long long int*)malloc(sizeof(long long int)*row);
    rchild->seq_no = (long long int*)malloc(sizeof(long long int)*row);

    /* If the nucleotide matches the position for a particular sequence, add it to left child else to the right child. */
    j = 0;
    k = 0;
    for(i=0; i<row; i++)
    {
        if(set->seq[i][pos] == nucleotide)
        {
            lflag = 1;
            lchild->seq[j] = (char*)malloc(col*sizeof(char));
            strcpy(lchild->seq[j], set->seq[i]);
	    lchild->seq_no[j] = set->seq_no[i];
            j++;

        }
        else
        {
            rflag = 1;
            rchild->seq[k] = (char*)malloc(col*sizeof(char));
            strcpy(rchild->seq[k], set->seq[i]);
	    rchild->seq_no[k] = set->seq_no[i];
            k++;
        }
    }

    lchild->right = NULL;
    lchild->left = NULL;

    rchild->right = NULL;
    rchild->left = NULL;

    if(lflag == 0)
    {
        free(lchild->seq);
        lchild->seq = NULL;
	free(lchild);
	lchild = NULL;
    }
    if(rflag == 0)
    {
        free(rchild->seq);
        rchild->seq = NULL;
	free(rchild);
	rchild = NULL;
    }

    if(is_root)
    {
	if((j != 0 && j <= limit) || (k != 0 && k <= limit))
		return set;
    }

    /* Attach the children to the parent and return the parent to the calling function. */
    set->left = lchild;
    set->right = rchild;

    set->nucleotide = max_nucltd;
    set->position = max_pos;

    /* Imp: ci_row and ci_inv_row are pointers because we need the no. of sequences in each
       of ci and ci_inv and the only way to get it is through this function. It is because,
       this is the function that creates these sets. */

    if(lchild != NULL)
    {
    	set->left->no_ele = j;
	set->left->nucleotide = '\0';
	set->left->position = -1;
    }
    if(rchild != NULL)
    {
	set->right->no_ele = k;
	set->right->nucleotide = '\0';
	set->right->position = -1;
    }
    *ci_row = j;
    *ci_inv_row = k;

    return set;
}

/* Rewrite traverse so that it scores the sequence */
float traverse(NODE root, char *sequence)
{
    float score = 0.0;
    int nuc = -1, i, j;
    char *temp_seq;
    if(root == NULL)
        return 0;

    temp_seq = (char *)malloc(sizeof(char)*(strlen(sequence)+1));
    strcpy(temp_seq, sequence);

    while(root != NULL)
    {
	for(j=0;j<no_groups;j++)
		if(root->nucleotide == groups[j])
			nuc = j;

        if(root->position != -1 && root->nucleotide == sequence[root->position])
	{
	    if(root->left != NULL && root->right != NULL)
	    	score = score + root->pwm[nuc][0];
	    else
		score = score + root->pwm[nuc][root->position];
	    sequence[root->position] = 'X';
	    root = root->left;
	}
        else
	{
		if(root->right == NULL || root->nucleotide == '\0')
		{
			for(i=0;i<strlen(sequence);i++)
			{
				if(sequence[i] == 'X')
				{
					//Do nothing
				}
				else
				{
					for(j=0;j<no_groups;j++)
						if(sequence[i] == groups[j])
							nuc = j;
					score = score + root->pwm[nuc][i];
				}
			}	
		}
		root = root->right;
	}
    }
    strcpy(sequence, temp_seq);
    free(temp_seq);
    return score;
}

/* formation() function will build a tree. */
void formation(NODE root, long long int n, float *row, float *N_X, float *f_X, char *consensus, float pseudo_n, float pseudo_d,long long int no_sequences, long long int N, long long int *ci_row, long long int *ci_inv_row, long long int r, long long int ri,long long int *prev,float *total_check, long long int *stack, long long int top, long long int dir, int df_ind)
{
    long long int pos_i = 0,pos_j = 0, max_index = 0;
    long long int i,j,k,l;
    float max = -1, total;
    long long int flag = 0;
    NODE cur = NULL;

    if(root == NULL || root->seq == NULL)
        return;
    if(root->no_ele <= limit)
    {
	prev[stack[top--]] = -1;
	return;
    }
    find_consensus(consensus, root->seq, no_sequences, n);
    //printf("Consensus:%s\n",consensus);

    /* Loop through till you fill the entire chi-squared table. */
    for(pos_i=0; pos_i<n; pos_i++)
    {
	cur = root;
        /* based on consensus nucleotides, split each time for each row and fill the table. */
        cur = split(root,consensus[pos_i],pos_i,no_sequences,(n+1),ci_row,ci_inv_row,'\0',0,0);

        /* Add the pseudocount (pseudo_n) to N. */
        N = *ci_row + pseudo_d;

        /* This loop will calculate the row of the table one by one. */
        for(pos_j=0; pos_j<n; pos_j++)
        {
            if(pos_i == pos_j)
            {
                row[pos_j] = 0.0;
            }
            else
            {
                for(i=0; i<no_groups; i++)
                {
                    N_X[i] = 0;
                    f_X[i] = 0;
                }
                /* To find N_X. */
                for(i=0; i<*ci_row; i++)
                {
		    for(l=0;l<no_groups;l++)
			if(cur->left->seq[i][pos_j] == groups[l])
				N_X[l]++;
                }
		for(l=0;l<no_groups;l++)
                	N_X[l] += pseudo_n;
                
                /* To find f_X. */
                for(i=0; i<*ci_inv_row; i++)
                {
		    for(l=0;l<no_groups;l++)
			if(cur->right->seq[i][pos_j] == groups[l])
				f_X[l]++;
                }
		for(l=0;l<no_groups;l++)
                	f_X[l] = (f_X[l] + pseudo_n)/(*ci_inv_row + pseudo_d);
                
                row[pos_j] = compute_chi(N,N_X,f_X);
                //printf("i = %lld j = %lld %f\t",pos_i,pos_j,row[pos_j]);
            }
        }
        total = 0.0;
        for(i=0; i<n; i++)
        {
            total += row[i];
        }
        //printf("total:%f\n",total);
        if(max - total < 0.0)
        {
            if(prev[pos_i] != 1)
            {
                max = total;
                max_index = pos_i;
            }
        }
        //printf("\n");
        free_children(cur);
	cur = NULL;
    }
    /*printf("Conserved nucleotide: %c in position %lld\n",consensus[max_index],max_index+1);
    printf("Max:%f\n",max);*/
    root = split(root,consensus[max_index],max_index,no_sequences,(n+1),ci_row,ci_inv_row,consensus[max_index],max_index,1);
    r = *ci_row;
    ri = *ci_inv_row;

    if((r <= limit && r != 0) || (ri <= limit && ri != 0))
	return;

    if(dir == 0)
    {
    	stack[++top] = max_index;
    }
    total_check[max_index] = max;
    prev[max_index] = 1;

    if(r <= limit && ri <= limit && (max - df[df_ind]) < 0.0 )
    {
	prev[stack[top--]] = -1;
        return;
    }
    else
    {
        no_sequences = *ci_row;

        formation(root->left, n, row, N_X, f_X, consensus, pseudo_n, pseudo_d, no_sequences, N, ci_row, ci_inv_row, r, ri, prev,total_check, stack, top, 0, df_ind);
        *ci_row = r;
        *ci_inv_row = ri;
    }

    if(r <= limit && ri <= limit && (max - df[df_ind]) < 0.0 )
        return;
    else
    {
        no_sequences = *ci_inv_row;

        formation(root->right, n, row, N_X, f_X, consensus, pseudo_n, pseudo_d, no_sequences, N, ci_row, ci_inv_row, r, ri, prev,total_check, stack, top, 1, df_ind);
        *ci_row = r;
        *ci_inv_row = ri;
    }
}
