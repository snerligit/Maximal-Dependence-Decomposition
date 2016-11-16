int no_groups;
char *groups;

long long int limit;
long long int pseudo_n, pseudo_d;

static float df[2] = {45.5579,65.4750};

struct node
{
    struct node *left;
    struct node *right;
    char **seq;
    long long int position;
    char nucleotide;
    float **pwm;
    long long int no_ele;
    long long int *seq_no;
};
typedef struct node * NODE;

/* Definition of frees in file free.c */
void freenode(NODE set);
void free_ptr(float *row, float *N_X, float *f_X, char *consensus, long long int *ci_row, long long int *ci_inv_row, long long int *prev, float *total_check, long long int *stack, char *line, char *k_mer);
void free_children(NODE cur);

/* Definition of computation in file computation.c */
float square(float x);
float compute_chi(float N,float *N_X,float *f_X);
void find_consensus(char *con, char **sequences, long long int row, long long int column);
long long int check(float *arr, long long int row);
void create_pwm(NODE root, int no_groups);
int nucleotide(char ch);
NODE compute_weights(NODE root, int col, int no_groups);

/* Definition of tree operations in tree.c */
NODE getnode();
void display(NODE tree);
NODE split(NODE set, char nucleotide, long long int pos, long long int row, long long int col, long long int *ci_row, long long int *ci_inv_row, char max_nucltd, long long int max_pos, int is_root);
float traverse(NODE root, char *sequence);
void formation(NODE root, long long int n, float *row, float *N_X, float *f_X, char *consensus, float pseudo_n, float pseudo_d,long long int no_sequences, long long int N, long long int *ci_row, long long int *ci_inv_row, long long int r, long long int ri,long long int *prev,float *total_check, long long int *stack, long long int top, long long int dir, int df_ind);

