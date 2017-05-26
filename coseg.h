//
// Prototypes
//
void usage(void);
void initialize_subfamily0();
void allocate_memoryS();
int build_sequence( char *filename );
void allocate_memory();
int count_seqs( char *filename );
double inverseNormalCDF( double x );
void build_eles( char *filename );
void compute_tri_bestmut();
void compute_bestmut();
void compute_mutfrac( int SS );
void compute_localbicount( int SS );
int intcmp( const void *v1, const void *v2 );
int getTRICount( int SS, int x, int a, int xx, int aa, int xxx, int aaa );
double compute_siegel_pvalue( int n, int n1, int n2, int n12, double epsilon );
double getDoubleEpsilon();
double compute_pvalue( int totalcount, int count1, int count2, 
                       int count12, double fudge );
double compute_score( int totalcount, int count1, int count2, int count12 );
void build_new_tri_subfamily() ;
void build_new_subfamily2();
void build_MST_full( char *filename );
double sigmage_to_logpvalue( double sigmage );
void build_singlemut_MST();
void build_new_singlemut_subfamily();
void build_new_subfamily();
void run_em();
void assign_to_pattern();
void assign_to_pattern_singlemut();
void assign_to_pattern_mark();
void pattern_to_assign_mark( int thisx );
void pattern_to_assign();
void build_local( int thisx, int thisa, int thisxx, int thisaa, int SS );
void localpattern_to_localassign( int SS );
void localassign_to_localpattern( int SS );
void print_singlemut_subfamilies( int S );
void print_subfamilies( int S );
char char_to_num( char c );
char num_to_char( char z );
double split_pvaluelocal( int SS, int s, int w, double pvaluehope );
void compute_bestmut1();
void print_assign( char *filename );
void prune_subfamilies();
void compute_numerdenom( int s );
double union_tri_pvalue( int label );
double union_pvalue( int label );
void merge_subfamilies( int label );
double compute_sigma( int totalcount, int count1, double emutfrac );
int compute_distance( int s, int t, int insPenalty );
void build_MST_scaffold( char *filename );
void print_subfamily( int s );
void build_tri_local( int thisx, int thisa, int thisxx, int thisaa,
                      int thisxxx, int thisaaa, int SS );
double compute_siegel_tri_pvalue( int n, int n1, int n2, int n3,
                                  int n12, int n13, int n23,
                                  int n123, double epsilon );
double compute_tri_pvalue( int totalcount, int count1, int count2,
                           int count3, int count12, int count13,
                           int count23, int count123, double fudge );
double compute_tri_score( int totalcount, int count1, int count2, int count3,
                          int count123 );
void computeMutationRatesKimura(double *mutrate, double *mutrate_CpGMod);

int
build_local_global(int thisx, int thisa, int thisxx, int thisaa, int SS);

 
