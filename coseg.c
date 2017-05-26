/**
 ** COSEG
 **
 ** A program to identify repeat subfamilies using significant
 ** co-segregating ( 2-3 bp ) mutations.  
 **
 ** This program is derived from three C programs and several perl
 ** scripts written by Alkes Price as part of an analysis of Alu 
 ** elements in the human genome ( "Whole-genome analysis of Alu repeat
 ** elements reveals complex evolutionary history", Alkes L. Price, 
 ** Eleazar Eskin, and Pavel Pevzner, 2004 Genome Research ). The program 
 ** was first adapted for use with other repeat families and then 
 ** extended to support consideration of three co-segregating mutations
 ** using Alkes statistical model.  In 2008 with the help of Andy Siegel
 ** an alternative statistical model was developed and the codebase
 ** repackaged into the single source file.
 **
 ** Major changes from the original codebase include:
 **    - Removal of ALU-specific hardcoded parameters
 **           o Ignore consensus positions bp 116-135 ( ALU variable region )
 **             NOTE: This functionality was generalized and any region(s)
 **                   can be disabled in this analysis.
 **           o Constant sigma threshold for single point mutations is
 **             is now calculated given the desired pValue, conLen, and
 **             a conservative choice for the number of tests being performed.
 **           o The consensus length is calculated based on the input file
 **    - Command line parameters rather than hard-coded filenames.
 **    - Repbase comparison analysis removed.
 **    - Output files streamlined.
 **
 ** Current Authors: Robert Hubley
 **                  Andy Siegel
 **                  Arian Smit
 **
 ** Original Author: Alkes Price
 **
 **/
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <ctype.h>
#include "coseg.h"

extern const char *Version;

//
// Critical Parameters
//
#define MAXS 400                //  Max number of subfamilies to build
#define DISALLOW_DASHDASH 1     //  Do not allow two indels to be diagnostic
#define MINDISTDEFAULT 10       //  Default Min distance between diagnostic sites
#define SCORETHRESH  3.0        //  Minimum number of standard deviations away
                                //    from mean before mutation will get scored.
#define SCAFFOLD_PVALUE_THRESH .001     //  The pValue significance threshold for
                                        //    scaffolds ( bi/tri mutations )
#define DISALLOW_CG  1          //  Filter out CpG's in first diagnostic site
#define DISALLOW_CG2 1          //  Filter out CpG's in second diagnostic site
#define DISALLOW_CG3 1          //  Filter out CpG's in third diagnostic site
#define SINGLE_MUTATION_INSPENALTY 11   // Penalty for a single insertion
                                        //   during the single mutation 
                                        //   search phase.
#define MULTIPLE_MUTATION_INSPENALTY 1  // Penalty for a single insertion
                                        //   during the bi/tri mutation
                                        //   search phase.

#define PVALUETHRESH log(SCAFFOLD_PVALUE_THRESH)        // Convert threshold to
                                                 //     logspace

#define MAXITER 100             // Maximum iterations for scaffold loop
#define BIGITER 10
#define VERBOSE 0               // Deprecated: Use DEBUG now

// Single mutation specific parameters
#define SINGLE_PVALUE_THRESH .001       //  The pValue significance threshold for
                                        //      single point mutations
#define NEW_CG 0                // Experimental filter of CG sites when
                                //   considering single mutations. TODO:
                                //   Test this more and document.
#define DISALLOW_INDEL 1        // flag to disallow indel in single diagnostic test
#define NEW_MUT_ONLY 1          // flag to disallow mut to value already present
                                //    -- only used in adding new singlmut families
                                //       Double and triple do this by default.


// 
// GLOBAL VARIABLES
//
int useTRI = 0;
int DEBUG = 0;                  // DEBUG Levels:  0 - Quiet
                                //                1 - Overall state
                                //                2 - Detailed state
                                //                3 - Subroutine Calls
char **ele;                     // maxN*L  0=A 1=C 2=G 3=T 4=- 
char **elei;                    // maxN*L  0=  1=+ 
int N;
int maxN;                       // The total number of sequences in 
                                //    the input file. 
int assigncount[MAXS];          // Number of elements in a subfamily
int trimut[MAXS];
double mstLogPValues[MAXS];     // pvalue of each subfamily relative 
                                //    to [one] of it's closest neighbors
                                //    ( consensus distance ).
int *assign;                    // Element to subfamily assignment array
int *sortedAssign;
int sortedAssignIndices[MAXS];
int *localassign;
double *logfactorial;
int localassigncount[2];
char **pattern;                 // MAXS*L : The most frequent char at pos x in m.a.
char **patterni;                // MAXS*L 
char **localpattern;            // 2*L 
char **localpatterni;           // 2*L 
int parent[MAXS];               // Subfamily to parent assignment array
int bestmutSS, bestmutx, bestmuta, bestmutxx, bestmutaa;
int bestmutxxx, bestmutaaa;
double bestmutpvalue;
char *Sxsequence;               // Consensus sequence passed to program
int S, MALLOCS;
int totdist, oldtotdist, localtotdist, localoldtotdist;
int ***count;                   /* MAXS*L*5, count of position x value a */
int ***counti;                  /* MAXS*L*2, count of position x value b */
int *****bicount;               /* MAXS*L*5*L*5, count of pos x val a pos xx val aa */
int ***localcount;              /* 2*L*5, count of position x value a */
int ***localcounti;             /* 2*L*2, count of position x value b */
int *****localbicount;          /* 2*L*5*L*5, count of pos x val a pos xx val aa */
double sigmaThreshold;          // Number of std devs away for single point mutation
                                // to be significant ( pVal 0.001 )
double sngBonferroni;           // Single mutation bonferroni
double dblBonferroni;           // Double mutation bonferroni
double triBonferroni;           // Tripple mutation bonferroni
int ***numer, ***denom;         /* MAXS*L*5 */
double **mutfrac;               /* L*5 */
int labels[MAXS];
int **distance;                 /* MAXS*MAXS */
int we_want_best_pvalue;
int **existingval;              /* L*5 */
char *seqFile;                  /* Input sequences file */
char *conFile;                  /* Consensus sequence file */
char *logFile;
int minCount;                   /* min count of a new subfamily */
int minDist;                    /* min distance between diagnostic sites. min = 1 */
int conLen;
FILE *outFP;
double dblEpsilon;              // Minimum number that a double can hold
int useSiegelPvalue = 1;        // Flag to change statistical models
int useOriginalGraphColors = 0; // Turn on the original graph heatmap colors.
int useDisabledSites = 0;       // Turn on special meaning of lower/upper
                                //   case characters in the consensus file.
int disabledPosCount = 0;       // Number of positions in cosensus which are
                                //   disabled ( see below ).
char *disabledSites;            // An array to hold a flag per consensus position
                                //   which indicates that the site has been 
                                //   disabled by the user. I.e there may be a
                                //   highly variable region inside an consensus
                                //   which shouldn't be used in this type of
                                //   analysis.
int numScaffolds = 0;           // The number of total subfamilies that were
                                //   built using co-segregating mutations.

// For singlemut
int lastScaffoldIndex = 0;
int mark[MAXS];
double bestmutsigma;
double mutsigma[MAXS];          // Holds bestmutsigma for all 1-mut families
double ***profile;              /* S*L*5 0=A 1=C 2=G 3=T 4=- */
double ***profilei;             /* S*L*2 0=  1=+ */
double age[MAXS];
double ***mutperage;            /* 5*L*5 */
double ***mutperagei;           /* 2*L*2 */
double ***mut;                  /* 5*L*5 */
double ***muti;                 /* 2*L*2 */
double **perage;                /* 5*L */
double **peragei;               /* 2*L */
int **edges;                    /* MAXS*MAXS* */
char *graphVizFile;
char *outAssignFile;
char *subfamFile;
// end for singlemut

// 
// EXPERIMENTAL : New EM for intermediate check of sites
//
  int *** l_count;
  int *** l_counti;
  int * l_assign;
  char ** l_pattern;
  char ** l_patterni;
//
// End Experimental
//



static char *options[] = {
  "\nusage: %s [-m #] [-t] [-k] [-d] [-u #] -c <consensus file>\n",
  "          -s <sequence file> \n",
  " Options\n",
  " -c <file> - Fasta file containing the consensus sequence.\n",
  " -s <file> - Sequence file containing one sequence per line\n",
  " -t        - Use tri-segregating mutations in calculation\n",
  " -m #      - The minimum number of elements in a subfamily ( default 50 )\n",
  " -k        - Use original Alkes Price statistical model to calculate pValues\n",
  "             Default: Use Andy Siegel's pvalue calculation instead.\n",
  " -d        - Use character case in consensus to flag disabled regions\n",
  " -o        - Use original colormap for graphviz/svg graphs: warm=young,\n",
  "             cold=older more highly diverged subfamiles.  The new default\n",
  "             reverses this scheme: cold=young, warm=older.\n",
  " -u #      - Minimum distance between diagnostic sites ( min 1, default 10 )\n",
  " -v #      - Verbosity. Use values 0-3 for increasing verbosity ( default 0 )\n\n",
  "   A program to identify repeat subfamilies using significant\n",
  " co-segregating ( 2-3 bp ) mutations.\n\n",
  "   This program is derived from three C programs and several perl\n",
  " scripts written by Alkes Price as part of an analysis of Alu \n",
  " elements in the human genome ( Whole-genome analysis of Alu repeat\n",
  " elements reveals complex evolutionary history, Alkes L. Price,\n",
  " Eleazar Eskin, and Pavel Pevzner, 2004 Genome Research ). The program\n",
  " was first adapted for use with other repeat families and then\n",
  " extended to support consideration of three co-segregating mutations\n",
  " using Alkes statistical model.  In 2008 with the help of Andy Siegel\n",
  " an alternative statistical model was developed and the codebase\n",
  " repackaged into the single source file.\n\n",
  " Current Authors: Robert Hubley\n",
  "                  Andy Siegel\n",
  "                  Arian Smit\n\n",
  " Original Code Author: Alkes Price\n",
  0,
};


//
// MAIN
//
int
main(int argc, char **argv)
{
  int origMinCount;
  time_t start, finish;
  int n, s, t, x, a, xx, aa, S_BEFORE_PRUNE, S_AFTER_PRUNE[BIGITER];
  int iter, length, B;
  double duration;
  double minCountFactor;


  // 
  // Option Processing
  //
  int option = 0;
  minCount = 50;
  minDist = MINDISTDEFAULT;
  while ((option = getopt(argc, argv, "c:s:m:u:v:tkd")) != EOF)
  {
    switch (option)
    {
    case 'c':
      conFile = (char *) malloc((strlen(optarg) + 1) * sizeof(char));
      strcpy(conFile, optarg);
      break;
    case 's':
      seqFile = (char *) malloc((strlen(optarg) + 1) * sizeof(char));
      strcpy(seqFile, optarg);
      // Log file -- TODO: Consider working this into other files
      logFile = (char *) malloc((strlen(optarg) + 10) * sizeof(char));
      strcpy(logFile, optarg);
      strcat(logFile, ".log");
      // Main output file
      subfamFile = (char *) malloc((strlen(optarg) + 13) * sizeof(char));
      strcpy(subfamFile, optarg);
      strcat(subfamFile, ".subfamilies");
      // Main graph file
      graphVizFile = (char *) malloc((strlen(optarg) + 10) * sizeof(char));
      strcpy(graphVizFile, optarg);
      strcat(graphVizFile, ".tree.viz");
      // Main assignment file
      outAssignFile = (char *) malloc((strlen(optarg) + 20) * sizeof(char));
      strcpy(outAssignFile, optarg);
      strcat(outAssignFile, ".assign");
      break;
    case 'o':
      useOriginalGraphColors = 1;
      break;
    case 'm':
      sscanf(optarg, "%d", &minCount);
      if (minCount <= 2)
      {
        printf("-m ( minCount ) must be >= 3\n");
        exit(1);
      }
      break;
    case 'u':
      sscanf(optarg, "%d", &minDist);
      if (minDist < 1)
      {
        printf("-u ( minDist ) must be >= 1\n");
        exit(1);
      }
      break;
    case 'v':
      sscanf(optarg, "%d", &DEBUG);
      break;
    case 'k':
      useSiegelPvalue = 0;
      break;
    case 'd':
      useDisabledSites = 1;
      break;
    case 't':
      useTRI = 1;
      break;
    default:
      usage();
      break;
    }
  }
  if (!conFile || !seqFile)
  {
    usage();
  }

  //
  if ((outFP = fopen(logFile, "w")) == NULL)
  {
    printf("Cannot open %s for writing!\n", logFile);
    exit(1);
  }

  start = time(0);
  length = build_sequence(conFile);
  if (disabledPosCount == length)
  {
    printf
      ("ERROR: The consensus supplied to the program has all lowercase\n"
       "       (disabled) bases.  No positions are available for\n"
       "       analysis.  Please edit the consensus file and change\n"
       "       some bases to upper case and resubmit.\n");
    exit(1);
  }
  conLen = length;

  maxN = count_seqs(seqFile);
  printf("COSEG - %s\n", Version);
  printf("Number of sequences = %d\n", maxN);
  if (useDisabledSites)
    printf("Consensus length = %d  ( %d disabled positions )\n",
           conLen, disabledPosCount);
  else
    printf("Consensus length = %d\n", conLen);

  printf("Min dist between sites = %d\n", minDist);
  printf("Min subfamily size = %d\n", minCount);

  allocate_memory();
  build_eles(seqFile);
  initialize_subfamily0();

  origMinCount = minCount;
  minCountFactor = 0.20;

  //
  // Calculate machine epsilon for doubles
  //
  dblEpsilon = getDoubleEpsilon();
  printf("Machine epsilon (double) = %le\n\n\n", dblEpsilon);

  //
  //  Main loops
  //
  for (B = 0; B < BIGITER; B++)
  {
    fprintf(outFP, "ADDING SUBFAMILIES\n");

    // TRI Mutations
    if (useTRI)
    {
      if ((((double) assigncount[S - 1]) * minCountFactor) > origMinCount)
      {
        minCount = (int) (((double) assigncount[S - 1]) * minCountFactor);
      }
      while (S < MAXS)
      {
        if (S == MALLOCS)
          allocate_memoryS();
        if (S > MALLOCS)
        {
          printf("OOPS S=%d MALLOCS=%d\n", S, MALLOCS);
          exit(1);
        }

        printf("\nComputing tri-mutations S=%d minCount=%d out of %d:\n",
               S, minCount, assigncount[S - 1]);
        compute_tri_bestmut();
        if (bestmutpvalue >= PVALUETHRESH)
        {
          printf("  No significant mutations at this minCount.\n");
          if (minCount == origMinCount)
            break;
          minCount *= 0.60;
          if (minCount < origMinCount)
            minCount = origMinCount;
          printf(" Trying new minCount = %d\n", minCount);
        }
        else
        {
          printf
            ("  Found significant mutation pvalue = %lf ( < %lf )\n",
             bestmutpvalue, PVALUETHRESH);

          build_new_tri_subfamily();    /* parent, S, assign, pattern_to_assign */
          if ((((double) assigncount[S - 1]) * minCountFactor) > origMinCount)
          {
            minCount = (int) (((double) assigncount[S - 1]) * minCountFactor);
            printf(" Trying new minCount = %d out of N=%d\n",
                   minCount, assigncount[S - 1]);
          }
        }

      }
      minCount = origMinCount;
    }

    // BI Mutations
    while (S < MAXS)
    {
      if (S == MALLOCS)
        allocate_memoryS();
      if (S > MALLOCS)
      {
        printf("OOPS S=%d MALLOCS=%d\n", S, MALLOCS);
        exit(1);
      }

      printf("\nComputing bi-mutations:\n");
      compute_bestmut();
      if (bestmutpvalue >= PVALUETHRESH)
      {
        printf("  No significant mutations\n");
        break;
      }
      else
      {
        printf("  Found significant mutation pvalue = %lf ( < %lf )\n",
               bestmutpvalue, PVALUETHRESH);
        build_new_subfamily();  /* parent, S, assign, pattern_to_assign */
      }
      fflush(outFP);
      sync();
    }

    fprintf(outFP, "VERIFYING SUBFAMILIES (PRUNE IF NECESSARY)\n");
    printf("\nVERIFYING %d SUBFAMILIES (PRUNE IF NECESSARY)\n", S);

    S_BEFORE_PRUNE = S;
    prune_subfamilies();
    S_AFTER_PRUNE[B] = S;
    if (S == S_BEFORE_PRUNE)
    {
      fprintf(outFP, "S=%d BEFORE/AFTER PRUNE\n", S);
      break;
    }
    else
    {
      printf("%d Subfamilies pruned %d remaining\n", (S_BEFORE_PRUNE - S), S);
    }
    if ((B > 0) && (S == S_AFTER_PRUNE[B - 1]))
    {
      fprintf(outFP, "S_AFTER_PRUNE=%d AFTER BIG ITER\n", B, S);
      break;
    }
  }

  fprintf(outFP, "THERE ARE %d SUBFAMILIES IN SCAFFOLD\n", S);
  printf("\nTHERE ARE %d SUBFAMILIES IN SCAFFOLD\n", S);

  we_want_best_pvalue = 1;

  // Build tree for scaffold subfamilies only
  //   - This calculates the p_value for scafffold members distinctly
  //     from the one contianing the scaffold and single point mutations.
  //     Saves results in mstLogPValues() for use by build_MST_full().
  build_MST_scaffold(NULL);

  numScaffolds = S;

  duration = difftime(time(0), start);
  printf("Scaffold duration is %.1f sec = %.1f min = %.1f hr\n",
         duration, duration / 60.0, duration / 3600.0);

  // Singlemut stuff
  printf("\nRunning single mutation algorithm...\n");
  assign_to_pattern_singlemut(); // M-step first?
  lastScaffoldIndex = S - 1;

  for (s = 0; s < S; s++)
  {
    for (t = 0; t < S; t++)
      edges[s][t] = 0;
  }

  // 
  // add internal 1-mutations (along existing edge of tree)  
  // 
  build_singlemut_MST();


  /*
     add new 1-mutations 
   */
  while (S < MAXS)
  {
    compute_bestmut1();         /* bestmutx, bestmuta, bestmutSS */
    if (bestmutsigma <= sigmaThreshold)
      break;
    mutsigma[S] = bestmutsigma;
    build_new_subfamily2();
  }
  build_MST_full(graphVizFile);
  print_subfamilies(S);

  print_assign(outAssignFile);

  // END Singlemut stuff

  finish = time(0);
  duration = difftime(finish, start);
  fprintf(outFP, "Program duration is %.1f sec = %.1f min = %.1f hr\n",
          duration, duration / 60.0, duration / 3600.0);
  printf("Program duration is %.1f sec = %.1f min = %.1f hr\n\n",
         duration, duration / 60.0, duration / 3600.0);
}




/************************ S U B R O U T I N E S ************************/


//
// Print program usage
// 
void
usage(void)
{
  int i;
  char *s;
  for (i = 0; (s = options[i]); ++i)
  {
    if (i == 0)
    {
      printf(s, "coseg");
    }
    else
    {
      printf("%s", s);
    }
  }
  printf("\nVersion: %s\n\n", Version);
  exit(1);
}


//
//
//
void
initialize_subfamily0()
{
  int x, a, xx, xxx, aa, aaa, n;

  if (DEBUG == 3)
    printf("initialize_subfamily0(): Called\n");

  parent[0] = -1;
  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      count[0][x][a] = 0;
    // RMH added 11/8/12
    for (a = 0; a < 2; a++)
      counti[0][x][a] = 0;
  }

  // Clear bicount array
  for (x = minDist; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
    {
      for (xx = 0; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
          bicount[0][x][a][xx][aa] = 0;
      }
    }
  }
  // Initialize count/bicount/assign arrays
  for (n = 0; n < N; n++)
  {
    for (x = 0; x < conLen; x++)
    {
      count[0][x][ele[n][x]]++;
      // RMH added 11/8/12
      counti[0][x][elei[n][x]]++;
    }
    for (x = minDist; x < conLen; x++)
    {
      for (xx = 0; xx <= x - minDist; xx++)
        bicount[0][x][ele[n][x]][xx][ele[n][xx]]++;
    }
    assign[n] = 0;
  }

  assigncount[0] = N;
  S = 1;
  assign_to_pattern();
}                               // initialize_subfamily0(...


//
// Allocate memory assoicated with adding a new subfamily
//
void
allocate_memoryS()
{
  int x, a, xx, xxx, aa, aaa;

  if (DEBUG == 3)
    printf("allocate_memoryS(): Called\n");

  /*
     allocate bicount[S] 
   */
  if ((bicount[S] = (int ****) malloc(conLen * sizeof(*bicount[S]))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (x = minDist; x < conLen; x++)
  {
    if ((bicount[S][x] =
         (int ***) malloc(5 * sizeof(*bicount[S][x]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (a = 0; a < 5; a++)
    {
      if ((bicount[S][x][a] =
           (int **) malloc((x - minDist + 1) *
                           sizeof(*bicount[S][x][a]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
      for (xx = 0; xx <= x - minDist; xx++)
      {
        if ((bicount[S][x][a][xx] =
             (int *) malloc(5 * sizeof(*bicount[S][x][a][xx]))) == NULL)
        {
          printf("Out of memory\n");
          exit(1);
        }
      }
    }
  }

  MALLOCS = S + 1;
}                               // allocate_memoryS(...


//
//
//
void
allocate_memory()
{
  int x, s, a, xx, aa, xxx, aaa, n, w, cc, dd;
  double sum;
  char c;
  FILE *fp;

  if (DEBUG == 3)
    printf("allocate_memory(): Called\n");

  // Based on N
  if ((assign = (int *) malloc(maxN * sizeof(int))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((localassign = (int *) malloc(maxN * sizeof(int))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((logfactorial = (double *) malloc((maxN + 1) * sizeof(double))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }

  // compute logfactorial
  sum = 0.0;
  logfactorial[0] = 0.0;
  for (n = 1; n <= maxN; n++)
  {
    sum += log(((double) n));
    logfactorial[n] = sum;
  }
  if (DEBUG >= 1)
    printf("logFactorial computed up to %d\n", maxN);

  // Compute bonferroni & sigmaThreshold
  //   TODO: Consider better way to adjust when disabled positions
  //         are in play and for ignored CG sites
  dblBonferroni = 0.5 * ((double) (conLen - disabledPosCount - minDist)) *
    ((double) (conLen - disabledPosCount - minDist + 1)) * 4.0 * 4.0;
  dblBonferroni = log(dblBonferroni);

  triBonferroni = (((double) (conLen - disabledPosCount - minDist))
                   *
                   ((double) (conLen - disabledPosCount - minDist - minDist))
                   *
                   ((double)
                    (conLen - disabledPosCount - minDist - minDist +
                     1)) * 4.0 * 4.0 * 4.0) / ((double) 6);
  triBonferroni = log(triBonferroni);
  sngBonferroni = (conLen - disabledPosCount) * 4 * MAXS;
  sigmaThreshold = -inverseNormalCDF(SINGLE_PVALUE_THRESH / sngBonferroni);

  printf("Bonferroni = %lf, triBonferroni = %lf, sigmaThreshold = %lf\n",
         dblBonferroni, triBonferroni, sigmaThreshold);

  // Allocate the sortedAssign list for optimizing the getTRICount
  // function.
  if (useTRI)
  {
    if ((sortedAssign = (int *) malloc(maxN * sizeof(int))) == NULL)
    {
      printf("Could not allocate sortedAssign: Out of memory?\n");
      exit(1);
    }
  }

  if ((ele = (char **) malloc(maxN * sizeof(*ele))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((elei = (char **) malloc(maxN * sizeof(*elei))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }

  for (N = 0; N < maxN; N++)
  {
    if ((ele[N] = (char *) malloc(conLen * sizeof(*ele[N]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    if ((elei[N] = (char *) malloc(conLen * sizeof(*elei[N]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((pattern = (char **) malloc(MAXS * sizeof(*pattern))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((patterni = (char **) malloc(MAXS * sizeof(*patterni))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((pattern[s] = (char *) malloc(conLen * sizeof(*pattern[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    if ((patterni[s] =
         (char *) malloc(conLen * sizeof(*patterni[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((localpattern = (char **) malloc(2 * sizeof(*localpattern))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((localpatterni = (char **) malloc(2 * sizeof(*localpatterni))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < 2; s++)
  {
    if ((localpattern[s] =
         (char *) malloc(conLen * sizeof(*localpattern[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    if ((localpatterni[s] =
         (char *) malloc(conLen * sizeof(*localpatterni[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((count = (int ***) malloc(MAXS * sizeof(*count))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((count[s] = (int **) malloc(conLen * sizeof(*count[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((count[s][x] = (int *) malloc(5 * sizeof(*count[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((counti = (int ***) malloc(MAXS * sizeof(*counti))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((counti[s] = (int **) malloc(conLen * sizeof(*counti[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((counti[s][x] = (int *) malloc(2 * sizeof(*counti[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((bicount = (int *****) malloc(MAXS * sizeof(*bicount))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  S = 0;
  if ((bicount[S] = (int ****) malloc(conLen * sizeof(*bicount[S]))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (x = minDist; x < conLen; x++)
  {
    if ((bicount[S][x] =
         (int ***) malloc(5 * sizeof(*bicount[S][x]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (a = 0; a < 5; a++)
    {
      if ((bicount[S][x][a] =
           (int **) malloc((x - minDist + 1) *
                           sizeof(*bicount[S][x][a]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
      for (xx = 0; xx <= x - minDist; xx++)
      {
        if ((bicount[S][x][a][xx] =
             (int *) malloc(5 * sizeof(*bicount[S][x][a][xx]))) == NULL)
        {
          printf("Out of memory\n");
          exit(1);
        }
      }
    }
  }

  MALLOCS = 1;

  if ((localcount = (int ***) malloc(2 * sizeof(*localcount))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (w = 0; w < 2; w++)
  {
    if ((localcount[w] =
         (int **) malloc(conLen * sizeof(*localcount[w]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((localcount[w][x] =
           (int *) malloc(5 * sizeof(*localcount[w][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((localcounti = (int ***) malloc(2 * sizeof(*localcounti))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (w = 0; w < 2; w++)
  {
    if ((localcounti[w] =
         (int **) malloc(conLen * sizeof(*localcounti[w]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((localcounti[w][x] =
           (int *) malloc(2 * sizeof(*localcounti[w][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((localbicount =
       (int *****) malloc(2 * conLen * sizeof(*localbicount))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (w = 0; w < 2; w++)
  {
    if ((localbicount[w] =
         (int ****) malloc(conLen * sizeof(*localbicount[w]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = minDist; x < conLen; x++)
    {
      if ((localbicount[w][x] =
           (int ***) malloc(5 * sizeof(*localbicount[w][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
      for (a = 0; a < 5; a++)
      {
        if ((localbicount[w][x][a] =
             (int **) malloc((x - minDist + 1) *
                             sizeof(*localbicount[w][x][a]))) == NULL)
        {
          printf("Out of memory\n");
          exit(1);
        }
        for (xx = 0; xx <= x - minDist; xx++)
        {
          if ((localbicount[w][x][a][xx] =
               (int *) malloc(5 *
                              sizeof(*localbicount[w][x][a][xx]))) == NULL)
          {
            printf("Out of memory\n");
            exit(1);
          }
        }
      }
    }
  }

  if ((numer = (int ***) malloc(MAXS * sizeof(*numer))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((numer[s] = (int **) malloc(conLen * sizeof(*numer[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((numer[s][x] = (int *) malloc(5 * sizeof(*numer[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((denom = (int ***) malloc(MAXS * sizeof(*denom))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((denom[s] = (int **) malloc(conLen * sizeof(*denom[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((denom[s][x] = (int *) malloc(5 * sizeof(*denom[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((mutfrac = (double **) malloc(conLen * sizeof(*mutfrac))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (x = 0; x < conLen; x++)
  {
    if ((mutfrac[x] = (double *) malloc(5 * sizeof(*mutfrac[x]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  // For singlemut begin
  if ((edges = (int **) malloc(MAXS * sizeof(*edges))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((edges[s] = (int *) malloc(MAXS * sizeof(*edges[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((profile = (double ***) malloc(MAXS * sizeof(*profile))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((profile[s] =
         (double **) malloc(conLen * sizeof(*profile[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((profile[s][x] =
           (double *) malloc(5 * sizeof(*profile[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((profilei = (double ***) malloc(MAXS * sizeof(*profilei))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((profilei[s] =
         (double **) malloc(conLen * sizeof(*profilei[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((profilei[s][x] =
           (double *) malloc(2 * sizeof(*profilei[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((mutperage = (double ***) malloc(5 * sizeof(*mutperage))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (cc = 0; cc < 5; cc++)
  {
    if ((mutperage[cc] =
         (double **) malloc(conLen * sizeof(*mutperage[cc]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((mutperage[cc][x] =
           (double *) malloc(5 * sizeof(*mutperage[cc][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((mutperagei = (double ***) malloc(2 * sizeof(*mutperagei))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (dd = 0; dd < 2; dd++)
  {
    if ((mutperagei[dd] =
         (double **) malloc(conLen * sizeof(*mutperagei[dd]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((mutperagei[dd][x] =
           (double *) malloc(2 * sizeof(*mutperagei[dd][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }


  if ((mut = (double ***) malloc(5 * sizeof(*mut))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (cc = 0; cc < 5; cc++)
  {
    if ((mut[cc] = (double **) malloc(conLen * sizeof(*mut[cc]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((mut[cc][x] = (double *) malloc(5 * sizeof(*mut[cc][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((muti = (double ***) malloc(2 * sizeof(*muti))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (dd = 0; dd < 2; dd++)
  {
    if ((muti[dd] = (double **) malloc(conLen * sizeof(*muti[dd]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((muti[dd][x] = (double *) malloc(2 * sizeof(*muti[dd][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((perage = (double **) malloc(5 * sizeof(*perage))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (cc = 0; cc < 5; cc++)
  {
    if ((perage[cc] =
         (double *) malloc(conLen * sizeof(*perage[cc]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((peragei = (double **) malloc(2 * sizeof(*peragei))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (dd = 0; dd < 2; dd++)
  {
    if ((peragei[dd] =
         (double *) malloc(conLen * sizeof(*peragei[dd]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((existingval = (int **) malloc(conLen * sizeof(*existingval))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (x = 0; x < conLen; x++)
  {
    if ((existingval[x] =
         (int *) malloc(5 * sizeof(*existingval[x]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  // For singlemut end

  if ((distance = (int **) malloc(MAXS * sizeof(*distance))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((distance[s] = (int *) malloc(MAXS * sizeof(*distance[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

  if ((existingval = (int **) malloc(conLen * sizeof(*existingval))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (x = 0; x < conLen; x++)
  {
    if ((existingval[x] =
         (int *) malloc(5 * sizeof(*existingval[x]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }

// 
// EXPERIMENTAL : New EM for intermediate check of sites
//
  if ((l_count = (int ***) malloc(MAXS * sizeof(*l_count))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((l_count[s] = (int **) malloc(conLen * sizeof(*l_count[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((l_count[s][x] = (int *) malloc(5 * sizeof(*l_count[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }

  if ((l_counti = (int ***) malloc(MAXS * sizeof(*l_counti))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((l_counti[s] = (int **) malloc(conLen * sizeof(*l_counti[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      if ((l_counti[s][x] = (int *) malloc(2 * sizeof(*l_counti[s][x]))) == NULL)
      {
        printf("Out of memory\n");
        exit(1);
      }
    }
  }
  if ((l_assign = (int *) malloc(N * sizeof(int))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((l_pattern = (char **) malloc(MAXS * sizeof(*l_pattern))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  if ((l_patterni = (char **) malloc(MAXS * sizeof(*l_patterni))) == NULL)
  {
    printf("Out of memory\n");
    exit(1);
  }
  for (s = 0; s < MAXS; s++)
  {
    if ((l_pattern[s] = (char *) malloc(conLen * sizeof(*l_pattern[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
    if ((l_patterni[s] = (char *) malloc(conLen * sizeof(*l_patterni[s]))) == NULL)
    {
      printf("Out of memory\n");
      exit(1);
    }
  }
 
 
// 
// Done EXPERIMENTAL
//


}                               // void allocate_memory(...


// Ported from Andy Siegel's VBA Function
double
inverseNormalCDF(double x)
{
  // Given x < 0.5, finds the value of z such that a standard normal Z
  // satisfies P(Z<z) = x.
  // Adapted from:
  //   http://home.online.no/~pjacklam/notes/invnorm/impl/misra/normsinv.html
  // with claimed relative error rate less than 1.15*10-9
  double t, numerator, denominator;
  if (!((0 < x) && (x <= 0.5)))
  {
    printf("Warning: x = %lf, x must be between 0 and 0.5\n");
    return 0;
  }
  t = sqrt(log(1 / (x * x)));
  numerator = -7.78489400243029E-03 * t * t * t * t * t -
    0.322396458041136 * t * t * t * t -
    2.40075827716184 * t * t * t -
    2.54973253934373 * t * t + 4.37466414146497 * t + 2.93816398269878;
  denominator = 7.78469570904146E-03 * t * t * t * t +
    0.32246712907004 * t * t * t +
    2.445134137143 * t * t + 3.75440866190742 * t + 1;
  return (numerator / denominator);
}


//
// Read in the consensus sequence and initialize the disabledSites array.
//
int
build_sequence(char *filename)
{
  int i, j;
  char c;
  FILE *fp;
  int sequenceSize = 0;

  if ((fp = fopen(filename, "r")) == NULL)
  {
    printf("Could not open input file %s\n", filename);
    exit(1);
  }
  i = 0;
  while (!feof(fp))
  {                             /* process one line of file */
    c = getc(fp);
    if (c == EOF)
      continue;
    if (c == '\n')
      continue;
    if (c == '>')
    {
      for (j = 0; (getc(fp) != '\n') && !feof(fp); j++)
        ;
    }
    else
    {
      if (i == sequenceSize)
      {
        sequenceSize += 50;
        Sxsequence = realloc(Sxsequence, (sequenceSize + 1) * sizeof(char));
      }
      //if(c > 64) Sxsequence[i++] = char_to_num(c);
      if (c > 64)
        Sxsequence[i++] = c;
      for (j = 0; ((c = getc(fp)) != '\n') && !feof(fp); j++)
      {
        if (i == sequenceSize)
        {
          sequenceSize += 50;
          Sxsequence = realloc(Sxsequence, sequenceSize * sizeof(char));
        }
        //if(c > 64) Sxsequence[i++] = char_to_num(c); 
        if (c > 64)
          Sxsequence[i++] = c;
      }
    }
  }
  fclose(fp);

  // Allocate disabledSites array
  if ((disabledSites = (char *) malloc(i * sizeof(char))) == NULL)
  {
    printf("Out of memory -- allocating disabledSites.\n");
    exit(1);
  }

  for (j = 0; j < i; j++)
  {
    if (useDisabledSites && islower(Sxsequence[j]))
    {
      disabledSites[j] = 1;
      disabledPosCount++;
    }
    else
    {
      disabledSites[j] = 0;
    }
    Sxsequence[j] = char_to_num(Sxsequence[j]);
  }

  return i;                     // length of consensus 
}                               // build_sequence(...


//
// Count the number of lines in the .seqs file.
//
int
count_seqs(char *filename)
{
  FILE *fp;
  char c;
  int numSeqs = 0;
  int hasData = 0;

  if (DEBUG >= 3)
    printf("count_seqs(): Called...\n");

  if ((fp = fopen(filename, "r")) == NULL)
  {
    printf("Could not open input file %s\n", filename);
    exit(1);
  }
  while (!feof(fp))
  {
    c = getc(fp);
    if (c == EOF)
      continue;
    if (c == '\n')
    {
      if (hasData)
        numSeqs++;
      hasData = 0;
    }
    else if (c != ' ')
    {
      hasData = 1;
    }
  }
  fclose(fp);
  return (numSeqs);
}                               // int count_seqs(...


//
// Read in the .seqs file.
//
void
build_eles(char *filename)
{
  int x, s, a, xx, aa, n, w;
  double sum;
  char c;
  FILE *fp;

  if (DEBUG >= 3)
    printf("build_eles(): Called...\n");

  N = 0;
  x = 0;

  if ((fp = fopen(filename, "r")) == NULL)
  {
    printf("Could not open input file %s\n", filename);
    exit(1);
  }
  while (!feof(fp))
  {                             /* process one line of file */
    c = getc(fp);
    if (c == EOF)
      continue;
    if (c == '\n')
    {
      N++;
      if (x != conLen)
      {
        printf
          ("OOPS sequence # %d from %s is %d bp long while the consensus sequence used is %d bp long.  Please rectify this mismatch and rerun coseg.\n",
           N, filename, x, conLen);
        exit(1);
      }
      x = 0;
      continue;
    }
    if (c == '+')
    {
      elei[N][x - 1] = 1;
    }
    else
    {
      ele[N][x] = char_to_num(c);
      elei[N][x] = 0;
      x++;
    }
  }
  fclose(fp);

  return;
}                               // void build_eles(...


//
// Search all scaffolds for most significant triple mutation
//
void
compute_tri_bestmut()
{
  int n, x, a, xx, aa, totalcount, SS, w, consensus, consensuscount, t;
  int xxx, aaa;
  double pvalue23, pvalue13, pvalue123;
  int bestxx, bestaa;
  int count1, count2, count3, count12, count13, count23, count123, s;
  double pvalue, pvaluelocal, fudge, fudge2, fudge3;
  int indices[MAXS];


  if (DEBUG >= 3)
    printf("compute_tri_bestmut(): Called\n");

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      existingval[x][a] = 0;
  }

  // Setting existing val for all subfamilies
  //printf("Setting existingval\n");
  for (s = 0; s < S; s++)
    for (x = 0; x < conLen; x++)
      existingval[x][pattern[s][x]] = 1;

  // debug
  if ( 0 ) 
  {
    for (s = 0; s < S; s++)
    {
      printf("Subfamily %d: ", s );
      for (x = 0; x < conLen; x++)
      {
        if ( pattern[s][x] == 0 )
          printf("A");
        else if ( pattern[s][x] == 1 )
          printf("C");
        else if ( pattern[s][x] == 2 )
          printf("G");
        else if ( pattern[s][x] == 3 )
          printf("T");
        else 
          printf("-");
      }
      printf("\n");
    }
  }

  // Create a sorted ( by subfamily ) list of elements.
  for (x = 0; x < N; x++)
  {
    sortedAssign[x] = x;
  }
  qsort(sortedAssign, N, sizeof(sortedAssign[0]), intcmp);
  // Create a set of pointers to the start of each subfamily 
  // block within the sorted list.
  for (x = N - 1; x >= 0; x--)
  {
    sortedAssignIndices[assign[sortedAssign[x]]] = x;
  }

  bestmutpvalue = PVALUETHRESH;
  long pvalueCount = 0;
  long triCountCalcs = 0;
  printf("  Searching for significant triple mutation:");
  for (SS = 0; SS < S; SS++)    /* try splitting subfamily SS */
  {
    if (DEBUG >= 2)
      printf("Trying subfamily %d..\n", SS);
    else
      printf(".");

    // Only needed for A. Price's method
    if (!useSiegelPvalue)
      compute_mutfrac(SS);

    for (x = (minDist * 2); x < conLen; x++)
    {
      // Do not consider this site for mutation if disabled.
      if (useDisabledSites && disabledSites[x])
        continue;

      for (a = 0; a < 5; a++)   /* try 1-mutation pos x val a */
      {
        if (existingval[x][a] == 1)
          continue;
        if (count[SS][x][a] < minCount)
          continue;
        /*
           new code to disallow CG->TG etc 
         */
        if (DISALLOW_CG)
        {
          if ((pattern[SS][x] == 1) && (a == 3) && (x < conLen - 1)
              && (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 3) && (a == 1) && (x < conLen - 1)
              && (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 2) && (a == 0) && (x > 0)
              && (pattern[SS][x - 1] == 1))
            continue;
          if ((pattern[SS][x] == 0) && (a == 2) && (x > 0)
              && (pattern[SS][x - 1] == 1))
            continue;
        }

        for (xx = minDist; xx <= x - minDist; xx++)
        {

          // Do not consider this site for mutation if disabled.
          if (useDisabledSites && disabledSites[xx])
            continue;

          for (aa = 0; aa < 5; aa++)
          {
            if (existingval[xx][aa] == 1)
              continue;
            if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4))
              continue;
            /*
               new code to disallow CG->TG etc 
             */
            if (DISALLOW_CG2)
            {
              // CG -> TG
              if ((pattern[SS][xx] == 1) && (aa == 3)
                  && (xx < conLen - 1) && (pattern[SS][xx + 1] == 2))
                continue;
              // TG -> CG
              if ((pattern[SS][xx] == 3) && (aa == 1)
                  && (xx < conLen - 1) && (pattern[SS][xx + 1] == 2))
                continue;
              // CG -> CA 
              if ((pattern[SS][xx] == 2) && (aa == 0) && (xx > 0)
                  && (pattern[SS][xx - 1] == 1))
                continue;
              // CA -> CG
              if ((pattern[SS][xx] == 0) && (aa == 2) && (xx > 0)
                  && (pattern[SS][xx - 1] == 1))
                continue;
            }
            // Don't consider this pair if it can't qualify as a group.
            if (bicount[SS][x][a][xx][aa] < minCount)
              continue;

            for (xxx = 0; xxx <= xx - minDist; xxx++)
            {
              // Do not consider this site for mutation if disabled.
              if (useDisabledSites && disabledSites[xxx])
                continue;

              for (aaa = 0; aaa < 5; aaa++)
              {
                if (existingval[xxx][aaa] == 1)
                  continue;
                if ((DISALLOW_DASHDASH)
                    && (((a == 4) && (aa == 4) && (aaa == 4))
                        || (aa == 4 && aaa == 4)))
                  continue;
                if (DISALLOW_CG3)
                {
                  // CG -> TG
                  if ((pattern[SS][xxx] == 1) && (aaa == 3)
                      && (xxx < conLen - 1) && (pattern[SS][xxx + 1] == 2))
                    continue;
                  // TG -> CG
                  if ((pattern[SS][xxx] == 3) && (aaa == 1)
                      && (xxx < conLen - 1) && (pattern[SS][xxx + 1] == 2))
                    continue;
                  // CG -> CA 
                  if ((pattern[SS][xxx] == 2) && (aaa == 0)
                      && (xxx > 0) && (pattern[SS][xxx - 1] == 1))
                    continue;
                  // CA -> CG
                  if ((pattern[SS][xxx] == 0) && (aaa == 2)
                      && (xxx > 0) && (pattern[SS][xxx - 1] == 1))
                    continue;
                }
// TODO: Consider why we don't run EM to check that these three mutations
//       stay after rebuilding the local consensus pair....
//
                // Don't consider this pair if it can't qualify as a group.
                if (bicount[SS][xx][aa][xxx][aaa] < minCount ||
                    bicount[SS][x][a][xxx][aaa] < minCount)
                  continue;
                totalcount = assigncount[SS];
                count1 = count[SS][x][a];
                count2 = count[SS][xx][aa];
                count3 = count[SS][xxx][aaa];
                count12 = bicount[SS][x][a][xx][aa];
                count13 = bicount[SS][x][a][xxx][aaa];
                count23 = bicount[SS][xx][aa][xxx][aaa];
                // Expensive
                count123 = getTRICount(SS, x, a, xx, aa, xxx, aaa);
                triCountCalcs++;

                // New shortcut...why make extra calcs and proc calls?
                if (count123 < minCount)
                  continue;


                if (useSiegelPvalue)
                {
                  pvalue123 =
                    compute_siegel_tri_pvalue(totalcount,
                                              count1, count2,
                                              count3,
                                              count12,
                                              count13,
                                              count23,
                                              count123,
                                              dblEpsilon) + triBonferroni;
                }
                else
                {
                  fudge =
                    mutfrac[x][a] * ((double) totalcount) / ((double) count1);
                  fudge2 =
                    mutfrac[xx][aa] * ((double) totalcount) /
                    ((double) count2);
                  fudge3 =
                    mutfrac[xxx][aaa] *
                    ((double) totalcount) / ((double) count3);
                  if (fudge2 < fudge)
                    fudge = fudge2;
                  if (fudge3 < fudge)
                    fudge = fudge3;
                  pvalue123 =
                    compute_tri_pvalue(totalcount, count1,
                                       count2, count3,
                                       count12, count13,
                                       count23, count123,
                                       fudge) + triBonferroni;
                }

                // Add up the total real pValue calculations we did
                if (pvalue123 < 1000);
                pvalueCount++;

                if (pvalue123 >= bestmutpvalue)
                  continue;

                if (DEBUG >= 2)
                  printf("     \\--> new high pvalue = %lf\n", pvalue123);

                // Now separate into groups and recalculate consensus
                build_tri_local(x, a, xx, aa, xxx, aaa, SS);
                if (localpattern[1][x] != a)
                  continue;     // NO GOOD 
                if (localpattern[1][xx] != aa)
                  continue;     // NO GOOD 
                if (localpattern[1][xxx] != aaa)
                  continue;     // NO GOOD 

                if (DEBUG >= 2)
                  printf
                    ("           Passed cons recalc test .. Pos %d, "
                     "%c -> %c AND Pos %d, %c -> %c,"
                     " AND Pos %d, %c -> %c count123 = %d", x,
                     num_to_char(pattern[SS][x]),
                     num_to_char(a), xx,
                     num_to_char(pattern[SS][xx]),
                     num_to_char(aa), xxx,
                     num_to_char(pattern[SS][xxx]),
                     num_to_char(aaa), count123);

                bestmutpvalue = pvalue123;
                bestmutSS = SS;
                bestmutx = x;
                bestmuta = a;
                bestmutxx = xx;
                bestmutaa = aa;
                bestmutxxx = xxx;
                bestmutaaa = aaa;
              }
            }
          }
        }
      }
    }
  }
  printf("\n");
  if (DEBUG >= 1)
  {
    printf("  -- total pvalue calculations = %ld\n", pvalueCount);
    printf("  -- total tri count calculations = %ld\n", triCountCalcs);
  }
} // compute_tri_bestmut(...


//
// Search all scaffolds for most significant double mutation
//
void
compute_bestmut()
{
  int n, x, a, xx, aa, totalcount, SS, w, consensus, consensuscount, t;
  int xxx, aaa;
  double pvalue23, pvalue13, pvalue123;
  int count13, count23;
  int bestxx, bestaa;
  int count1, count2, count3, count12, s;
  double pvalue, pvalue2, pvaluelocal, fudge, fudge2;
  //int useOld = 0;

  long considered = 0;
  long pvalueCalcs = 0;

  if (DEBUG >= 3)
    printf("compute_bestmut(): Called\n");

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      existingval[x][a] = 0;
  }

  // Setting existing value
  for (s = 0; s < S; s++)
    for (x = 0; x < conLen; x++)
      existingval[x][pattern[s][x]] = 1;

if ( 0 ) 
{
  for (s = 0; s < S; s++)
  {
printf("Subfamily %d: ", s );
    for (x = 0; x < conLen; x++)
    {
if ( pattern[s][x] == 0 )
  printf("A");
else if ( pattern[s][x] == 1 )
  printf("C");
else if ( pattern[s][x] == 2 )
  printf("G");
else if ( pattern[s][x] == 3 )
  printf("T");
else 
  printf("-");
    }
printf("\n");
  }
}

  bestmutpvalue = PVALUETHRESH;
  printf("  Searching for significant double mutation:");
  for (SS = 0; SS < S; SS++)    /* try splitting subfamily SS */
  {
    considered = 0;
    pvalueCalcs = 0;

    if (DEBUG >= 2)
      printf("Trying subfamily %d..\n", SS);
    else
      printf(".");

    // Only needed for A. Price's method
    if (!useSiegelPvalue)
      compute_mutfrac(SS);

    for (x = minDist; x < conLen; x++)
    {
      // Do not consider this site for mutation if disabled.
      if (useDisabledSites && disabledSites[x])
        continue;

      for (a = 0; a < 5; a++)   /* try 1-mutation pos x val a */
      {
        if (DEBUG >= 3)
          printf(".... 1st Pos %d, %c -> %c : pattern %c%c%c\n",
                 x, num_to_char(pattern[SS][x]),
                 num_to_char(a),
                 num_to_char(pattern[SS][x - 1]),
                 num_to_char(pattern[SS][x]),
                 num_to_char(pattern[SS][x + 1]));

        //if ( useOld ) 
        //{
        if (existingval[x][a] == 1)
          continue;

        //}
        //else 
        //{
        //   if( pattern[SS][x] == a ||
        //       ( SS > 0 && pattern[parent[SS]][x] == a ) ) continue;
        //}

        if (count[SS][x][a] < minCount)
          continue;
        /*
           new code to disallow CG->TG etc 
         */
        if (DISALLOW_CG)
        {
          if ((pattern[SS][x] == 1) && (a == 3) && (x < conLen - 1) &&
              (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 3) && (a == 1) && (x < conLen - 1) &&
              (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 2) && (a == 0) && (x > 0) &&
              (pattern[SS][x - 1] == 1))
            continue;
          if ((pattern[SS][x] == 0) && (a == 2) && (x > 0) &&
              (pattern[SS][x - 1] == 1))
            continue;
        }
        for (xx = 0; xx <= x - minDist; xx++)
        {

          // Do not consider this site for mutation if disabled.
          if (useDisabledSites && disabledSites[xx])
            continue;

          for (aa = 0; aa < 5; aa++)
          {
            considered++;
            if (DEBUG >= 3)
              printf
                (".... 2nd Pos %d, %c -> %c : pattern %c%c%c\n", xx,
                 num_to_char(pattern[SS][xx]), num_to_char(aa),
                 num_to_char(pattern[SS][xx - 1]),
                 num_to_char(pattern[SS][xx]),
                 num_to_char(pattern[SS][xx + 1]));

            //if ( useOld ) 
            //{
            if (existingval[xx][aa] == 1)
              continue;
            //}
            //else 
            //{
            //   if( pattern[SS][xx] == aa ||
            //       ( SS > 0 && pattern[parent[SS]][xx] == aa ) ) continue;
            //}

            if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4))
              continue;
            /*
               new code to disallow CG->TG etc 
             */
            if (DISALLOW_CG2)
            {
              if ((pattern[SS][xx] == 1) && (aa == 3)
                  && (xx < conLen - 1) && (pattern[SS][xx + 1] == 2))
                continue;       /* CG -> TG */
              if ((pattern[SS][xx] == 3) && (aa == 1)
                  && (xx < conLen - 1) && (pattern[SS][xx + 1] == 2))
                continue;       /* TG -> CG */
              if ((pattern[SS][xx] == 2) && (aa == 0) && (xx > 0)
                  && (pattern[SS][xx - 1] == 1))
                continue;       /* CG -> CA */
              if ((pattern[SS][xx] == 0) && (aa == 2) && (xx > 0)
                  && (pattern[SS][xx - 1] == 1))
                continue;       /* CA -> CG */
            }
            totalcount = assigncount[SS];
            count1 = count[SS][x][a];
            count2 = count[SS][xx][aa];
            count12 = bicount[SS][x][a][xx][aa];
            if (count12 < minCount)
              continue;

            //printf("Computing Pvalue: n=%d, n1=%d, n2=%d, n12=%d\n", totalcount, count1, count2, count12 );
            if (useSiegelPvalue)
            {
              pvalue = compute_siegel_pvalue(totalcount, count1,
                                             count2, count12,
                                             dblEpsilon) + dblBonferroni;
            }
            else
            {
              fudge =
                mutfrac[x][a] * ((double) totalcount) / ((double) count1);
              fudge2 =
                mutfrac[xx][aa] * ((double) totalcount) / ((double) count2);
              if (fudge2 < fudge)
                fudge = fudge2;
              pvalue = compute_pvalue(totalcount, count1, count2,
                                      count12, fudge) + dblBonferroni;
            }
            pvalueCalcs++;

            if (DEBUG >= 3)
            {
              printf(".... Pos %d, %c -> %c AND Pos "
                     "%d, %c -> %c = %lf < %lf?\n", x,
                     num_to_char(pattern[SS][x]),
                     num_to_char(a), xx,
                     num_to_char(pattern[SS][xx]),
                     num_to_char(aa), pvalue, bestmutpvalue);
              printf
                ("     totalcount = %d, c1 = %d, c2 = %d, c12 = %d,"
                 " fudge = %lf\n", totalcount, count1, count2,
                 count12, fudge);
            }
            if (pvalue >= bestmutpvalue)
              continue;

            /*
             * Now, build local{count,assigncount,bicount,pattern,assign} 
             *
             * TODO: How does this differ from the build_subfamily EM 
             * algorithm?  Why does it resolve the two sites
             * here but often not at the end?
             *  
             *   - localpattern[0] original family consensus
             *   - localpattern[1] two diagnostic mutations
             *       - localpattern_to_localassign()
             *       - localasign_to_localpattern()
             */
            build_local(x, a, xx, aa, SS);      /* local{pattern,assign} */
            if (localpattern[1][x] != a)
              continue;         /* NO GOOD */
            if (localpattern[1][xx] != aa)
              continue;         /* NO GOOD */

            // RMH: Experimental
            // TODO: Optimize
            if ( S > 0 )
            {
              int status = build_local_global( x, a, xx, aa, SS );
              if ( status == 0 )
              {
                printf("Diagnostic mutations failed global check...continuing\n");
                continue;
              }
            }

            /*
               ALL NEW check pos x val a split against other subfamilies 
             */
            compute_localbicount(SS);   /* localcount and localbicount */
            for (s = 0; s < S; s++)
            {
              if (s == SS)
                continue;
              for (w = 0; w < 2; w++)   /* check local subfamily w U subfamily s */
              {
                pvaluelocal = split_pvaluelocal(SS, s, w, pvalue);
                if (pvaluelocal > pvalue)
                  pvalue = pvaluelocal;
                if (pvalue >= bestmutpvalue)
                  break;        /* inside w loop */
              }
              if (pvalue >= bestmutpvalue)
                break;          /* inside s loop */
            }
            if (pvalue >= bestmutpvalue)
              continue;         /* inside aa loop */
          
  
            /*
               WE HAVE A NEW WINNER! 
             */
            //printf("We have a winner: SS=%d Pos %d, %c -> %c AND Pos %d, %c -> %c :"
            //       " count12 = %d, totalcount=%d,  pvalue = %lf\n", SS, 
            //       x, num_to_char(pattern[SS][x]),num_to_char(a),
            //       xx, num_to_char(pattern[SS][xx]),num_to_char(aa), count12, totalcount, pvalue );
            //if ( x >= 3 && x < conLen - 3 && xx >= 3 && xx < conLen )
           // {
            //  printf("Parent1: %c%c%c %c %c%c%c\n", num_to_char(pattern[SS][x-3]), num_to_char(pattern[SS][x-2]), 
            //          num_to_char(pattern[SS][x-1]), num_to_char(pattern[SS][x]), num_to_char(pattern[SS][x+1]),
             //         num_to_char(pattern[SS][x+2]), num_to_char(pattern[SS][x+3]) );
             // printf("...avoid: a = %c = %d existingval[x][a] = %d\n", num_to_char(a), a, existingval[x][a]);
             // printf("Parent2: %c%c%c %c %c%c%c\n", num_to_char(pattern[SS][xx-3]), num_to_char(pattern[SS][xx-2]), 
             //         num_to_char(pattern[SS][xx-1]), num_to_char(pattern[SS][xx]), num_to_char(pattern[SS][xx+1]),
             //         num_to_char(pattern[SS][xx+2]), num_to_char(pattern[SS][xx+3]) );
              //printf("...avoid: aa = %c = %d existingval[xx][aa] = %d\n", num_to_char(aa), aa, existingval[xx][aa]);
           // }

            bestmutpvalue = pvalue;
            bestmutSS = SS;
            bestmutx = x;
            bestmuta = a;
            bestmutxx = xx;
            bestmutaa = aa;
          }
        }
      }
    }
    //printf(" --------> Considered = %ld, pvalueCalcs = %ld\n", considered, pvalueCalcs );
  }
  printf("\n");
}                               // void compute_bestmut() 



/* 
 *  Calculate mutfrac[SS][x][a]:
 *
 *  This matrix is useful for answering the question.  Given a 
 *  column 'x' and base 'a' what is the fraction of mutation ( in all
 *  positions x+/-mindist ) in all sequences containing an 'x' at 'a'.
 *
 *              mindist    mindist
 *              ----------v----------
 *cons ACCAACGTACGTACGATCGCTCGTACGTACGTAGCGCGGCGTACGA
 *
 *     ----C--------------A-------------------G--T---
 *     G-------T----------A-----A------T-------------
 *     ----G--------------A----------------G------G--
 *     -------------------A--T-----------------------
 *     G-----C------------A-------------C--------C---
 *cnt  2   2 1 1                       11  1  1  21   =  13 / 5 = 2.6
 *  
 *
 *     -------------------G-----------------C--------
 *     ----T----------G---G--T-----------------------
 *     -------------------G--------T-----------------
 *     -------------------G--------------------------
 *cnt      1                                1         = 2 / 4 = 2.0
 *                          
 *
 *                     No T mutations                = 1.0
 *    
 *
 *     -------------------C--------------------------
 *     -------------------C--------------------------
 *     -------------------C--------------------------
 *     ----T----------G---C--T-----------------------
 *     -------------------C--------------------------
 *     -------------------C--------------------------
 *cnt      1                                          = 1 / 6 = 0.17 
 *
 */
void
compute_mutfrac(int SS)
{
  int x, a, xx, aa;

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
    {
      numer[0][x][a] = 0;
      denom[0][x][a] = 0;
      for (xx = x + minDist; xx < conLen; xx++)
      {
        for (aa = 0; aa < 5; aa++)
        {
          if (pattern[SS][xx] == aa)
            continue;
          numer[0][x][a] += bicount[SS][xx][aa][x][a];  /* NOTE WELL */
          denom[0][x][a] += count[SS][xx][aa];
        }
      }
      for (xx = 0; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
        {
          if (pattern[SS][xx] == aa)
            continue;
          numer[0][x][a] += bicount[SS][x][a][xx][aa];
          denom[0][x][a] += count[SS][xx][aa];
        }
      }
      if (denom[0][x][a] == 0)
        mutfrac[x][a] = 1.0;
      else
        mutfrac[x][a] = ((double) numer[0][x][a]) / ((double) denom[0][x][a]);
    }
  }
}



void
compute_localbicount(int SS)    /* localcount and localbicount */
{
  int x, a, xx, aa, s, n;
  int *interesting, *dist2nextinteresting;
  int firstinteresting, firstinterestingminDist;

  interesting = malloc((conLen + 1) * sizeof(int));
  dist2nextinteresting = malloc((conLen + 1) * sizeof(int));


  for (x = 0; x < conLen; x++)
    interesting[x] = 0;
  for (s = 0; s < S; s++)
  {
    if (s == SS)
      continue;
    for (x = 0; x < conLen; x++)
    {
      if ((pattern[s][x] != localpattern[0][x]) ||
          (pattern[s][x] != localpattern[1][x]))
        interesting[x] = 1;
    }
  }
  for (x = 0; x < conLen; x++)
  {
    dist2nextinteresting[x] = conLen;
    for (xx = conLen - 1; xx > x; xx--)
    {
      if (interesting[xx] == 1)
        dist2nextinteresting[x] = xx - x;
    }
  }
  firstinteresting = 0 + dist2nextinteresting[0];
  firstinterestingminDist = minDist + dist2nextinteresting[minDist];

  for (x = firstinteresting; x < conLen; x += dist2nextinteresting[x])
  {
    for (a = 0; a < 5; a++)
    {
      localcount[0][x][a] = count[SS][x][a];
      localcount[1][x][a] = 0;
    }
  }
  for (x = firstinterestingminDist; x < conLen; x += dist2nextinteresting[x])
  {
    for (a = 0; a < 5; a++)
    {
      for (xx = firstinteresting; xx <= x - minDist;
           xx += dist2nextinteresting[xx])
      {
        for (aa = 0; aa < 5; aa++)
        {
          localbicount[0][x][a][xx][aa] = bicount[SS][x][a][xx][aa];
          localbicount[1][x][a][xx][aa] = 0;
        }
      }
    }
  }

  for (n = 0; n < N; n++)
  {
    if ((assign[n] == SS) && (localassign[n] == 1))     /* 0 -> 1 */
    {
      for (x = firstinteresting; x < conLen; x += dist2nextinteresting[x])
      {
        localcount[0][x][ele[n][x]]--;
        localcount[1][x][ele[n][x]]++;
      }
      for (x = firstinterestingminDist; x < conLen;
           x += dist2nextinteresting[x])
      {
        for (xx = firstinteresting; xx <= x - minDist;
             xx += dist2nextinteresting[xx])
        {
          localbicount[0][x][ele[n][x]][xx][ele[n][xx]]--;
          localbicount[1][x][ele[n][x]][xx][ele[n][xx]]++;
        }
      }
    }
  }
  free( interesting );
  free( dist2nextinteresting );
}


// 
// A comparator for use with qsort.  Compares
// two elements subfamily assignments.
//
int
intcmp(const void *v1, const void *v2)
{
  return (assign[*(int *) v1] - assign[*(int *) v2]);
}


//
// Get a count of all elements in subfamily SS which contain
// the co-occurance of mutation x:a, xx:aa, and xxx:aaa.   
//
int
getTRICount(int SS, int x, int a, int xx, int aa, int xxx, int aaa)
{

  int n = 0;
  int i = 0;
  int endIdx;
  int count = 0;

  // Lookup the end of the subfamily range within the
  // sortedAssign list.
  if (SS + 1 == S)
    endIdx = maxN;
  else
    endIdx = sortedAssignIndices[SS + 1];
  // Iterating over the subfamily block within the 
  // sortedAssign list count all co-occurances of
  // mutations.
  for (i = sortedAssignIndices[SS]; i < endIdx; i++)
  {
    n = sortedAssign[i];
    if (ele[n][x] == a && ele[n][xx] == aa && ele[n][xxx] == aaa)
      count++;
  }

  return (count);
}



double
compute_tri_pvalue(int totalcount, int count1, int count2,
                   int count3, int count12, int count13,
                   int count23, int count123, double fudge)
{
  int t, countmin;
  double logdenom, answer, logterm, fudge2, fudge3, logcorrection;
  int tt, c1, c2, c12;
  double expanswer;
  double score;

  //printf(" totalcount = %d, count1 = %d, count2=%d, count3=%d, count12=%d, count13=%d, count23=%d, count123=%d\n", totalcount, count1, count2, count3, count12, count13, count23, count123 );

  // This is a very important optimisation
  if ((score =
       compute_tri_score(totalcount, count1, count2, count3, count123)) <
      SCORETHRESH)
  {
    return 1000.0;
  }

  answer = 0.0;
  expanswer = 0.0;
  // Triple factor
  if (fudge < 1.0)
    fudge = 1.0;
  else if (fudge > 2.0)
    fudge = 2.0;
  fudge2 =
    1.0 - (fudge -
           1.0) * ((double) count1) / ((double) (totalcount - count1));
  // Hmmm
  fudge3 =
    1.0 - (fudge -
           1.0) * ((double) count2) / ((double) (totalcount - count2));

  logdenom = logfactorial[totalcount] - logfactorial[count1]
    - logfactorial[totalcount - count1]
    + logfactorial[totalcount] - logfactorial[count2]
    - logfactorial[totalcount - count2]
    + logfactorial[totalcount] - logfactorial[count3]
    - logfactorial[totalcount - count3];

  //printf("logdenom = %lf\n", logdenom );
  countmin = count1;
  if (count2 < countmin)
    countmin = count2;
  if (count3 < countmin)
    countmin = count3;
  if (count12 < countmin)
    countmin = count12;
  if (count13 < countmin)
    countmin = count13;
  if (count23 < countmin)
    countmin = count23;


  //printf("count123 = %d to countmin = %d or %d iterations\n", count123, countmin, (countmin - count123) );

  t = count123;
  {
    logterm = logfactorial[totalcount]
      - logfactorial[t]
      - logfactorial[count3 - count23 - count13 + t]
      - logfactorial[count2 - count12 - count23 + t]
      - logfactorial[count1 - count13 - count12 + t]
      - logfactorial[count23 - t]
      - logfactorial[count13 - t]
      - logfactorial[count12 - t]
      - logfactorial[totalcount - count1 - count2 - count3 + count13 +
                     count23 + count12 - t];

    if (logterm > logdenom)
      return 1000.0;
    answer = logterm - logdenom;
  }

  for (t = count123 + 1; t <= countmin; t++)
  {
    logterm = logfactorial[totalcount]
      - logfactorial[t]
      - logfactorial[count3 - count23 - count13 + t]
      - logfactorial[count2 - count12 - count23 + t]
      - logfactorial[count1 - count13 - count12 + t]
      - logfactorial[count23 - t]
      - logfactorial[count13 - t]
      - logfactorial[count12 - t]
      - logfactorial[totalcount - count1 - count2 - count3 + count13 +
                     count23 + count12 - t];

    //logcorrection = ((double)t)*log(fudge) + ((double)(count2-t))*log(fudge2) +
    //                ((double)(count3-t))*log(fudge3);
    //if(logcorrection > 0) logterm += logcorrection;
    if (logterm > logdenom)
      return 1000.0;
    answer += exp(logterm - logdenom - answer);
  }

  return answer;
}

double
compute_siegel_tri_pvalue(int n, int n1, int n2, int n3,
                          int n12, int n13, int n23, int n123, double epsilon)
{
  double p1_3sigma = 0;
  double p2_3sigma = 0;
  double p3_3sigma = 0;
  double mTerm = 0;
  double cTerm = 0;
  double twoThirds;
  double oneEighteenth;
  double numerator = 0;
  double denominator = 0;
  int nMin;
  int nMax;
  double minVal = 0;
  double maxVal = 0;
  double logPValue = 0;
  int m = 0;
  int i = 0;

  double mean, sdev, actual;

  if (n123 < minCount || n1 == 0 || n2 == 0 || n3 == 0)
  {
    return 1000.0;
  }
  else
  {
    mean =
      ((double) n1) * ((double) n2) * ((double) n3) / (((double) n) *
                                                       ((double) n) *
                                                       ((double) n));
    // binomial standard deviation
    sdev = mean * (1.0 - mean) / ((double) n);
    sdev = exp(0.5 * log(sdev));
    actual = ((double) n123) / ((double) n);
    if (((actual - mean) / sdev) < SCORETHRESH)
    {
      return 1000.0;
    }
  }

  // Conservative estimates of p1, p2, and p3 may be obtained by using 
  // the following three-sigma upper bounds instead:
  p1_3sigma = (n1 / n) + ((double) 3 * sqrt(((n1 / n) * (1 - n1 / n)) / n));
  p2_3sigma = (n2 / n) + ((double) 3 * sqrt(((n2 / n) * (1 - n2 / n)) / n));
  p3_3sigma = (n3 / n) + ((double) 3 * sqrt(((n3 / n) * (1 - n3 / n)) / n));

  twoThirds = (double) 2 / (double) 3;
  oneEighteenth = (double) 1 / (double) 18;

  // constant portion of numerator/denominator
  cTerm = ((((twoThirds - p2_3sigma) *
             (twoThirds - p3_3sigma)) + oneEighteenth) *
           (((twoThirds - p1_3sigma) *
             (twoThirds - p3_3sigma)) + oneEighteenth) *
           (((twoThirds - p1_3sigma) *
             (twoThirds - p2_3sigma)) + oneEighteenth)) /
    (((twoThirds - p3_3sigma) *
      (twoThirds - p2_3sigma) *
      (twoThirds - p1_3sigma)) *
     (((twoThirds - p1_3sigma) *
       (twoThirds - p2_3sigma) *
       (twoThirds - p3_3sigma)) -
      ((p1_3sigma + p2_3sigma + p3_3sigma) / (double) 18) +
      ((double) 11 / (double) 54)));

  // NEW: nMin = max(0,n12+n13-n1,n12+n23-n2,n13+n23-n3)
  nMin = 0;
  if (nMin < (n12 + n13 - n1))
    nMin = (n12 + n13 - n1);
  if (nMin < (n12 + n23 - n2))
    nMin = (n12 + n23 - n2);
  if (nMin < (n13 + n23 - n3))
    nMin = (n13 + n23 - n3);

  // NEW: nMax = min( n12, n13, n23, n-n1-n2-n3+n12+n13+n23-n123 )
  nMax = n12;
  if (nMax > n13)
    nMax = n13;
  if (nMax > n23)
    nMax = n23;
  if (nMax > (n - n1 - n2 - n3 + n12 + n13 + n23 - n123))
    nMax = (n - n1 - n2 - n3 + n12 + n13 + n23 - n123);

  // Given: We have a lookup table logfact where
  //        logfactorial[n] = ln( n! )
  // 
  // Calculate numerator of pValue:
  for (m = n123; m <= nMax; m++)
  {
    // combinatorial portion ( in logspace )
    mTerm = (logfactorial[n]
             - (logfactorial[m]
                + logfactorial[n12 - m]
                + logfactorial[n13 - m]
                + logfactorial[n23 - m]
                + logfactorial[n1 - n12 - n13 + m]
                + logfactorial[n2 - n12 - n23 + m]
                + logfactorial[n3 - n13 - n23 + m]
                + logfactorial[n - n1 - n2 - n3 + n12 + n13 + n23 - m]));

    // constant portion ( cTerm^m in logspace )
    mTerm += m * log(cTerm);

    // NOTE: We are summing numbers that have already been converted
    //       to logspace.
    //
    // Andy's algorithm:
    //    Given x1 = ln(p1), x2 = ln(p2)
    //    Then ln(x1+x2) = max(x1,x2) + ln{ 1 + exp(min(x1,x2)-max(x1,x2)) }
    //    This should be computed directly if exp(x1-x2) does not cause
    //    underflow.  If underflow is observed replace exp(x1-x2) with zero.
    //
    // Make sure we do not have underflow ( parameter epsilon holds
    // the smallest floating point number we can use safely -- to
    // be used as a threshold for the calculation ).
    //
    // Calculate numerator += mTerm in logspace:
    minVal = numerator;
    maxVal = numerator;
    if (minVal > mTerm)
      minVal = mTerm;
    if (maxVal < mTerm)
      maxVal = mTerm;

    //printf("tri numerator: n1=%d n2=%d n3=%d, n12=%d n13=%d n23=%d n123=%d n=%d\n", n1, n2, n3, n12, n13, n23, n123, n );
    //printf("tri numerator: Below epsilon (%lf - %lf) = %le epsilon = %le\n", minVal, maxVal, fabs(minVal - maxVal), epsilon );
    //printf("tri numerator: %lf\n", numerator );
    //printf("newcalc: %lf\n", log( 1 + exp( minVal - maxVal )  ));
    if (fabs(minVal - maxVal) > epsilon)
      numerator = maxVal + log(1 + exp(minVal - maxVal));
    else
    {
      numerator = maxVal;
    }
  }

  // Calculate the denominator similarly
  for (m = nMin; m <= nMax; m++)
  {
    // combinatorial portion ( in logspace )
    mTerm = (logfactorial[n]
             - (logfactorial[m]
                + logfactorial[n12 - m]
                + logfactorial[n13 - m]
                + logfactorial[n23 - m]
                + logfactorial[n1 - n12 - n13 + m]
                + logfactorial[n2 - n12 - n23 + m]
                + logfactorial[n3 - n13 - n23 + m]
                + logfactorial[n - n1 - n2 - n3 + n12 + n13 + n23 - m]));

    // constant portion ( cTerm^m in logspace )
    mTerm += m * log(cTerm);

    // Calculate denominator += mTerm in logspace:
    minVal = denominator;
    maxVal = denominator;
    if (minVal > mTerm)
      minVal = mTerm;
    if (maxVal < mTerm)
      maxVal = mTerm;

    if (fabs(minVal - maxVal) > epsilon)
      denominator = maxVal + log(1 + exp(minVal - maxVal));
    else
    {
      denominator = maxVal;
      //printf("tri denominator: Below epsilon (%lf - %lf) = %lf\n", minVal, maxVal, fabs(minVal - maxVal) );
      //exit(1);
    }
  }

  logPValue = numerator - denominator;

  return logPValue;
}



//
// Test of a pvalue function from Andy Siegel
//
double
compute_siegel_pvalue(int n, int n1, int n2, int n12, double epsilon)
{
  double p1_3sigma = 0;
  double p2_3sigma = 0;
  double mTerm = 0;
  double cTerm = 0;
  double numerator = 0;
  double denominator = 0;
  int minN = n1;
  double minVal = 0;
  double maxVal = 0;
  double logPValue = 0;
  int m = 0;
  int i = 0;

  double mean, sdev, actual;

  if (n12 < minCount || n1 == 0 || n2 == 0)
  {
    //printf("Returning early counts\n");
    return 1000.0;
  }
  else
  {
    mean = ((double) n1) * ((double) n2) / (((double) n) * ((double) n));
    sdev = mean * (1.0 - mean) / ((double) n);
    sdev = exp(0.5 * log(sdev));
    actual = ((double) n12) / ((double) n);
    if (((actual - mean) / sdev) < SCORETHRESH)
    {
      //printf("Returning early stddev\n");
      return 1000.0;
    }
  }

  // Conservative estimates of p1 and p2 may be obtained by using 
  // the following three-sigma upper bounds instead:
  p1_3sigma = (n1 / n) + (3 * sqrt(((n1 / n) * (1 - n1 / n)) / n));
  p2_3sigma = (n2 / n) + (3 * sqrt(((n2 / n) * (1 - n2 / n)) / n));

  // constant portion of numerator/denominator
  cTerm =
    (4 / 3 * (1 - p1_3sigma - p2_3sigma + (4 / 3) * p1_3sigma * p2_3sigma)) /
    ((1 - 4 / 3 * p1_3sigma) * (1 - 4 / 3 * p2_3sigma));

  // Determine min(n1,n2)
  if (minN > n2)
    minN = n2;

  // Given: We have a lookup table logfact where
  //        logfactorial[n] = ln( n! )
  // 
  // Calculate numerator of pValue:
  for (m = n12; m <= minN; m++)
  {
    // combinatorial portion ( in logspace )
    mTerm = (logfactorial[n]
             - (logfactorial[m]
                + logfactorial[n1 - m]
                + logfactorial[n2 - m] + logfactorial[n - n1 - n2 + m]));
    // constant portion ( cTerm^m in logspace )
    mTerm += m * log(cTerm);


    // NOTE: We are summing numbers that have already been converted
    //       to logspace.
    //
    // Andy's algorithm:
    //    Given x1 = ln(p1), x2 = ln(p2)
    //    Then ln(x1+x2) = max(x1,x2) + ln{ 1 + exp(min(x1,x2)-max(x1,x2)) }
    //    This should be computed directly if exp(x1-x2) does not cause
    //    underflow.  If underflow is observed replace exp(x1-x2) with zero.
    //
    // Make sure we do not have underflow ( parameter epsilon holds
    // the smallest floating point number we can use safely -- to
    // be used as a threshold for the calculation ).
    //
    // Calculate numerator += mTerm in logspace:
    minVal = numerator;
    maxVal = numerator;
    if (minVal > mTerm)
      minVal = mTerm;
    if (maxVal < mTerm)
      maxVal = mTerm;

    if (fabs(minVal - maxVal) > epsilon)
      numerator = maxVal + log(1 + exp(minVal - maxVal));
    else
    {
      numerator = maxVal;
      printf("************ Below epsilon (%lf - %lf) = %lf\n", minVal,
             maxVal, fabs(minVal - maxVal));
    }
  }

  // Calculate the denominator similarly
  for (m = 0; m <= minN; m++)
  {
    // combinatorial portion ( in logspace )
    mTerm = (logfactorial[n]
             - (logfactorial[m]
                + logfactorial[n1 - m]
                + logfactorial[n2 - m] + logfactorial[n - n1 - n2 + m]));
    // constant portion ( cTerm^m in logspace )
    mTerm += m * log(cTerm);

    // Calculate denominator += mTerm in logspace:
    minVal = denominator;
    maxVal = denominator;
    if (minVal > mTerm)
      minVal = mTerm;
    if (maxVal < mTerm)
      maxVal = mTerm;

    if (fabs(minVal - maxVal) > epsilon)
      denominator = maxVal + log(1 + exp(minVal - maxVal));
    else
    {
      denominator = maxVal;
      printf("Below epsilon (%lf - %lf) = %lf\n", minVal, maxVal,
             fabs(minVal - maxVal));
    }


  }

  logPValue = numerator - denominator;

  return logPValue;

}


double
getDoubleEpsilon()
{
  double epsilon = 1;
  do
  {
    epsilon = epsilon / (double) 2;
  }
  while ((double) 1 + epsilon > (double) 1);
  epsilon = (double) 2 *epsilon;

  return epsilon;
}

double
compute_pvalue(int totalcount, int count1, int count2, int count12,
               double fudge)
{
  int t, countmin;
  double logdenom, answer, logterm, fudge2, logcorrection;
  int tt, c1, c2, c12;


  if (compute_score(totalcount, count1, count2, count12) < SCORETHRESH)
  {
    return 1000.0;
  }

  answer = 0.0;
  if (fudge < 1.0)
    fudge = 1.0;
  else if (fudge > 1.333333)
    fudge = 4.0 / 3.0;
  fudge2 =
    1.0 - (fudge -
           1.0) * ((double) count1) / ((double) (totalcount - count1));

  logdenom =
    logfactorial[totalcount] - logfactorial[count1] -
    logfactorial[totalcount - count1] + logfactorial[totalcount] -
    logfactorial[count2] - logfactorial[totalcount - count2];
  countmin = count1;
  if (count2 < countmin)
    countmin = count2;
  t = count12;
  {
    logterm =
      logfactorial[totalcount] - logfactorial[t] - logfactorial[count1 - t] -
      logfactorial[count2 - t] - logfactorial[totalcount -
                                              (count1 + count2 - t)];
    logcorrection =
      ((double) t) * log(fudge) + ((double) (count2 - t)) * log(fudge2);
    if (logcorrection > 0)
      logterm += logcorrection;
    if (logterm > logdenom)
      return 1000.0;
    answer = logterm - logdenom;
  }
  for (t = count12 + 1; t <= countmin; t++)
  {
    logterm =
      logfactorial[totalcount] - logfactorial[t] - logfactorial[count1 -
                                                                t] -
      logfactorial[count2 - t] - logfactorial[totalcount -
                                              (count1 + count2 - t)];
    logcorrection =
      ((double) t) * log(fudge) + ((double) (count2 - t)) * log(fudge2);
    if (logcorrection > 0)
      logterm += logcorrection;
    if (logterm > logdenom)
      return 1000.0;
    answer += exp(logterm - logdenom - answer);
  }

  return answer;
}



double
compute_score(int totalcount, int count1, int count2, int count12)
{
  double mean, sdev, actual;

  if (count12 < minCount)
    return 0.0;

  mean =
    ((double) count1) * ((double) count2) / (((double) totalcount) *
                                             ((double) totalcount));
  if (mean == 0.0)
    return -1000000.0;          /* count1 or count2 is 0 */
  sdev = mean * (1.0 - mean) / ((double) totalcount);
  sdev = exp(0.5 * log(sdev));
  actual = ((double) count12) / ((double) totalcount);

  return (actual - mean) / sdev;
}



double
compute_tri_score(int totalcount, int count1, int count2, int count3,
                  int count123)
{
  double mean, sdev, actual;

  if (count123 < minCount)
    return 0.0;

  mean = ((double) count1) * ((double) count2) * ((double) count3) /
    (((double) totalcount) * ((double) totalcount) * ((double) totalcount));
  if (mean == 0.0)
    return -1000000.0;          /* count1 or count2 or count3 is 0 */
  sdev = mean * (1.0 - mean) / ((double) totalcount);
  sdev = exp(0.5 * log(sdev));
  actual = ((double) count123) / ((double) totalcount);

  //printf("compute_tri_score = %lf\n", (actual - mean)/sdev );

  return (actual - mean) / sdev;
}



void
build_new_tri_subfamily()       /* parent, S, assign, pattern_to_assign */
{
  int n, s, x, a, xx, aa, xxx, aaa, SS, besta, bestb, b;
  int ccount[5], ccounti[2];
  int oldcount, newcount;
  double bestscore, countratio;

  x = bestmutx;
  a = bestmuta;
  xx = bestmutxx;
  aa = bestmutaa;
  xxx = bestmutxxx;
  aaa = bestmutaaa;

  /*
   *  we are given bestmutSS, bestmutpos, bestmutval 
   */
  SS = bestmutSS;
  parent[S] = SS;

  build_tri_local(x, a, xx, aa, xxx, aaa, SS);

  fprintf(outFP, "Building subfamily %d (parent %d, logpvalue %e): "
          "pos %d %c to %c and pos %d %c to %c and pos %d %c to %c\n",
          S, SS, bestmutpvalue, x, num_to_char(pattern[SS][x]),
          num_to_char(a), xx, num_to_char(pattern[SS][xx]),
          num_to_char(aa), xxx, num_to_char(pattern[SS][xxx]),
          num_to_char(aaa));
  fprintf(stdout,
          "  Building subfamily %d (parent %d, logpvalue %e): "
          "pos %d %c to %c and pos %d %c to %c and pos %d %c to %c\n", S, SS,
          bestmutpvalue, x, num_to_char(pattern[SS][x]), num_to_char(a),
          xx, num_to_char(pattern[SS][xx]), num_to_char(aa), xxx,
          num_to_char(pattern[SS][xxx]), num_to_char(aaa));

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      count[S][x][a] = 0;
    // RMH Added 11/8/12
    for (a = 0; a < 2; a++)
      counti[S][x][a] = 0;
  }
  for (x = minDist; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
    {
      for (xx = 0; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
          bicount[S][x][a][xx][aa] = 0;
      }
    }
  }

// Does nothing......for EM 
if ( 0 ) 
{
  // Initialy breakup up subfamily SS into SS & S.
  for (n = 0; n < N; n++)
  {
    if ((assign[n] == SS) && (localassign[n] == 1))
    {
      for (x = 0; x < conLen; x++)
      {
        count[S][x][ele[n][x]]++;
        count[SS][x][ele[n][x]]--;
        // RMH Added 11/8/12
        counti[S][x][elei[n][x]]++;
        counti[SS][x][elei[n][x]]--;
      }
      for (x = minDist; x < conLen; x++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
        {
          bicount[S][x][ele[n][x]][xx][ele[n][xx]]++;
          bicount[SS][x][ele[n][x]][xx][ele[n][xx]]--;
        }
      }
      assign[n] = S;
    }
  }
  assigncount[SS] = localassigncount[0];
  assigncount[S] = localassigncount[1];
}

  for (x = 0; x < conLen; x++)
  {
    pattern[SS][x] = localpattern[0][x];
    pattern[S][x] = localpattern[1][x];
    patterni[SS][x] = localpatterni[0][x];
    patterni[S][x] = localpatterni[1][x];
  }
  S++;

  // This calls the EM routine which works directly
  // on pattern/patterni[][] and assign[].  It 
  // adjusts all elements and all consensi.
  run_em();


if ( 0 ) 
{
  printf("Validating Initial Diagnostic Sites: %d: ",bestmutx);
  if(pattern[S-1][bestmutx] != bestmuta) 
    printf("No  ");
  else
    printf("Yes ");
  printf("%d: ", bestmutxx);
  if(pattern[S-1][bestmutxx] != bestmutaa) 
    printf("No\n");
  else
    printf("Yes\n");
  printf("%d: ", bestmutxxx);
  if(pattern[S-1][bestmutxxx] != bestmutaaa) 
    printf("No\n");
  else
    printf("Yes\n");

  printf("Validation of diagnostic distribution:\n");
  printf("   - pos=%d: ", bestmutx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), count[S-1][bestmutx][x] );
  printf("\n");
  printf("   - pos=%d: ", bestmutxx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), count[S-1][bestmutxx][x] );
  printf("\n");
  printf("   - pos=%d: ", bestmutxxx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), count[S-1][bestmutxxx][x] );
  printf("\n");
}
 
} // build_new_tri_subfamily


//
// build_new_subfamily()
//   Create a datastructure to hold new subfamily count/bicounts
//   Move ( subtract/add ) counts from old subfamily to new datastructure
//   Initialize new subfamily pattern datastructure
//
//  Modifies:
//      count[], counti[], and bicount[]
void
build_new_subfamily()           /* parent, S, assign, pattern_to_assign */
{
  int n, s, x, a, xx, aa, SS, besta, bestb, b;
  int ccount[5], ccounti[2];
  int oldcount, newcount;
  int totalAssigned = 0;
  double bestscore, countratio;

  x = bestmutx;
  a = bestmuta;
  xx = bestmutxx;
  aa = bestmutaa;

  /*
     we are given bestmutSS, bestmutpos, bestmutval 
   */
  SS = bestmutSS;
  parent[S] = SS;
  build_local(x, a, xx, aa, SS);        /* local{pattern,assign} */

  fprintf(outFP,
          "Building subfamily %d (parent %d, logpvalue %e): pos %d %c to %c and pos %d %c to %c\n",
          S, SS, bestmutpvalue, x, num_to_char(pattern[SS][x]),
          num_to_char(a), xx, num_to_char(pattern[SS][xx]), num_to_char(aa));
  fprintf(stdout,
          "  Building subfamily %d (parent %d, size %d, logpvalue %e): pos %d %c to %c and pos %d %c to %c\n",
          S, SS, assigncount[SS], bestmutpvalue, x,
          num_to_char(pattern[SS][x]), num_to_char(a), xx,
          num_to_char(pattern[SS][xx]), num_to_char(aa));


  // 
  // Very noisy debuging
  //
  if (VERBOSE & 4)
  {
    fprintf(stdout,
            "Breaking subfamily %d using: pos %d %c to %c and pos %d %c to %c\n",
            SS, bestmutx, num_to_char(pattern[SS][bestmutx]),
            num_to_char(bestmuta), bestmutxx,
            num_to_char(pattern[SS][bestmutxx]), num_to_char(bestmutaa));
    fprintf(stdout, "Current Cons: ");
    for (x = 0; x < conLen; x++)
    {
      if (x == bestmutx || x == bestmutxx)
        fprintf(stdout, "%c", toupper(num_to_char(pattern[SS][x])));
      else
        fprintf(stdout, "%c", num_to_char(pattern[SS][x]));
    }
    fprintf(stdout, "\nStarting Sequences: ");
    for (n = 0; n < N; n++)
    {
      if (assign[n] == SS)
        fprintf(stdout, "%d,", n);
    }
    fprintf(stdout, "\nInitial Family Sizes: ");
    for (s = 0; s < S; s++)
      fprintf(stdout, "%d,", assigncount[s]);
    fprintf(stdout, "\n");
  }

// Special DEBUG
if ( 0 ) 
{
    fprintf(stdout, "Current Cons: ");
    for (x = 0; x < conLen; x++)
    {
      if (x == bestmutx || x == bestmutxx)
        if ( pattern[SS][x] == 4 )
          fprintf(stdout, "^");
        else
          fprintf(stdout, "%c", toupper(num_to_char(pattern[SS][x])));
      else
        fprintf(stdout, "%c", num_to_char(pattern[SS][x]));
    }
    int numWithMut = 0;
    int numTotal = 0;
    for (n = 0; n < N; n++)
    {
      if (assign[n] == SS)
      {
        numTotal++;
        if ( ele[n][bestmutx] == bestmuta && ele[n][bestmutxx] == bestmutaa )
          numWithMut++;
      }
    }
    fprintf(stdout, "Parent size %d subset with both muts %d,", numTotal, numWithMut );
}
// End special debug


  // Initialize the new subfamily count/bicount datastructure
  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      count[S][x][a] = 0;
    for (a = 0; a < 2; a++ )
      counti[S][x][a] = 0;
  }
  for (x = minDist; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
    {
      for (xx = 0; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
          bicount[S][x][a][xx][aa] = 0;
      }
    }
  }

  // Subtract counts/bicounts for elements moving over 
  // to the new subfamily.
  for (n = 0; n < N; n++)
  {
    // Guesing that localassign is set when the element
    // should move over to the new subfamily.
    if ((assign[n] == SS) && (localassign[n] == 1))
    {
if ( 0 ) {
      // If so...move it's counts
      for (x = 0; x < conLen; x++)
      {
        count[S][x][ele[n][x]]++;
        count[SS][x][ele[n][x]]--;
        // RMH added 11/8/12
        counti[S][x][elei[n][x]]++;
        counti[SS][x][elei[n][x]]--;
      }
      // If so...move it's bicounts
      for (x = minDist; x < conLen; x++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
        {
          bicount[S][x][ele[n][x]][xx][ele[n][xx]]++;
          bicount[SS][x][ele[n][x]][xx][ele[n][xx]]--;
        }
      }
      assign[n] = S;
}
      totalAssigned++;
    }
  }
//  assigncount[SS] = localassigncount[0];
//  assigncount[S] = localassigncount[1];

  //printf("build_new_subfamily(): totalAssigned = %d assigncount[%d] = %d, assigncount[%d] = %d\n", totalAssigned, SS, assigncount[SS], S, assigncount[S] );
  //printf("build_new_subfamily(): totalAssigned = %d assigncount[%d] = %d, assigncount[%d] = %d\n", totalAssigned, SS, localassigncount[0], S, localassigncount[1] );

  for (x = 0; x < conLen; x++)
  {
    pattern[SS][x] = localpattern[0][x];
    pattern[S][x] = localpattern[1][x];
    patterni[SS][x] = localpatterni[0][x];
    patterni[S][x] = localpatterni[1][x];
  }
  S++;

  // Now pattern[] and assign[] both are up-to-date
  //
  // Reconsider the assignment of all sequences to subfamily consensus
  // patterns 0-S, redevelop cosensus subfamily consensus sequences 0-S,
  // reconsider the assignment of all sequences.......basically cycle
  // until everything settles down.
  run_em();


if ( 0 ) 
{
  printf("Validating Initial Diagnostic Sites: %d: ",bestmutx);
  if(pattern[S-1][bestmutx] != bestmuta) 
    printf("No  ");
  else
    printf("Yes ");
  printf("%d: ", bestmutxx);
  if(pattern[S-1][bestmutxx] != bestmutaa) 
    printf("No\n");
  else
    printf("Yes\n");
  printf("Validation of diagnostic distribution:\n");
  printf("   - pos=%d: ", bestmutx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), count[S-1][bestmutx][x] );
  printf("\n");
  printf("   - pos=%d: ", bestmutxx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), count[S-1][bestmutxx][x] );
  printf("\n");
}
 
  if (VERBOSE > 4)
  {
    fprintf(stdout, "\nFinal Family Sizes: ");
    for (s = 0; s < S; s++)
      fprintf(stdout, "%d,", assigncount[s]);
    fprintf(stdout, "\n");
  }


}

/*
 * EM Algorithm for to-date built subfamily consensi
 *
 * This EM run is done after a subfamily has been 
 * added to the set and involves all subfamilies in
 * the competition.  This means that the new subfamily
 * just added can steal members from other subfamilies
 * other than it's parent.
 *
 * NOTE: After this is run there is currently no check to
 * see if the last added subfamily did maintain the two 
 * co-segregating sites.
 */
void
run_em()
{
  int n, s, x, SS;
  int iter;

  SS = bestmutSS;

  // TODO: Hmmmm...why why why? I wish I had documented this!
  // RMH added 11/7/2012
  //
  // Take element assignments and rebuild consensi stored
  // in pattern[][] and patterni[][]
  //
  //  This shouldn't do anything because assignments and pattern 
  //  are up-to-date.
  //assign_to_pattern();          /* M-step */

  if (VERBOSE & 4)
  {
    fprintf(stdout, "EM Initial Step %d\n", iter);
    fprintf(stdout, "Cluster 0 Sequences: ");
    for (n = 0; n < N; n++)
    {
      if (assign[n] == SS)
        fprintf(stdout, "%d,", n);
    }
    fprintf(stdout, "\nCluster 1 Sequences: ");
    for (n = 0; n < N; n++)
    {
      if (assign[n] == (S - 1))
        fprintf(stdout, "%d,", n);
    }
    fprintf(stdout, "\n");
  }

  totdist = 1000000000;
  for (iter = 0; iter < MAXITER; iter++)
  {
    // Using the patterns determine which pattern
    // each element is closest to ( edit distance ) and
    // alter the assign[] array accordingly.
    pattern_to_assign();        /* E-step */

    if (VERBOSE & 4)
    {
      fprintf(stdout, "EM Optimisation Step %d\n", iter);
      fprintf(stdout, "Cluster 0 Sequences: ");
      for (n = 0; n < N; n++)
      {
        if (assign[n] == SS)
          fprintf(stdout, "%d,", n);
      }
      fprintf(stdout, "\nCluster 1 Sequences: ");
      for (n = 0; n < N; n++)
      {
        if (assign[n] == (S - 1))
          fprintf(stdout, "%d,", n);
      }
      fprintf(stdout, "\n");
    }

    // TODO: Doesn't this lead to the counts being out-of-step with the patterns?  Should we roll-back the assignment changes?
    //       Or...probably Alkes considered the pattern authoritative and the assignment should be derived whenever needed.
    if (totdist >= oldtotdist)
      break;

    assign_to_pattern();        /* M-step */

    if (VERBOSE & 4)
    {
      fprintf(stdout, "Cluster 0 Cons: ");
      for (x = 0; x < conLen; x++)
      {
        if (pattern[SS][x] != pattern[S - 1][x])
          fprintf(stdout, "%c", toupper(num_to_char(pattern[SS][x])));
        else
          fprintf(stdout, "%c", num_to_char(pattern[SS][x]));
      }
      fprintf(stdout, "\n   patternI: ");
      for (x = 0; x < conLen; x++)
        fprintf(stdout, "%d", patterni[SS][x]);
      fprintf(stdout, "\nCluster 1 Cons: ");
      for (x = 0; x < conLen; x++)
      {
        if (pattern[SS][x] != pattern[S - 1][x])
          fprintf(stdout, "%c", toupper(num_to_char(pattern[S - 1][x])));
        else
          fprintf(stdout, "%c", num_to_char(pattern[S - 1][x]));
      }
      fprintf(stdout, "\n");
      fprintf(stdout, "   patternI: ");
      for (x = 0; x < conLen; x++)
        fprintf(stdout, "%d", patterni[S - 1][x]);
      fprintf(stdout, "\n");
    }

    if (VERBOSE & 1)
    {
      fprintf(outFP, "S=%d run_em iter %d, totdist=%d\n", S, iter, totdist);
      printf("S=%d run_em iter %d, totdist=%d\n", S, iter, totdist);
    }
  }
  // RMH: 11/2/2012 : Shouldn't the pattern now reflect the last step
  //                  of reassignment???  Adding this
  //TODO: RMH:....um..no
  //assign_to_pattern();          /* M-step */

  if (VERBOSE & 4)
  {
    fprintf(stdout, "Final 0 Cons: ");
    for (x = 0; x < conLen; x++)
    {
      if (pattern[SS][x] != pattern[S - 1][x])
        fprintf(stdout, "%c", toupper(num_to_char(pattern[SS][x])));
      else
        fprintf(stdout, "%c", num_to_char(pattern[SS][x]));
    }
    fprintf(stdout, "\nFinal 1 Cons: ");
    for (x = 0; x < conLen; x++)
    {
      if (pattern[SS][x] != pattern[S - 1][x])
        fprintf(stdout, "%c", toupper(num_to_char(pattern[S - 1][x])));
      else
        fprintf(stdout, "%c", num_to_char(pattern[S - 1][x]));
    }
    fprintf(stdout, "\n");
  }
}


/* 
 * M step for EM algorithm
 *
 * For all subfamilies rebuild consensus based on 
 * current assignments.  
 * 
 * Uses: count[s][x][a]
 *       pattern[s][x]
 */
void
assign_to_pattern()             /* M-step */
{
  int n, s, x, a, b, maxcount;

  for (s = 0; s < S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      maxcount = -1;
      /*
         a,c,g,t max 
       */
      for (a = 0; a < 5; a++)
      {
        if (count[s][x][a] > maxcount)
        {
          maxcount = count[s][x][a];
          pattern[s][x] = a;
        }
      }

      maxcount = -1;
      /*
         0 = base, 1 = insert max 
       */
      for (b = 0; b < 2; b++)
      {
        if (counti[s][x][b] > maxcount)
        {
          maxcount = counti[s][x][b];
          patterni[s][x] = b;
        }
      }
    }
  }
  // The above is identical to the singlmut routines 
  // singlemut then goes on to calc many more things.
}


/*
 * E step of EM algorithm
 *
 * Alter element assignments based on edit distance from
 * any of the ( possibly altered ) subfamily consensi.
 */
void
pattern_to_assign()             /* E-step */
{
  int n, s, x, a, xx, aa, b, mindist, dist, ok, oldassign, newassign;

  oldtotdist = totdist;
  totdist = 0;
  for (n = 0; n < N; n++)
  {
    oldassign = assign[n];
    ok = 0;
    mindist = 1000000000;
    for (s = 0; s < S; s++)
    {
      dist = 0;
      for (x = 0; x < conLen; x++)
      {
        a = ele[n][x];
        if (a != pattern[s][x])
          dist++;
        b = elei[n][x];
        if (b != patterni[s][x])
          dist++;
      }
      /*
         Subtract off CpG mutations from dist 
       */
      if (DISALLOW_CG)
      {
        for (x = 0; x < conLen - 1; x++)
        {
          if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 2)) /* CpG */
          {
            if (ele[n][x] == 3)
              dist--;
            else if (ele[n][x + 1] == 0)
              dist--;
          }
          else if ((pattern[s][x] == 3) && (pattern[s][x + 1] == 2))    /* TpG */
          {
            if (ele[n][x] == 1)
              dist--;
          }
          else if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 0))    /* TpG */
          {
            if (ele[n][x + 1] == 2)
              dist--;
          }
        }
      }
      if (dist < mindist)
      {
        mindist = dist;
        assign[n] = s;
        ok = 1;
      }
      else if (dist == mindist)
        ok = 0;
    }
    totdist += mindist;

    newassign = assign[n];
    // TODO: Adjust count/counti here so we don't have to recalculate
    //       the entire matrix in assign_to_pattern
    if (newassign != oldassign)
    {
      assigncount[oldassign]--;
      for (x = 0; x < conLen; x++)
      {
        count[oldassign][x][ele[n][x]]--;
        // RMH Added 11/8/12
        counti[oldassign][x][elei[n][x]]--;
      }
      for (x = minDist; x < conLen; x++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
          bicount[oldassign][x][ele[n][x]][xx][ele[n][xx]]--;
      }

      assigncount[newassign]++;
      for (x = 0; x < conLen; x++)
      {
        count[newassign][x][ele[n][x]]++;
        // RMH Added 11/8/12
        counti[newassign][x][elei[n][x]]++;
      }
      for (x = minDist; x < conLen; x++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
          bicount[newassign][x][ele[n][x]][xx][ele[n][xx]]++;
      }
    }
  }
}

void
buildCpGAdjustedConsensus()
{

  int n, s, x, k, dnScore, cgScore;

  int CGParam = 14;
  int TAParam = -5;

  char const linupMatrix[] = {
    9, -8, -15, -17, -6,
    -4, 10, -14, -15, -6,
    -15, -14, 10, -4, -6,
    -16, -15, 1, 2, -6,
    -11, -2, -11, -3, -6,
    -6, -6, -6, -6, 3
  };
  /*
     if(c == 'a') return 0;
     if(c == 'c') return 1;
     if(c == 'g') return 2;
     if(c == 't') return 3;
     if(c == '-') return 4;
     if(c == 'N') return 99;
     A    G   C    T   -
   */

  // Load in insertions file
  // and build maxleni[s][x] = len
  //int n = 0;
  //readline{
  //  p = pos;
  //  l = len( );
  //  if ( maxleni[assign[n]][p] < l )
  //    maxleni[assign[n]][p] = l;
  //  n++;
  //}
  //seek back to begining

  for (s = 0; s < S; s++)
  {

    // Calculate the total cons len
    // int totalConsLen = conLen;
    // for ( x = 0; x < conLen; x++ )
    //   totalConsLen += maxleni[s][x];
    // malloc conscount array
    //
    // int ccIdx = 0;
    // for ( n = 0; n < N; n++ )
    // { 
    //   int inAligned = 0;
    //   for (x = 0; x < conLen; x++)
    //   {
    //     // add ele[n][x] to counts
    //     if ( elei[n][x] == 1 )
    //     {
    //       // seek to row of file
    //     }
    //   }
    // }


    for (x = 0; x < conLen - 2; x++)
    {
      // Go through the derived consensus and consider changing
      // each dinucleotide to a 'CG'
      cgScore = 0;
      dnScore = 0;

      // Ignore gaps
      if (pattern[s][x] == 4)
        continue;

      // Skip adacent gap characters
      k = x + 1;
      while (k < conLen && pattern[s][k] == 4)
        k++;

      // Now x = left di-nucleotide base and k = right di-nucleotide base
      for (n = 0; n < N; n++)
      {
        // Is this sequence in this subfamily?
        if (assign[n] == s)
        {
          // Scored as-is
          dnScore += linupMatrix[(pattern[s][x] * 4) +
                                 ele[n][x]] +
            linupMatrix[(pattern[s][k] * 4) + ele[n][k]];

          // Ele di-nucleotide is "CA" or "TG"
          if ((ele[n][x] == 1 && ele[n][k] == 0) ||
              (ele[n][x] == 3 && ele[n][k] == 2))
          {
            cgScore += CGParam;
            // Ele di nucleotide is "TA"
          }
          else if (ele[n][x] == 3 && ele[n][k] == 0)
          {
            cgScore += TAParam;
            // Ele di nucleotide is "TC" or "TT"
          }
          else if (ele[n][x] == 3 && (ele[n][k] == 1 || ele[n][k] == 3))
          {
            // in other words; c->t transition scores +2
            // transversion scored normally.
            cgScore += 2 + linupMatrix[8 + ele[n][k]];
            // Ele di nucleotide is "AA" or "GA"
          }
          else if ((ele[n][x] == 0 || ele[n][x] == 2) && ele[n][k] == 0)
          {
            // same as above
            cgScore += 2 + linupMatrix[4 + ele[n][x]];
          }
          else
          {
            cgScore += linupMatrix[4 + ele[n][x]] +
              linupMatrix[8 + ele[n][k]];
          }
          if (cgScore > dnScore)
          {
            // Pattern should be "C*G"
            pattern[s][x] = 1;
            pattern[s][k] = 2;
          }
        }                       // if ( assign[n] == s...
      }                         // for ( n = 0; n < N...
    }                           // for ( x = 0; x<conLen...
  }                             // for ( s = 0; s<S...


  // Create a DS indicating if a "s" has one or more insertions in patterni
  // so....like inserts[s] = [ { int pos,
  //                             char * seqStruct == null }, ... ]
  //                           
  //malloc inserts[S]
  // set all inserts[0-S] to null;
  // for each patterni[s][0-conLen]
  //   if ( patterni[s][x] == 1 )
  //     add to inserts[s] ( malloc etc )

  // open insert file
  //   for each element n
  //     if ( inserts[ assign[ n ] ] != null )
  //     {
  //       foreach insert position in inserts[assign[n]][]->pos
  //         add sequence to inserts[assign[n][]->seqStruct[].seq/count 
  //     }
  //  close insert file

  // Go through all insert_seqs[s][x].keys and pick highest value
  //

}                               // void buildCpGAdjustedConsensus()




void
build_tri_local(int thisx, int thisa, int thisxx, int thisaa,
                int thisxxx, int thisaaa, int SS)
{
  int x, iter;


  /*
     Step 1: build localpattern 
   */
  for (x = 0; x < conLen; x++)
  {
    localpattern[0][x] = pattern[SS][x];
    localpatterni[0][x] = patterni[SS][x];
    localpattern[1][x] = pattern[SS][x];
    localpatterni[1][x] = patterni[SS][x];
  }

  // Create initial subfamily separation
  localpattern[1][thisx] = thisa;
  localpattern[1][thisxx] = thisaa;
  localpattern[1][thisxxx] = thisaaa;

  /*
     NOW, run local em 
   */
  localtotdist = 1000000000;
  localpattern_to_localassign(SS);      /* E-step */
  for (iter = 0; iter < MAXITER; iter++)
  {
    //printf( "build_local:  EM Iteration %d\n", iter );
    localassign_to_localpattern(SS);    /* M-step */
    localpattern_to_localassign(SS);    /* E-step */
    if (localtotdist >= localoldtotdist)
      break;
  }
}

// EXPERIMENTAL
void blc_pattern_to_assign()
{
  int mindist, n, s, x, a, b, dist;
  // blc_pattern_to_assign
  oldtotdist = totdist;
  totdist = 0;
  for (n = 0; n < N; n++)
  {
    mindist = 1000000000;
    // s = S is the new subfamily being considered
    for (s = 0; s <= S; s++)
    {
      dist = 0;
      for (x = 0; x < conLen; x++)
      {
        a = ele[n][x];
        if ( a != l_pattern[s][x] )
          dist++;
        b = elei[n][x];
        if ( b != l_patterni[s][x] )
          dist++;
      }
      // Subtract off CpG mutations from dist 
      if (DISALLOW_CG)
      {
        for (x = 0; x < conLen - 1; x++)
        { 
          // CpG
          if ((l_pattern[s][x] == 1) && (l_pattern[s][x + 1] == 2))
          { 
            if (ele[n][x] == 3)
              dist--;
            else if (ele[n][x + 1] == 0)
              dist--;
          }
          // TpG
          else if ((l_pattern[s][x] == 3) && (l_pattern[s][x + 1] == 2)) 
          { 
            if (ele[n][x] == 1)
              dist--;
          }
          // TpG
          else if ((l_pattern[s][x] == 1) && (l_pattern[s][x + 1] == 0))
          { 
            if (ele[n][x + 1] == 2)
              dist--;
          }
        }
      }
//printf("Distance ele=%d S=%d: %d\n", n, s, dist);
      if (dist < mindist)
      {
        mindist = dist;
        l_assign[n] = s;
      }
    } // for ( s = 0
    totdist += mindist;
  } // for( n = 0
}

// EXPERIMENTAL
void
blc_assign_to_pattern(int thisx, int thisa, int thisxx, int thisaa, int SS)
{
  int maxcount,x, n, s, a, b;

  int ac = 0;

  // s=S is the new subfamily being proposed.
  for (s = 0; s <= S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
        l_count[s][x][a] = 0;
      for (b = 0; b < 2; b++)
        l_counti[s][x][b] = 0;
    }
  }

  // Count each cluster's columns
  for (n = 0; n < N; n++)
  {
    s=l_assign[n];
    if ( s == S )
      ac++;
    for (x = 0; x < conLen; x++)
    {
      a = ele[n][x];
      b = elei[n][x];
      l_count[s][x][a]++;
      l_counti[s][x][b]++;
    }
  }
/*
printf("Assigned count = %d out of N = %d \n",ac, N);
printf("Updated diagnostic distribution:\n");
printf("   - pos=%d: ", thisx );
for ( x = 0; x < 5; x++ )
  printf("%c=%d, ", num_to_char(x), l_count[S][thisx][x] );
printf("\n");
printf("   - pos=%d: ", thisxx );
for ( x = 0; x < 5; x++ )
  printf("%c=%d, ", num_to_char(x), l_count[S][thisxx][x] );
printf("\n");
*/
 
  for (s = 0; s <= S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      maxcount = -1;
      for (a = 0; a < 5; a++)
      {
        if (l_count[s][x][a] > maxcount)
        {
          maxcount = l_count[s][x][a];
          l_pattern[s][x] = a;
        }
      }

      maxcount = -1;
      for (b = 0; b < 2; b++)
      {
        if (l_counti[s][x][b] > maxcount)
        {
          maxcount = l_counti[s][x][b];
          l_patterni[s][x] = b;
        }
      }
    }
  }
}



// 
// Experimental: The existing method for picking a significant
//               diagnostic site set includes a test to see if
//               the sites remain in the final consensus after
//               running EM on the members of the parent subfamily.
//               After that process the code commits to those sites
//               and runs EM over all elements and all subfamilies.
//               In this EM implementation we also check ( before
//               committing ) that the new subfamily consensus
//               keeps the diagnostic sites after doing a global
//               pattern_to_assign/assign_to_pattern ( ie. one round
//               of EM ) on only the new consensus and by not altering
//               the actual global state datastructures.
//
int
build_local_global(int thisx, int thisa, int thisxx, int thisaa, int SS)
{
  int n, s, w, x, a, b, mindist, dist;
  int iter;
  int maxcount = 0;


  // Assume that localpattern/i[1] contains the proposed
  // subfamily consensus and localpattern/i[0] contains
  // the revised parent family consensus.
  for (s = 0; s <= S; s++)
  {
    //if ( s == S )
    //  printf("Proposed Cons:\n");
    //else
    //  printf("sub %d:\n", s);
    for (x = 0; x < conLen; x++)
    {
      if ( s == S )
      {
        l_pattern[s][x] = localpattern[1][x];
        l_patterni[s][x] = localpatterni[1][x];
        //printf("%c",num_to_char(l_pattern[s][x]));
      }else if ( s == SS )
      {
        l_pattern[s][x] = localpattern[0][x];
        l_patterni[s][x] = localpatterni[0][x];
        //printf("%c",num_to_char(l_pattern[s][x]));
      }else
      {
        l_pattern[s][x] = pattern[s][x];
        l_patterni[s][x] = patterni[s][x];
        //printf("%c",num_to_char(l_pattern[s][x]));
      }
    }
    //printf("\n");
  }

  // Clear the l_assign[]
  for (n = 0; n < N; n++)
    l_assign[n] = 0;

  totdist = 100000000;
  blc_pattern_to_assign();

  for (iter = 0; iter < MAXITER; iter++)
  {
    blc_assign_to_pattern(thisx, thisa, thisxx, thisaa, SS);
    blc_pattern_to_assign();
    if (totdist >= oldtotdist)
        break;
  }

  //printf("iterations performed = %d\n", iter );

  /*
  if (l_pattern[S][thisx] != thisa)
    printf("Missing site %d is %c should be %c\n", thisx, num_to_char(l_pattern[S][thisx]), num_to_char(thisa));
  if (l_pattern[S][thisxx] != thisaa)
    printf("Missing site %d is %c should be %c\n", thisxx, num_to_char(l_pattern[S][thisxx]), num_to_char(thisaa));
  */

  /*
  printf("Final: ");
  for (x = 0; x < conLen; x++)
    printf("%c",num_to_char(l_pattern[S][x]) );
  printf("\n");
  */
  
  int retVal = 0;
  if (l_pattern[S][thisx] == thisa && l_pattern[S][thisxx] == thisaa)
    retVal = 1;

  return( retVal );
}




/* local{pattern,assign} */
void
build_local(int thisx, int thisa, int thisxx, int thisaa, int SS)
{
  int x, iter;

//printf( "****Starting build_local***** : Seeding split with %d:%c and "
//       "%d:%c\n", thisx, num_to_char(thisa), 
//      thisxx, num_to_char(thisaa) );

  /*
     Step 1: build localpattern 
   */
  for (x = 0; x < conLen; x++)
  {
    localpattern[0][x] = pattern[SS][x];
    localpatterni[0][x] = patterni[SS][x];
    localpattern[1][x] = pattern[SS][x];
    localpatterni[1][x] = patterni[SS][x];
  }

  //printf( "build_local: Seeding split with %d:%c and "
  //       "%d:%c\n", thisx, num_to_char(thisa), 
  //      thisxx, num_to_char(thisaa) );

  // Create initial subfamily separation
  localpattern[1][thisx] = thisa;
  localpattern[1][thisxx] = thisaa;

  /*
     NOW, run local em 
   */
  localtotdist = 1000000000;
  localpattern_to_localassign(SS);      /* E-step */
  for (iter = 0; iter < MAXITER; iter++)
  {
    //printf( "build_local:  EM Iteration %d\n", iter );
    localassign_to_localpattern(SS);    /* M-step */
    localpattern_to_localassign(SS);    /* E-step */
    if (localtotdist >= localoldtotdist)
      break;
  }
  //printf("build_local: em iter reached %d out of max %d\n",iter, MAXITER);
  //printf( "build_local: Final %d:%c and "
  //        "%d:%c\n", thisx, num_to_char(localpattern[1][thisx]), 
  //        thisxx, num_to_char(localpattern[1][thisxx]) );

  /*
  int s;
  for (s = 0; s < 2; s++)
  {
    printf("localpattern[%d] count = %d\n", s, localassigncount[s] );
    for (x = 0; x < conLen; x++)
      printf("%c",num_to_char(localpattern[s][x]));
    printf("\n");
  }
  
  printf("Diagnostic distribution:\n");
  printf("   - pos = %d: ", thisx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), localcount[1][thisx][x] );
  printf("\n");
  printf("   - pos = %d: ", thisxx );
  for ( x = 0; x < 5; x++ )
    printf("%c=%d, ", num_to_char(x), localcount[1][thisxx][x] );
  printf("\n");
  */
 

}



//
// localpattern_to_localassign( int SS )
//
//      SS : The parent (source) of the new subfamily
//   
//      Expectation function of the EM algorithm.  This takes
//      two consensus patterns and assigns all elements of the
//      parent subfamily SS to one or the other pattern based
//      on a distance score.  The score is +1/-1 (match/mismatch/deletion)
//      and -1 for insertion ( any length ).  CpG sites are ignored.
//
// Globals Used:
//       ele[0-N][0-conLen]
//       elei[0-N][0-conLen]
//       assign[0-N]
//       localpattern[0/1]
//       localassign[0/1]
//       localassigncount[0/1]
//       localoldtotdist
//       localtotdist
//
// Used By: build_local()
//
void
localpattern_to_localassign(int SS)
{
  int n, w, x, a, b, mindist, dist;

  localoldtotdist = localtotdist;
  localtotdist = 0;

  // reset the assignment DS
  for (w = 0; w < 2; w++)
    localassigncount[w] = 0;
  for (n = 0; n < N; n++)
  {
    if (assign[n] != SS)
      continue;
    mindist = 1000000000;
    for (w = 0; w < 2; w++)
    {
      dist = 0;
      for (x = 0; x < conLen; x++)
      {
        a = ele[n][x];
        if (a != localpattern[w][x])
          dist++;
        b = elei[n][x];
        if (b != localpatterni[w][x])
          dist++;
      }
      // Subtract off CpG mutations from dist 
      if (DISALLOW_CG)
      {
        for (x = 0; x < conLen - 1; x++)
        {
          if ((localpattern[w][x] == 1) && (localpattern[w][x + 1] == 2))       // CpG 
          {
            if (ele[n][x] == 3)
              dist--;
            else if (ele[n][x + 1] == 0)
              dist--;
          }
          else if ((localpattern[w][x] == 3) && (localpattern[w][x + 1] == 2))  // TpG 
          {
            if (ele[n][x] == 1)
              dist--;
          }
          else if ((localpattern[w][x] == 1) && (localpattern[w][x + 1] == 0))  // TpG 
          {
            if (ele[n][x + 1] == 2)
              dist--;
          }
        }
      }
      if (dist < mindist)
      {
        mindist = dist;
        localassign[n] = w;
      }
    }
    localassigncount[localassign[n]] += 1;
    localtotdist += mindist;
  }
}                               // void localpattern_to_localassign( ...



//
// localassign_to_localpattern( int SS )
//
//      SS : The parent (source) of the new subfamily
//   
//      Maximization function of the EM algorithm.  This takes
//      the assignment of SS subfamliy members ( into 2 groups )
//      and recalculates the consensus patterns for each group.
//
// Globals Used:
//       ele[0-N][0-conLen]
//       elei[0-N][0-conLen]
//       assign[0-N]
//       localpattern[0/1]
//       localassign[0/1]
//       localassigncount[0/1]
//
// Used By: build_local()
//
void
localassign_to_localpattern(int SS)
{
  int n, w, x, a, b, maxcount;

  // Clear the counts
  for (w = 0; w < 2; w++)
  {
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
        localcount[w][x][a] = 0;
      for (b = 0; b < 2; b++)
        localcounti[w][x][b] = 0;
    }
  }

  // Count each cluster's columns
  for (n = 0; n < N; n++)
  {
    if (assign[n] != SS)
      continue;
    w = localassign[n];
    for (x = 0; x < conLen; x++)
    {
      a = ele[n][x];
      b = elei[n][x];
      localcount[w][x][a]++;
      localcounti[w][x][b]++;
    }
  }

  for (w = 0; w < 2; w++)
  {
    for (x = 0; x < conLen; x++)
    {
      maxcount = -1;
      for (a = 0; a < 5; a++)
      {
        if (localcount[w][x][a] > maxcount)
        {
          maxcount = localcount[w][x][a];
          localpattern[w][x] = a;
        }
      }

      maxcount = -1;
      for (b = 0; b < 2; b++)
      {
        if (localcounti[w][x][b] > maxcount)
        {
          maxcount = localcounti[w][x][b];
          localpatterni[w][x] = b;
        }
      }
    }
  }
}                               // void localassign_to_localpattern(... 

char
char_to_num(char c)
{
  if (c == 'A')
    return 0;
  if (c == 'C')
    return 1;
  if (c == 'G')
    return 2;
  if (c == 'T')
    return 3;
  if (c == 'a')
    return 0;
  if (c == 'c')
    return 1;
  if (c == 'g')
    return 2;
  if (c == 't')
    return 3;
  if (c == '-')
    return 4;
  if (c == 'N')
    return 99;
  if (c == 'n')
    return 99;
  if (c == 'B')
    return 99;
  if (c == 'b')
    return 99;
  if (c == 'D')
    return 99;
  if (c == 'd')
    return 99;
  if (c == 'H')
    return 99;
  if (c == 'h')
    return 99;
  if (c == 'V')
    return 99;
  if (c == 'v')
    return 99;
  if (c == 'R')
    return 99;
  if (c == 'r')
    return 99;
  if (c == 'K')
    return 99;
  if (c == 'k')
    return 99;
  if (c == 'M')
    return 99;
  if (c == 'm')
    return 99;
  if (c == 'S')
    return 99;
  if (c == 's')
    return 99;
  if (c == 'W')
    return 99;
  if (c == 'w')
    return 99;
  printf("char_to_num(): Invalid DNA character c=%d=%c\n", c, c);
  printf("               Check your sequence file for non DNA letters.\n");
  printf
    ("               NOTE: The sequence file is generated by preprocessAlighments.pl.\n");
  printf("                     and is not in FASTA format.\n");
  exit(1);
}


char
num_to_char(char z)
{
  if (z == 0)
    return (char) 'a';
  if (z == 1)
    return (char) 'c';
  if (z == 2)
    return (char) 'g';
  if (z == 3)
    return (char) 't';
  if (z == 4)
    return (char) '-';
  return (char) 'N';
}


double
split_pvaluelocal(int SS, int s, int w, double pvaluehope)
{
  int x, a, olda, xx, aa, oldaa, totalcount, count1, count2, count12, flip;
  char *consensus;
  double answer, pvalue, fudge;
  double localbonferroni;

  consensus = malloc((conLen + 1) * sizeof(char));

  for (x = 0; x < conLen; x++)
  {
    // Do not consider this site for mutation if disabled.
    if (useDisabledSites && disabledSites[x])
      continue;
    olda = pattern[s][x];
    a = localpattern[w][x];
    if ((olda == a)
        || (count[s][x][olda] + localcount[w][x][olda] >=
            count[s][x][a] + localcount[w][x][a]))
      consensus[x] = olda;
    else
      consensus[x] = a;
  }

  localbonferroni = 0.0;
  for (x = minDist; x < conLen; x++)
  {
    // Do not consider this site for mutation if disabled.
    if (useDisabledSites && disabledSites[x])
      continue;
    if (pattern[s][x] != localpattern[w][x])
      continue;
    for (xx = 0; xx <= x - minDist; xx++)
    {
      // Do not consider this site for mutation if disabled.
      if (useDisabledSites && disabledSites[xx])
        continue;
      if (pattern[s][xx] != localpattern[w][xx])
        continue;
      localbonferroni += 1.0;
    }
  }
  if (localbonferroni == 0.0)
    return 0.0;
  localbonferroni = log(localbonferroni);

  answer = 0.0;
  for (x = minDist; x < conLen; x++)    /* hope pattern[s][x] vs. localpattern[w][x] */
  {
    // Do not consider this site for mutation if disabled.
    if (useDisabledSites && disabledSites[x])
      continue;
    if (pattern[s][x] == localpattern[w][x])
      continue;
    olda = consensus[x];
    a = pattern[s][x];
    if (a == olda)
      a = localpattern[w][x];
    if (count[s][x][a] + localcount[w][x][a] < minCount)
      continue;
    /*
       new code to disallow CG->TG etc 
     */
    if (DISALLOW_CG)
    {
      if ((olda == 1) && (a == 3) && (x < conLen - 1)
          && (consensus[x + 1] == 2))
        continue;               /* CG -> TG */
      if ((olda == 3) && (a == 1) && (x < conLen - 1)
          && (consensus[x + 1] == 2))
        continue;               /* TG -> CG */
      if ((olda == 2) && (a == 0) && (x > 0) && (consensus[x - 1] == 1))
        continue;               /* CG -> CA */
      if ((olda == 0) && (a == 2) && (x > 0) && (consensus[x - 1] == 1))
        continue;               /* CA -> CG */
    }
    for (xx = 0; xx <= x - minDist; xx++)
    {
      // Do not consider this site for mutation if disabled.
      if (useDisabledSites && disabledSites[xx])
        continue;
      if (pattern[s][xx] == localpattern[w][xx])
        continue;
      oldaa = consensus[xx];
      aa = pattern[s][xx];
      if (aa == oldaa)
        aa = localpattern[w][xx];
      if (count[s][xx][aa] + localcount[w][xx][aa] < minCount)
        continue;
      if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4))
        continue;
      /*
         new code to disallow CG->TG etc 
       */
      if (DISALLOW_CG2)
      {
        if ((oldaa == 1) && (aa == 3) && (xx < conLen - 1)
            && (consensus[xx + 1] == 2))
          continue;             /* CG -> TG */
        if ((oldaa == 3) && (aa == 1) && (xx < conLen - 1)
            && (consensus[xx + 1] == 2))
          continue;             /* TG -> CG */
        if ((oldaa == 2) && (aa == 0) && (xx > 0) && (consensus[xx - 1] == 1))
          continue;             /* CG -> CA */
        if ((oldaa == 0) && (aa == 2) && (xx > 0) && (consensus[xx - 1] == 1))
          continue;             /* CA -> CG */
      }
      totalcount = assigncount[s] + localassigncount[w];
      count1 = count[s][x][a] + localcount[w][x][a];
      count2 = count[s][xx][aa] + localcount[w][xx][aa];
      count12 = bicount[s][x][a][xx][aa] + localbicount[w][x][a][xx][aa];
      if (count12 == 0)
        continue;
      if (useSiegelPvalue)
      {
        pvalue = compute_siegel_pvalue(totalcount, count1, count2,
                                       count12, dblEpsilon) + localbonferroni;
      }
      else
      {
        pvalue =
          compute_pvalue(totalcount, count1, count2, count12,
                         1.0) + localbonferroni;
      }
      //printf( "**pvalue = %lf\n", pvalue);

      if (pvalue <= pvaluehope)
        return pvaluehope;
      if (pvalue < answer)
        answer = pvalue;
    }
  }

  free ( consensus );

  return answer;
}


void
print_assign(char *filename)
{
  int n;
  FILE *fp;

  if ((fp = fopen(filename, "w")) == NULL)
  {
    printf("Could not open input file %s\n", filename);
    exit(1);
  }

  for (n = 0; n < N; n++)
    fprintf(fp, "%d %d\n", n, assign[n]);

  fclose(fp);
}

// prune_subfamilies()
//
// Input: assigncount[]
// Output: distance[][] 
void
prune_subfamilies()
{
  int s, t, s2, thisdist, maxdist, label1, label2, STREE, sunite, tunite,
    mindist;
  double pvalue;
  FILE *fp;

  // Run EM just to get everything squared away
  run_em();

  // This calculates numer/denom for thefactor mutfrac[] which is used only in 
  // A. Prices method.
  if (!useSiegelPvalue)
    for (s = 0; s < S; s++)
      compute_numerdenom(s);

  maxdist = 0;
  /*
     Compute pairwise distances, merge immediately if <2 or
     * if any subfamily now exhibits less than minCount membership 
   */
  for (s = 0; s < S; s++)
  {
    if (assigncount[s] < minCount)
    {
      t = parent[s];
      for (s2 = 0; s2 < S; s2++)
        labels[s2] = 0;
      labels[s] = 1;
      labels[t] = 1;
      fprintf(outFP,
              "Merging subfamilies %d and %d because subfamily %d is less than mincount.\n",
              s, t, s);
      printf
        ("Merging subfamilies %d and %d because subfamily %d is less than mincount.\n",
         s, t, s);
      merge_subfamilies(1);
      prune_subfamilies();
    }
    for (t = s + 1; t < S; t++)
    {
      distance[s][t] = compute_distance(s, t, MULTIPLE_MUTATION_INSPENALTY);
      distance[t][s] = distance[s][t];
      if (distance[s][t] > maxdist)
        maxdist = distance[s][t];
      if (distance[s][t] < 2)
      {
        for (s2 = 0; s2 < S; s2++)
          labels[s2] = 0;
        labels[s] = 1;
        labels[t] = 1;
        fprintf(outFP,
                "Merging subfamilies %d and %d because dist=%d\n", s,
                t, distance[s][t]);
        printf("Merging subfamilies %d and %d because dist=%d\n", s, t,
               distance[s][t]);
        merge_subfamilies(1);
        prune_subfamilies();
        return;
      }
    }
  }

  /*
     evaluate/merge pairwise in order of ascending distance 
   */
  for (thisdist = 2; thisdist <= maxdist; thisdist++)
  {
    for (s = 0; s < S; s++)
    {
      for (t = s + 1; t < S; t++)
      {
        if (distance[s][t] != thisdist)
          continue;
        for (s2 = 0; s2 < S; s2++)
          labels[s2] = 0;
        labels[s] = 1;
        labels[t] = 1;
        if (useTRI)
        {
          pvalue = union_tri_pvalue(1);
        }
        else
        {
          pvalue = union_pvalue(1);
        }
        //printf("UNION PVALUE %d to %d = %lf\n", s, t, pvalue );
        if (pvalue <= PVALUETHRESH)
          continue;

        fprintf(outFP,
                "Merging subfamilies %d and %d (dist=%d) because "
                "pvalue=%f\n", s, t, distance[s][t], pvalue);
        printf("Merging subfamilies %d and %d (dist=%d) because "
               "pvalue=%f\n", s, t, distance[s][t], pvalue);
        merge_subfamilies(1);
        prune_subfamilies();
        return;
      }
    }
  }

  // Appears to go with all pairs in ascending order.
  /*
     agglomerative evaluate/merge (Kruskal) 
   */
  for (s = 0; s < S; s++)
    labels[s] = s;
  STREE = S;
  while (STREE > 1)
  {
    /*
       unite closest s and t with different labels 
     */
    mindist = 1000;
    for (s = 0; s < S; s++)
    {
      for (t = s + 1; t < S; t++)
      {
        if (labels[s] == labels[t])
          continue;
        if (distance[s][t] < mindist)
        {
          sunite = s;
          tunite = t;
          mindist = distance[s][t];
        }
      }
    }

    if (labels[sunite] < labels[tunite])
    {
      label1 = labels[sunite];
      label2 = labels[tunite];
    }
    else
    {
      label1 = labels[tunite];
      label2 = labels[sunite];
    }

    for (s = 0; s < S; s++)
    {
      if (labels[s] == label2)
        labels[s] = label1;
    }

    if (useTRI)
    {
      pvalue = union_tri_pvalue(label1);
    }
    else
    {
      pvalue = union_pvalue(label1);
    }
    STREE--;
    if (pvalue <= PVALUETHRESH)
      continue;

    fprintf(outFP, "Merging subfamilies ");
    for (s = 0; s < S; s++)
    {
      if (labels[s] == label1)
        fprintf(outFP, "%d ", s);
    }
    fprintf(outFP, "because pvalue=%f\n", pvalue);
    printf("Merging subfamilies ");
    for (s = 0; s < S; s++)
    {
      if (labels[s] == label1)
        printf("%d ", s);
    }
    printf("because pvalue=%f\n", pvalue);
    merge_subfamilies(label1);
    prune_subfamilies();
    return;
  }
} // prune_subfamilies()

void
compute_numerdenom(int s)
{
  int x, a, xx, aa;

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
    {
      numer[s][x][a] = 0;
      denom[s][x][a] = 0;
      for (xx = x + minDist; xx < conLen; xx++)
      {
        for (aa = 0; aa < 5; aa++)
        {
          if (pattern[s][xx] == aa)
            continue;
          numer[s][x][a] += bicount[s][xx][aa][x][a];   /* NOTE WELL */
          denom[s][x][a] += count[s][xx][aa];
        }
      }
      for (xx = 0; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
        {
          if (pattern[s][xx] == aa)
            continue;
          numer[s][x][a] += bicount[s][x][a][xx][aa];
          denom[s][x][a] += count[s][xx][aa];
        }
      }
    }
  }
}

double
union_tri_pvalue(int label)
{
  int s, x, a, xx, aa, xxx, aaa, ok, totnumer1, totdenom1, totnumer2,
    totdenom2, firsts;
  int totdenom3, totnumer3;
  int totalcount, count1, count2, count3, count12, count13, count23, count123;
  int *interesting;
  int bestcount, thiscount;
  char *consensus;
  double mutfrac1, mutfrac2, mutfrac3, fudge, fudge2, fudge3, pvalue,
    bestpvalue;

  interesting = malloc((conLen + 1) * sizeof(int));
  consensus = malloc((conLen + 1) * sizeof(char));
  bestpvalue = 1.0;


  // find first subfamily with the label we are looking for
  for (s = 0; s < S; s++)
  {
    if (labels[s] == label)
    {
      firsts = s;
      break;
    }
  }

  // Create an array containing all positions in the consensus marked
  // with a "1" if the column contains consensus disagreement between
  // labeled subfamilies.
  for (x = 0; x < conLen; x++)
  {
    interesting[x] = 0;
    for (s = 0; s < S; s++)
    {
      if ((labels[s] == label) && (pattern[s][x] != pattern[firsts][x]))
        interesting[x] = 1;
    }
  }

  // Create a union consensus for the subfamilies with the label "label"
  for (x = 0; x < conLen; x++)
  {
    bestcount = 0;
    for (a = 0; a < 5; a++)
    {
      thiscount = 0;
      for (s = 0; s < S; s++)
      {
        if (labels[s] == label)
          thiscount += count[s][x][a];
      }
      if (thiscount > bestcount)
      {
        consensus[x] = a;
        bestcount = thiscount;
      }
    }
  }

  // Create a sorted ( by subfamily ) list of elements.
  for (x = 0; x < N; x++)
  {
    sortedAssign[x] = x;
  }
  qsort(sortedAssign, N, sizeof(sortedAssign[0]), intcmp);
  // Create a set of pointers to the start of each subfamily 
  // block within the sorted list.
  for (x = N - 1; x >= 0; x--)
  {
    sortedAssignIndices[assign[sortedAssign[x]]] = x;
  }


  /*
     Now, compute pvalues 
   */
  for (x = minDist * 2; x < conLen; x++)
  {
    if (interesting[x] == 0)
      continue;
    for (a = 0; a < 5; a++)     /* try 1-mutation pos x val a */
    {
      if (consensus[x] == a)
        continue;
      ok = 0;
      for (s = 0; s < S; s++)
      {
        if ((labels[s] == label) && (pattern[s][x] == a))
          ok = 1;
      }
      if (ok == 0)
        continue;
      if (DISALLOW_CG)
      {
        if ((consensus[x] == 1) && (a == 3) && (x < conLen - 1)
            && (consensus[x + 1] == 2))
          continue;
        if ((consensus[x] == 3) && (a == 1) && (x < conLen - 1)
            && (consensus[x + 1] == 2))
          continue;
        if ((consensus[x] == 2) && (a == 0) && (x > 0)
            && (consensus[x - 1] == 1))
          continue;
        if ((consensus[x] == 0) && (a == 2) && (x > 0)
            && (consensus[x - 1] == 1))
          continue;
      }
      for (xx = minDist; xx <= x - minDist; xx++)
      {
        for (aa = 0; aa < 5; aa++)
        {
          if (consensus[xx] == aa)
            continue;
          ok = 0;
          for (s = 0; s < S; s++)
          {
            if ((labels[s] == label) && (pattern[s][xx] == aa))
              ok = 1;
          }
          if (ok == 0)
            continue;
          if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4))
            continue;
          if (DISALLOW_CG2)
          {
            if ((consensus[xx] == 1) && (aa == 3)
                && (xx < conLen - 1) && (consensus[xx + 1] == 2))
              continue;         /* CG -> TG */
            if ((consensus[xx] == 3) && (aa == 1)
                && (xx < conLen - 1) && (consensus[xx + 1] == 2))
              continue;         /* TG -> CG */
            if ((consensus[xx] == 2) && (aa == 0) && (xx > 0)
                && (consensus[xx - 1] == 1))
              continue;         /* CG -> CA */
            if ((consensus[xx] == 0) && (aa == 2) && (xx > 0)
                && (consensus[xx - 1] == 1))
              continue;         /* CA -> CG */
          }

          totalcount = 0;
          count1 = 0;
          count2 = 0;
          count12 = 0;
          for (s = 0; s < S; s++)
          {
            if (labels[s] != label)
              continue;
            totalcount += assigncount[s];
            count1 += count[s][x][a];
            count2 += count[s][xx][aa];
            count12 += bicount[s][x][a][xx][aa];
            if (!useSiegelPvalue)
            {
              totnumer1 += numer[s][x][a];
              totdenom1 += denom[s][x][a];
              totnumer2 += numer[s][xx][aa];
              totdenom2 += denom[s][xx][aa];
            }
          }
          if (count12 == 0)
            continue;

          if (useSiegelPvalue)
          {
            pvalue =
              compute_siegel_pvalue(totalcount, count1, count2,
                                    count12, dblEpsilon) + triBonferroni;
          }
          else
          {
            if (totdenom1 == 0)
              mutfrac1 = 1.0;
            else
              mutfrac1 = ((double) totnumer1) / ((double) totdenom1);
            if (totdenom2 == 0)
              mutfrac2 = 1.0;
            else
              mutfrac2 = ((double) totnumer2) / ((double) totdenom2);
            fudge = mutfrac1 * ((double) totalcount) / ((double) count1);
            fudge2 = mutfrac2 * ((double) totalcount) / ((double) count2);
            if (fudge2 < fudge)
              fudge = fudge2;
            pvalue =
              compute_pvalue(totalcount, count1, count2, count12,
                             fudge) + triBonferroni;
          }
          //printf( "--pvalue = %lf\n", pvalue);

          if ((pvalue <= PVALUETHRESH) && (we_want_best_pvalue == 0))
          {
            //printf("Found good enough bi-pvalue = %lf\n", pvalue );
            return pvalue;
          }
          if (pvalue < bestpvalue)
            bestpvalue = pvalue;

          for (xxx = 0; xxx <= xx - minDist; xxx++)
          {
            for (aaa = 0; aaa < 5; aaa++)
            {
              if (consensus[xxx] == aaa)
                continue;
              ok = 0;
              for (s = 0; s < S; s++)
              {
                if ((labels[s] == label) && (pattern[s][xxx] == aaa))
                  ok = 1;
              }
              if (ok == 0)
                continue;
              if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4) && (aaa == 4))
                continue;
              if (DISALLOW_CG3)
              {
                // CG -> TG
                if ((consensus[xxx] == 1) && (aaa == 3)
                    && (xxx < conLen - 1) && (consensus[xxx + 1] == 2))
                  continue;
                // TG -> CG
                if ((consensus[xxx] == 3) && (aaa == 1)
                    && (xxx < conLen - 1) && (consensus[xxx + 1] == 2))
                  continue;
                // CG -> CA
                if ((consensus[xxx] == 2) && (aaa == 0)
                    && (xxx > 0) && (consensus[xxx - 1] == 1))
                  continue;
                // CA -> CG
                if ((consensus[xxx] == 0) && (aaa == 2)
                    && (xxx > 0) && (consensus[xxx - 1] == 1))
                  continue;
              }

              /*
                 OK, should compute pvalue 
               */
              totalcount = 0;
              count1 = 0;
              count2 = 0;
              count3 = 0;
              count12 = 0;
              count23 = 0;
              count13 = 0;
              totnumer1 = 0;
              totnumer2 = 0;
              totnumer3 = 0;
              totdenom1 = 0;
              totdenom2 = 0;
              totdenom3 = 0;
              count123 = 0;
              for (s = 0; s < S; s++)
              {
                if (labels[s] != label)
                  continue;
                totalcount += assigncount[s];
                count1 += count[s][x][a];
                count2 += count[s][xx][aa];
                count3 += count[s][xxx][aaa];
                count12 += bicount[s][x][a][xx][aa];
                count13 += bicount[s][x][a][xxx][aaa];
                count23 += bicount[s][xx][aa][xxx][aaa];
                count123 += getTRICount(s, x, a, xx, aa, xxx, aaa);
                if (!useSiegelPvalue)
                {
                  totnumer1 += numer[s][x][a];
                  totdenom1 += denom[s][x][a];
                  totnumer2 += numer[s][xx][aa];
                  totdenom2 += denom[s][xx][aa];
                  totnumer3 += numer[s][xxx][aaa];
                  totdenom3 += denom[s][xxx][aaa];
                }
              }
              if (count12 == 0 || count13 == 0 || count23 == 0)
                continue;
              if (useSiegelPvalue)
              {
                pvalue =
                  compute_siegel_tri_pvalue(totalcount, count1,
                                            count2, count3,
                                            count12, count13,
                                            count23, count123,
                                            dblEpsilon) + triBonferroni;
              }
              else
              {
                if (totdenom1 == 0)
                  mutfrac1 = 1.0;
                else
                  mutfrac1 = ((double) totnumer1) / ((double) totdenom1);
                if (totdenom2 == 0)
                  mutfrac2 = 1.0;
                else
                  mutfrac2 = ((double) totnumer2) / ((double) totdenom2);
                if (totdenom3 == 0)
                  mutfrac3 = 1.0;
                else
                  mutfrac3 = ((double) totnumer3) / ((double) totdenom3);
                fudge = mutfrac1 * ((double) totalcount) / ((double) count1);
                fudge2 = mutfrac2 * ((double) totalcount) / ((double) count2);
                fudge3 = mutfrac3 * ((double) totalcount) / ((double) count3);
                if (fudge2 < fudge)
                  fudge = fudge2;
                if (fudge3 < fudge)
                  fudge = fudge3;

                pvalue = compute_tri_pvalue(totalcount, count1,
                                            count2, count3,
                                            count12, count13,
                                            count23, count123,
                                            fudge) + triBonferroni;
              }

              if ((pvalue <= PVALUETHRESH) && (we_want_best_pvalue == 0))
              {
                printf("Found good enough tri-pvalue = %lf\n", pvalue);
                return pvalue;
              }
              if (pvalue < bestpvalue)
                bestpvalue = pvalue;
            }
          }
        }
      }
    }
  }

  free( interesting );
  free( consensus );

  if (we_want_best_pvalue == 0)
  {
    printf("Best pvalue found was too low ( %lf ) returning 1.0\n",
           bestpvalue);
    return 1.0;
  }
  else
  {
    return bestpvalue;
  }
}



// union_pvalue():
// A low ( <= PVALUETHRESHOLD ) indicates that these two or more subfamilies
// can stand on their own.  Higher pvalue indicates they can't.
//
//   NOTE: This routine uses a global array ( labels[] ) which holds
//         one label per subfamily.  It uses this along with the scalar "label"
//         to determine which subfamilies ( 2 or more ) will be used in the
//         calculation. 
//
double
union_pvalue(int label)
{
  int s, x, a, xx, aa, ok, totnumer1, totdenom1, totnumer2, totdenom2, firsts;
  int totalcount, count1, count2, count12;
  int *interesting;
  int bestcount, thiscount;
  char *consensus;
  double mutfrac1, mutfrac2, fudge, fudge2, pvalue, bestpvalue;

  interesting = malloc((conLen + 1) * sizeof(int));
  consensus = malloc((conLen + 1) * sizeof(char));
  bestpvalue = 1.0;

  for (s = 0; s < S; s++)
  {
    if (labels[s] == label)
    {
      firsts = s;
      break;
    }
  }

  // Label positions which differ between two ( or more ) subfamily consensi
  for (x = 0; x < conLen; x++)
  {
    interesting[x] = 0;
    for (s = 0; s < S; s++)
    {
      if ((labels[s] == label) && (pattern[s][x] != pattern[firsts][x]))
        interesting[x] = 1;
    }
  }

  // Build a consensus based on the elements in two ( or more ) subfamilies
  for (x = 0; x < conLen; x++)
  {
    bestcount = 0;
    for (a = 0; a < 5; a++)
    {
      thiscount = 0;
      for (s = 0; s < S; s++)
      {
        if (labels[s] == label)
          thiscount += count[s][x][a];
      }
      if (thiscount > bestcount)
      {
        consensus[x] = a;
        bestcount = thiscount;
      }
    }
  }

  /*
     Now, compute pvalues 
   */
  for (x = minDist; x < conLen; x++)
  {
    // Do not consider this site for mutation if disabled.
    if (useDisabledSites && disabledSites[x])
      continue;
    if (interesting[x] == 0)
      continue;
    for (a = 0; a < 5; a++)     /* try 1-mutation pos x val a */
    {
      if (consensus[x] == a)
        continue;
      ok = 0;
      // Make sure at least one subfamily has "a" at "x" as a consensus.
      for (s = 0; s < S; s++)
      {
        if ((labels[s] == label) && (pattern[s][x] == a))
          ok = 1;
      }
      if (ok == 0)
        continue;
      if (DISALLOW_CG)
      {
        if ((consensus[x] == 1) && (a == 3) && (x < conLen - 1)
            && (consensus[x + 1] == 2))
          continue;
        if ((consensus[x] == 3) && (a == 1) && (x < conLen - 1)
            && (consensus[x + 1] == 2))
          continue;
        if ((consensus[x] == 2) && (a == 0) && (x > 0)
            && (consensus[x - 1] == 1))
          continue;
        if ((consensus[x] == 0) && (a == 2) && (x > 0)
            && (consensus[x - 1] == 1))
          continue;
      }
      if (DEBUG >= 3)
        printf("***UNION pos %d, %c -> %c\n", x,
               num_to_char(consensus[x]), num_to_char(a));
      for (xx = 0; xx <= x - minDist; xx++)
      {
        // Do not consider this site for mutation if disabled.
        if (useDisabledSites && disabledSites[xx])
          continue;
        for (aa = 0; aa < 5; aa++)
        {
          if (consensus[xx] == aa)
            continue;
          ok = 0;
          for (s = 0; s < S; s++)
          {
            if ((labels[s] == label) && (pattern[s][xx] == aa))
              ok = 1;
          }
          if (ok == 0)
            continue;
          if (DEBUG >= 3)
            printf("  ++UNION pos %d, %c -> %c\n", xx,
                   num_to_char(consensus[xx]), num_to_char(aa));
          if ((DISALLOW_DASHDASH) && (a == 4) && (aa == 4))
            continue;
          if (DISALLOW_CG2)
          {
            if ((consensus[xx] == 1) && (aa == 3)
                && (xx < conLen - 1) && (consensus[xx + 1] == 2))
              continue;         /* CG -> TG */
            if ((consensus[xx] == 3) && (aa == 1)
                && (xx < conLen - 1) && (consensus[xx + 1] == 2))
              continue;         /* TG -> CG */
            if ((consensus[xx] == 2) && (aa == 0) && (xx > 0)
                && (consensus[xx - 1] == 1))
              continue;         /* CG -> CA */
            if ((consensus[xx] == 0) && (aa == 2) && (xx > 0)
                && (consensus[xx - 1] == 1))
              continue;         /* CA -> CG */
          }
          /*
             OK, should compute pvalue 
           */
          totalcount = 0;
          count1 = 0;
          count2 = 0;
          count12 = 0;
          for (s = 0; s < S; s++)
          {
            if (labels[s] != label)
              continue;
            totalcount += assigncount[s];
            count1 += count[s][x][a];
            count2 += count[s][xx][aa];
            count12 += bicount[s][x][a][xx][aa];
            if (!useSiegelPvalue)
            {
              totnumer1 += numer[s][x][a];
              totdenom1 += denom[s][x][a];
              totnumer2 += numer[s][xx][aa];
              totdenom2 += denom[s][xx][aa];
            }
          }
          // used to be == 0
          if (count12 < minCount)
            continue;
          if (useSiegelPvalue)
          {
            pvalue =
              compute_siegel_pvalue(totalcount, count1, count2,
                                    count12, dblEpsilon) + dblBonferroni;
          }
          else
          {
            if (totdenom1 == 0)
              mutfrac1 = 1.0;
            else
              mutfrac1 = ((double) totnumer1) / ((double) totdenom1);
            if (totdenom2 == 0)
              mutfrac2 = 1.0;
            else
              mutfrac2 = ((double) totnumer2) / ((double) totdenom2);
            fudge = mutfrac1 * ((double) totalcount) / ((double) count1);
            fudge2 = mutfrac2 * ((double) totalcount) / ((double) count2);
            if (fudge2 < fudge)
              fudge = fudge2;
            pvalue =
              compute_pvalue(totalcount, count1, count2, count12,
                             fudge) + dblBonferroni;
          }
          //printf( "++pvalue = %lf\n", pvalue);
          if (DEBUG >= 3)
            printf("UNION pos %d, %c -> %c  & pos %d, %c -> %c = "
                   "pvalue = %lf\n", x, num_to_char(consensus[x]),
                   num_to_char(a), xx, num_to_char(consensus[xx]),
                   num_to_char(aa), pvalue);
          if ((pvalue <= PVALUETHRESH) && (we_want_best_pvalue == 0))
            return pvalue;
          if (pvalue < bestpvalue)
            bestpvalue = pvalue;
        }
      }
    }
  }

  free( interesting );
  free( consensus );

  if (we_want_best_pvalue == 0)
    return 1.0;
  else
    return bestpvalue;
}



void
merge_subfamilies(int label)
{
  int firsts, x, a, xx, aa, s, n;

  for (s = 0; s < S; s++)
  {
    if (labels[s] == label)
    {
      firsts = s;
      break;
    }
  }

  for (s = S - 1; s > firsts; s--)
  {
    if (labels[s] != label)
      continue;
    /*
       Now, merge s into firsts 
     */
    for (n = 0; n < N; n++)
    {
      if (assign[n] != s)
        continue;
      assign[n] = firsts;
      assigncount[firsts]++;
      assigncount[s]--;
      for (x = 0; x < conLen; x++)
      {
        count[firsts][x][ele[n][x]]++;
        count[s][x][ele[n][x]]--;
        // RMH added 11/8/12
        counti[firsts][x][elei[n][x]]++;
        counti[s][x][elei[n][x]]--;
      }
      for (x = minDist; x < conLen; x++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
        {
          bicount[firsts][x][ele[n][x]][xx][ele[n][xx]]++;
          bicount[s][x][ele[n][x]][xx][ele[n][xx]]--;
        }
      }
    }
    /*
       Now reassign S-1 to s 
     */
    for (n = 0; n < N; n++)
    {
      if (assign[n] == S - 1)
        assign[n] = s;
    }
    assigncount[s] = assigncount[S - 1];
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
        count[s][x][a] = count[S - 1][x][a];
      // RMH addd 11/8/12
      for (a = 0; a < 2; a++)
        counti[s][x][a] = counti[S - 1][x][a];
    }
    for (x = minDist; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
      {
        for (xx = 0; xx <= x - minDist; xx++)
        {
          for (aa = 0; aa < 5; aa++)
            bicount[s][x][a][xx][aa] = bicount[S - 1][x][a][xx][aa];
        }
      }
    }
    S--;
  }
  assign_to_pattern();
}



int
compute_distance(int s, int t, int insPenalty)
{
  int x, a, b, dist;

  dist = 0;

  // OLD Code Difference:
  // Alkes originally used a greater distance for indels than for mismatches.
  // The strange thing is that the penalty was really high ( 100 ) for 
  // each deletion ( ie "--" = 202 ) and a constant ( also 101 ) for any
  // length insertion. It only seems to be used to calculate the distance 
  // for the the final pvalue scaffold calculations and not in the building
  // process.
  //
  // Scaffold build process: mismatch/complete_insertion/deletion=-1
  //   for(x=0; x<L; x++)
  //   {
  //     if(pattern[s][x] != pattern[t][x]) dist++;
  //     if(patterni[s][x] != patterni[t][x]) dist++;
  //   }
  //
  // Fill Scaffold and Compute_logpvalues processes: mismatch = -1,
  //                                                 complete_insertion = 101, 
  //                                                 deletion(each) = 101,
  //  for(x=0; x<conLen; x++)
  //  {
  //    if(pattern[s][x] != pattern[t][x]) dist++;
  //    if(patterni[s][x] != patterni[t][x]) dist++;
  //
  //    /* ignore freq indels in middle A-rich region */
  //    if((x > 115) && (x <= 135)) continue; 
  //    if(patterni[s][x] != patterni[t][x]) dist += 100;
  //    if(pattern[s][x] == pattern[t][x]) continue;
  //    if((pattern[s][x]==4) || (pattern[t][x]==4)) dist += 100;
  //  }
  // 
  //
  // NOW: 
  //       Mismatch = -1,
  //       Insertion ( any length ) = insPenalty ( 1 for building and
  //                                               bi/tri and 11 for
  //                                               building singlemutation 
  //                                               and tree building )
  //       Deletion = -1
  for (x = 0; x < conLen; x++)
  {
    // Do not consider this site for mutation if disabled.
    if (useDisabledSites && disabledSites[x])
      continue;

    // Default mismatch
    if (pattern[s][x] != pattern[t][x])
      dist++;

    // Default insertion penalty ( 1  or 11 currently )
    if (patterni[s][x] != patterni[t][x])
      dist += insPenalty;

    // Subtract off CpG mutations from dist 
    if (DISALLOW_CG)
    {
      if (x < conLen - 1)
      {
        // SEQ1:CG->SEQ2:T* or SEQ1:CG->SEQ2:*A
        if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 2))
        {
          if (pattern[t][x] == 3)
            dist--;
          else if (pattern[t][x + 1] == 0)
            dist--;
        }
        // SEQ2:CG->SEQ1:T* or SEQ2:CG->SEQ1:*A
        else if ((pattern[t][x] == 1) && (pattern[t][x + 1] == 2))
        {
          if (pattern[s][x] == 3)
            dist--;
          else if (pattern[s][x + 1] == 0)
            dist--;
        }
        // SEQ1:TG->SEQ2:CA
        else if ((pattern[s][x] == 3) && (pattern[s][x + 1] == 2))
        {
          if ((pattern[t][x] == 1) && (pattern[t][x + 1] == 0))
            dist--;
        }
        // SEQ1:CA->SEQ2:TG
        else if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 0))
        {
          if ((pattern[t][x] == 3) && (pattern[t][x + 1] == 2))
            dist--;
        }
      }
    }
  }

  return dist;
}


// 
// build_MST_scaffold() - originally "build_MST"
//
// Build the Minimum Spanning Tree ( MST ) of the scaffold subfamilies
// Seed the building processes with the oldest subfamily ( cluster with 
// the most average divergence of it's members from it's consensus ) and
// build a MST using consensus hamming (almost) distances. 
//
//  Input:  distance[][] ( set by prune_subfamilies ) 
//
//  Reuse:  Clears and re-defines parent[]
//
//  Output: parent[]
//          mstLogPValues[]
//          Optionaly saves MST to a file.
//
void
build_MST_scaffold(char *filename)
{
  int s, t, s0, done[S], numdone, thiss, thist, mindist, x;
  int root;
  double sage, s0age, pvalue;
  FILE *fp = NULL;

  if (VERBOSE)
    fprintf(outFP, "MINIMUM SPANNING TREE:\n");
  if (VERBOSE)
    printf("MINIMUM SPANNING TREE:\n");

  // Only save tree data if we need to
  if (filename)
  {
    // Create the viz file 
    if ((fp = fopen(filename, "w")) == NULL)
    {
      printf("Could not open input file %s\n", filename);
      exit(1);
    }
    fprintf(fp, "digraph G {\n");
    fprintf(fp, "        size=\"8,10\";\n");
  }

  // Choose the oldest subfamily s0, with age s0age 
  //   Age is determined by the sum of the consensus position percentage
  //   mismatches.  Ie. the higher the "sage" the higher the overall divergence
  //   in the cluster.
  for (s = 0; s < S; s++)
    done[s] = 0;

  s0 = 0;
  s0age = 0.0;
  for (s = 0; s < S; s++)
  {
    sage = 0.0;
    for (x = 0; x < conLen; x++)
    {
      sage += 1.0 -
        ((double) count[s][x][pattern[s][x]]) / ((double) assigncount[s]);
    }
    if (sage > s0age)
    {
      s0age = sage;
      s0 = s;
    }
  }
  done[s0] = 1;

  // RMH:
  // parent[] array is being reused here.  It was previously used to hold 
  //          the source subfamily for each derived subfamily. Here it 
  //          appears to be used to hold the minimum spanning tree structure.
  //parent[s0] = -1;
  //  RMH: Why was this full-clear disabed in my version?
  for(s=0; s<S; s++) parent[s] = -1;

  root = s0;
  mstLogPValues[root] = 0;

  if (VERBOSE)
    print_subfamily(s0);

  numdone = 1;
  while (numdone < S)
  {
    thiss = -1;
    thist = -1;
    mindist = 1000000000.0;
    for (s = 0; s < S; s++)     
    {
      // Only want existing nodes in the tree
      if (done[s] == 0)
        continue;
      for (t = 0; t < S; t++)  
      {
        // Only use unassigned nodes
        if (done[t] == 1)
          continue;             /* guarantees t != s */
        if (distance[s][t] < mindist)
        {
          mindist = distance[s][t];
          thiss = s;
          thist = t;
        }
      }
    }
    if (thiss < 0)
    {
      printf("OOPS numdone=%d thiss<0\n", numdone);
      exit(1);
    }
    s = thiss;
    t = thist;

    done[t] = 1;
    numdone++;
    parent[t] = s;

    // Setup array used by untion_tri_pvalue() and union_pvalue()
    for (s0 = 0; s0 < S; s0++)
      labels[s0] = 0;
    labels[s] = 1;
    labels[t] = 1;

    // So...this pvalue represents the most significant pvalue in 
    // the union of subfamily t with one of it's closest (distance)
    // subfamilies.
    if (useTRI)
    {
      pvalue = union_tri_pvalue(1);
    }
    else
    {
      pvalue = union_pvalue(1);
    }
    mstLogPValues[t] = pvalue;
    //printf( "Calculating union pvalue s=%d t=%d: %lf\n", s, t, pvalue );
    if ((s == root) && (pvalue < mstLogPValues[root]))
      mstLogPValues[root] = pvalue;
    if (fp)
      fprintf(fp, "        %d -> %d [label = \"p=%e,c=%d\"];\n", s, t, 
              exp(pvalue), assigncount[t]);
    if (VERBOSE)
      print_subfamily(t);
  }

  // Only save MST data if requested.
  if (fp)
  {
    fprintf(fp, "}\n");
    fclose(fp);
  }

  // RMH: Sanity Check
  int rootCnt = 0;
  for(s=0; s<S; s++) 
  {
    if ( parent[s] == -1 )
    {
      printf("Tree root = subfamily-%d\n", s );
      rootCnt++;
    }
  }
  if ( rootCnt > 1 )
  {
    printf("Error: tree has more than one root!!!\n");
    exit(1);
  }

}


void
print_subfamily(int s)
{
  int x;
  double sage, sage2, sagedenom, sagedenom2;

  sage = 0.0;
  sage2 = 0.0;
  sagedenom = 0.0;
  sagedenom2 = 0.0;
  for (x = 0; x < conLen; x++)
  {
    sage += 1.0 - ((double) count[s][x][pattern[s][x]]) /
      ((double) assigncount[s]);
    sagedenom += 1.0;
    if ((x < conLen - 1) && (pattern[s][x] == 1) && (pattern[s][x + 1] == 2))
    {
      x++;
      continue;
    }
    if ((x < conLen - 1) && (pattern[s][x] == 1) && (pattern[s][x + 1] == 0))
    {
      x++;
      continue;
    }
    if ((x < conLen - 1) && (pattern[s][x] == 3) && (pattern[s][x + 1] == 2))
    {
      x++;
      continue;
    }
    sage2 += 1.0 - ((double) count[s][x][pattern[s][x]]) /
      ((double) assigncount[s]);
    sagedenom2 += 1.0;
  }
  sage /= sagedenom;
  sage2 /= sagedenom2;

  fprintf(outFP, "Subfamily %d (parent %d count %d age %.03f age2 %.03f):\n",
          s, parent[s], assigncount[s], sage, sage2);
  for (x = 0; x < conLen; x++)
  {
    if (pattern[s][x] != Sxsequence[x])
      fprintf(outFP, "%d:%c ", x, num_to_char(pattern[s][x]));
    if (patterni[s][x] == 1)
      fprintf(outFP, "%d:+ ", x);
  }
  fprintf(outFP, "\n");
  fprintf(outFP, ">Subfamily%d (parent %d count %d age %.03f age2 %.03f):\n",
          s, parent[s], assigncount[s], sage, sage2);
  for (x = 0; x < conLen; x++)
  {
    fprintf(outFP, "%c", x, num_to_char(pattern[s][x]));
  }
  fprintf(outFP, "\n");

  printf("Subfamily %d (parent %d count %d age %.03f age2 %.03f):\n", s,
         parent[s], assigncount[s], sage, sage2);
  for (x = 0; x < conLen; x++)
  {
    if (pattern[s][x] != Sxsequence[x])
      printf("%d:%c ", x, num_to_char(pattern[s][x]));
    if (patterni[s][x] == 1)
      printf("%d:+ ", x);
  }
  printf("\n");
}

/*************************************************************
 *************  Single Mutation Subfamily Code ***************
 *************************************************************/

// From singlemut
double
compute_sigma(int totalcount, int count1, double emutfrac)
{
  double mean, sdev;

  mean = emutfrac * ((double) totalcount);
  if (((double) count1) <= mean)
    return 0.0;
  if ((mean <= 0) || (1.0 - emutfrac <= 0))
    return 0.0;
  sdev = exp(0.5 * log(mean * (1.0 - emutfrac)));

  return (((double) count1) - mean) / sdev;
}

// From singlemut
void
compute_bestmut1()
{
  int SS, x, a, s, x0, thisdist, nogood;
  double emutfrac, sigma;

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      existingval[x][a] = 0;
  }
  for (s = 0; s < S; s++)
  {
    for (x = 0; x < conLen; x++)
      existingval[x][pattern[s][x]] = 1;
  }

  bestmutsigma = sigmaThreshold;
  for (SS = 0; SS < S; SS++)    /* try splitting subfamily SS */
  {
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)   /* try 1-mutation pos x val a */
      {
        if ((NEW_MUT_ONLY) && (existingval[x][a] == 1))
          continue;
        if (DISALLOW_INDEL)
        {
          if ((pattern[SS][x] == 4) || (a == 4))
            continue;
        }
        if (pattern[SS][x] == a)
          continue;
        if (count[SS][x][a] < minCount)
          continue;
        /*
           new code to disallow CG->TG etc 
         */
        if (DISALLOW_CG)
        {
          if ((pattern[SS][x] == 1) && (a == 3) &&
              (x < conLen - 1) && (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 3) && (a == 1) &&
              (x < conLen - 1) && (pattern[SS][x + 1] == 2))
            continue;
          if ((pattern[SS][x] == 2) && (a == 0) &&
              (x > 0) && (pattern[SS][x - 1] == 1))
            continue;
          if ((pattern[SS][x] == 0) && (a == 2) &&
              (x > 0) && (pattern[SS][x - 1] == 1))
            continue;
          if ((pattern[SS][x] == 1) && (a == 3) &&
              (x < conLen - 1) && (pattern[SS][x + 1] == 0))
            continue;
          if ((pattern[SS][x] == 2) && (a == 0) &&
              (x > 0) && (pattern[SS][x - 1] == 3))
            continue;
        }

        emutfrac = age[SS] * mutperage[pattern[SS][x]][x][a];
        if ((emutfrac < 0.0) || (emutfrac > 1.0))
        {
          printf("OOPS emutfrac=%f\n", emutfrac);
          exit(1);
        }
        sigma = compute_sigma(assigncount[SS], count[SS][x][a], emutfrac);
        // DEBUG
        //{ printf("SS = %d Pos=%d %c -> %c %lf ( have to beat %lf )\n", SS, x, num_to_char( pattern[SS][x] ), num_to_char( a ), sigma, bestmutsigma ); }

        if (sigma <= bestmutsigma)
          continue;

        /*
           One last check to make sure result is not distance 0 from any subf 
         */
        for (x0 = 0; x0 < conLen; x0++)
        {
          pattern[S][x0] = pattern[SS][x0];
          patterni[S][x0] = patterni[SS][x0];
        }
        pattern[S][x] = a;
        nogood = 0;
        for (s = 0; s < S; s++)
        {
          if (compute_distance(S, s, SINGLE_MUTATION_INSPENALTY) == 0)
            nogood = 1;
        }
        if (nogood == 1)
          continue;

        /*
           WE HAVE A NEW WINNER! 
         */
        bestmutsigma = sigma;
        bestmutSS = SS;
        bestmutx = x;
        bestmuta = a;
      }
    }
  }
}

// from singlemut obviously
// TODO: merge with old routine?
void
assign_to_pattern_singlemut()
{
  int n, s, x, a, b, maxcount, c, d;

  for (s = 0; s < S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
        count[s][x][a] = 0;
      for (b = 0; b < 2; b++)
        counti[s][x][b] = 0;
    }
  }

  for (n = 0; n < N; n++)
  {
    s = assign[n];
    if (s < 0)
      continue;
    for (x = 0; x < conLen; x++)
    {
      a = ele[n][x];
      b = elei[n][x];
      count[s][x][a]++;
      counti[s][x][b]++;
    }
  }

  for (s = 0; s < S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      maxcount = -1;
      for (a = 0; a < 5; a++)
      {
        if (count[s][x][a] > maxcount)
        {
          maxcount = count[s][x][a];
          pattern[s][x] = a;
        }
      }
      maxcount = -1;
      for (b = 0; b < 2; b++)
      {
        if (counti[s][x][b] > maxcount)
        {
          maxcount = counti[s][x][b];
          patterni[s][x] = b;
        }
      }
    }
  }

  /*
     next, compute profiles 
   */
  //printf("assign_to_pattern_singlemut() - compute profiles\n");
  for (s = 0; s < S; s++)
  {
    if (assigncount[s] == 0)
    {
      printf("OOPS s=%d assigncount=0\n", s);
      exit(1);
    }
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
      {
        profile[s][x][a] =
          ((double) count[s][x][a]) / ((double) assigncount[s]);
      }
      for (b = 0; b < 2; b++)
      {
        profilei[s][x][b] =
          ((double) counti[s][x][b]) / ((double) assigncount[s]);
      }
    }
  }

  /*
     Next, compute age 
   */
  //printf
  //  ("assign_to_pattern_singlemut() - compute age MAXS = %d, S = %d, conLen = %d\n",
  //   MAXS, S, conLen);
  for (s = 0; s < MAXS; s++)
    age[s] = 0.0;

  for (s = 0; s < S; s++)
  {
    for (x = 0; x < conLen; x++)
    {
      // By this calculation age[s] will be a number between 0-200.  
      //   - I suspect Alkes assumed that for Alu it probably wouldn't
      //     exceed 100.  
      age[s] += ((double) 1.0) - profile[s][x][pattern[s][x]];
      age[s] += ((double) 1.0) - profilei[s][x][patterni[s][x]];
    }
    // New: Scale result back to 0-100
    age[s] = (age[s] / 200) * 100;
    //printf("assign_pat_to_sing_mut: age[%d] = %f , patterni[43][0] = %d\n", s,
    //       age[s], patterni[43][0]);
  }

  /*
     Next, compute mutperage 
   */
  for (x = 0; x < conLen; x++)
  {
    for (c = 0; c < 5; c++)
    {
      perage[c][x] = 0.0;
      for (a = 0; a < 5; a++)
        mut[c][x][a] = 0.0;
    }
    for (d = 0; d < 2; d++)
    {
      peragei[d][x] = 0.0;
      for (b = 0; b < 2; b++)
        muti[d][x][b] = 0.0;
    }
  }
  for (n = 0; n < N; n++)
  {
    s = assign[n];
    for (x = 0; x < conLen; x++)
    {
      a = ele[n][x];
      b = elei[n][x];
      c = pattern[s][x];
      if (a != c)
        mut[c][x][a] += 1.0;
      perage[c][x] += age[s];

      d = patterni[s][x];
      if (b != d)
        muti[d][x][b] += 1.0;
      peragei[d][x] += age[s];
    }
  }
  for (x = 0; x < conLen; x++)
  {
    for (c = 0; c < 5; c++)
    {
      if (perage[c][x] > 0.0)
      {
        for (a = 0; a < 5; a++)
        {
          if (a == c)
            continue;
          mutperage[c][x][a] = mut[c][x][a] / perage[c][x];
          //if ( x == 411 )
          //  printf("mutperage[%d][411][%d] = mut[c][x][a]=%f perage[c][x]=%f\n", c, a, mut[c][x][a], perage[c][x] );
        }
      }
      else
      {
        for (a = 0; a < 5; a++)
        {
          if (a == c)
            continue;
          mutperage[c][x][a] = 1.0 / (((double) conLen) * 5.0);
          //if ( x == 411 )
          //  printf("mutperage[%d][411][%d] = 1/conLen*5\n", c, a,mut[c][x][a], perage[c][x] );
        }
      }
    }
    for (d = 0; d < 2; d++)
    {
      if (peragei[d][x] > 0.0)
      {
        for (b = 0; b < 2; b++)
        {
          if (b == d)
            continue;
          mutperagei[d][x][b] = muti[d][x][b] / peragei[d][x];
        }
      }
      else
      {
        for (b = 0; b < 2; b++)
        {
          if (b == d)
            continue;
          mutperagei[d][x][b] = 1.0 / (((double) conLen) * 5.0);
        }
      }
    }
  }
}


// for Singlemut
// Same as assign_to_pattern except that it appears that
// specific subfamilies can be left out.
void
assign_to_pattern_mark()        /* M-step */
{
  int n, s, x, a, b, maxcount, c, d;

  for (s = 0; s < S; s++)
  {
    if (mark[s] == 0)
      continue;
    for (x = 0; x < conLen; x++)
    {
      for (a = 0; a < 5; a++)
        count[s][x][a] = 0;
      for (b = 0; b < 2; b++)
        counti[s][x][b] = 0;
    }
  }

  for (n = 0; n < N; n++)
  {
    s = assign[n];
    if (mark[s] == 0)
      continue;
    for (x = 0; x < conLen; x++)
    {
      a = ele[n][x];
      b = elei[n][x];
      count[s][x][a]++;
      counti[s][x][b]++;
    }
  }

  for (s = 0; s < S; s++)
  {
    if (mark[s] == 0)
      continue;
    for (x = 0; x < conLen; x++)
    {
      maxcount = -1;
      for (a = 0; a < 5; a++)
      {
        if (count[s][x][a] > maxcount)
        {
          maxcount = count[s][x][a];
          pattern[s][x] = a;
        }
      }
      maxcount = -1;
      for (b = 0; b < 2; b++)
      {
        if (counti[s][x][b] > maxcount)
        {
          maxcount = counti[s][x][b];
          patterni[s][x] = b;
        }
      }
    }
  }
}


// for Singlemut
void
pattern_to_assign_mark(int thisx)       /* E-step */
{
  int n, s, x, a, xx, aa, b, mindist, dist, ok, oldassign, newassign;

  //printf("pattern_to_assign_mark\n");
  for (n = 0; n < N; n++)
  {
    oldassign = assign[n];
    if (mark[oldassign] == 0)
      continue;
    ok = 0;
    mindist = 1000000000;
    for (s = 0; s < S; s++)
    {
      if (mark[s] == 0)
        continue;
      dist = 0;
      for (x = 0; x < conLen; x++)
      {
        a = ele[n][x];
        if (a != pattern[s][x])
          dist++;
        b = elei[n][x];
        if (b != patterni[s][x])
          dist++;
      }
      //printf("thisx = %d n=%d s=%d mutdist=%d",thisx,n,s,dist);
      /*
         Subtract off CpG mutations from dist 
       */
      if (DISALLOW_CG)
      {
        for (x = 0; x < conLen - 1; x++)
        {
          // RMH: This seems to be modified by me to 
          //      ignore CpG mutations in the diagnostic position
          //      when computing the distance.
          if (x == thisx)
            continue;
          if (x == thisx - 1)
            continue;
          if (x == thisx + 1)
            continue;
          if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 2)) /* CpG */
          {
            if (ele[n][x] == 3)
            {
              dist--;
            }
            else if (ele[n][x + 1] == 0)
            {
              dist--;
            }
          }
          else if ((pattern[s][x] == 3) && (pattern[s][x + 1] == 2))    /* TpG */
          {
            if (ele[n][x] == 1)
            {
              dist--;
            }
          }
          else if ((pattern[s][x] == 1) && (pattern[s][x + 1] == 0))    /* CpA */
          {
            if (ele[n][x + 1] == 2)
            {
              dist--;
            }
          }
        }
      }
      //printf(" newdist=%d element=%c%c%c pattern=%c%c%c\n",
      //       dist, num_to_char( (char)ele[n][thisx-1] ),
      //       num_to_char( (char)ele[n][thisx] ),
      //       num_to_char( (char)ele[n][thisx+1] ),
      //       num_to_char( (char)pattern[s][thisx-1] ),
      //       num_to_char( (char)pattern[s][thisx] ),
      //       num_to_char( (char)pattern[s][thisx+1] ) );
      if (dist < mindist)
      {
        mindist = dist;
        assign[n] = s;
        ok = 1;
      }
      else if (dist == mindist)
        ok = 0;
    }

    newassign = assign[n];
    if (newassign != oldassign)
    {
      //printf("Shifting from sub %d (%d remain) to %d (%d remain)\n", oldassign, (assigncount[oldassign] - 1), newassign, (assigncount[newassign] + 1 ) );
      assigncount[oldassign]--;
      for (x = 0; x < conLen; x++)
        count[oldassign][x][ele[n][x]]--;

      assigncount[newassign]++;
      for (x = 0; x < conLen; x++)
        count[newassign][x][ele[n][x]]++;
    }
  }
}

// from singlemut
// TODO: Possibly merge with build_new_subfamily
void
build_new_subfamily2()
{
  int n, s, thisx, thisa, x, a, SS, besta, bestb, b;
  int ccount[5], ccounti[2];
  int oldcount, newcount;
  int totalcount, count1;
  double emutfrac;
  double bestscore, countratio;

  /*
     we are given bestmutSS, bestmutx, bestmuta 
   */
  SS = bestmutSS;
  thisx = bestmutx;
  thisa = bestmuta;
  parent[S] = SS;

  totalcount = assigncount[SS];
  count1 = count[SS][thisx][thisa];
  emutfrac = age[SS] * mutperage[pattern[SS][thisx]][thisx][thisa];
  fprintf(outFP,
          "Building subfamily %d (parent %d, logpvalue %f, sigma %.2f): pos %d %c to %c (%d %d %f)\n",
          S, SS, sigmage_to_logpvalue(bestmutsigma), bestmutsigma, thisx,
          num_to_char(pattern[SS][thisx]), num_to_char(thisa), totalcount,
          count1, emutfrac);
  printf
    ("  Building subfamily %d (parent %d, logpvalue %f, sigma %.2f): pos %d %c to %c (%d %d %f)\n",
     S, SS, sigmage_to_logpvalue(bestmutsigma), bestmutsigma, thisx,
     num_to_char(pattern[SS][thisx]), num_to_char(thisa), totalcount,
     count1, emutfrac);

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      count[S][x][a] = 0;
  }
  for (n = 0; n < N; n++)
  {

    if ((assign[n] == SS) && (ele[n][thisx] == thisa))
    {
      for (x = 0; x < conLen; x++)
      {
        count[S][x][ele[n][x]]++;
        count[SS][x][ele[n][x]]--;
      }
      assign[n] = S;
      assigncount[S]++;
      assigncount[SS]--;
    }
  }
  for (x = 0; x < conLen; x++)
  {
    pattern[S][x] = pattern[SS][x];
    patterni[S][x] = patterni[SS][x];
  }
  pattern[S][thisx] = thisa;
  S++;

  for (s = 0; s < S; s++)
    mark[s] = 0;
  mark[SS] = 1;
  mark[S - 1] = 1;
  assign_to_pattern_mark();
}

static char *
pValueStr(double logpvalue)
{
  double aa, c, d;
  int b, x, y;
  static char pvStr[255];

  aa = fabs(logpvalue * (double) 0.434294481);
  b = (int) aa;
  c = aa - b;
  d = exp(-c / (double) 0.434294481);
  if (d >= 0.95)
  {
    x = 1;
    y = -b;
  }
  else
  {
    x = (int) ((double) 10 * d + (double) 0.5);
    y = -(b + 1);
  }
  sprintf(pvStr, "%de%d", x, y);
  return pvStr;
}


//
// build_MST_full() - originaly build_MST2 from fill_scaffold.c
//
//   Build a minimum spanning tree (MST) using scaffold +
//   single point mutation derived subfamilies.  
//
//   In this tree build the insertion penalty is changes
//   to SINGLE_MUTATION_INSPENALTY ( currently ) which differs
//   from the scaffold only MST build which uses 
//   MULTIPLE_MUTATION_INSPENALTY
//
//   
//
//   Recomputes distance[][] using different insertion pentalty
//
//  Uses:
//      edges[][]    : overwrites previous values
//      distance[][] : overwrites previous values
//                        NOTE: It is currently hardcoded to use
//                              SINGLE_MUTATION_INSPENALTY.  If this
//                              is dropped back into service *before*
//                              single mutations are calculated reconsider
//                              making this a parameter.
//      parent[]     : overwrites previous values and redefines
//                     the meaning of "parent".
//
//  OUTPUT
//      parent[] 
//      mstLogPValues[]
//
void
build_MST_full(char *filename)
{
  int s, t, s0, done[S], numdone, thiss, thist, mindist, x;
  int isleaf[S], nleaves;
  double sage, s0age, pvalue;
  FILE *fp;
  double mutrate[MAXS], mutrate_CpGMod[MAXS];
  double diameter, hue;

  // Get mutation rates
  computeMutationRatesKimura(mutrate, mutrate_CpGMod);
  double minMutRate = 1.0;
  double maxMutRate = 0.0;
  int maxCount = 0;
  double mutRateRange = 0.0;
  for (s = 0; s < S; s++)
  {
    if (mutrate_CpGMod[s] > maxMutRate)
      maxMutRate = mutrate_CpGMod[s];
    if (mutrate_CpGMod[s] < minMutRate)
      minMutRate = mutrate_CpGMod[s];
    if (assigncount[s] > maxCount)
      maxCount = assigncount[s];
  }
  mutRateRange = maxMutRate - minMutRate;
  if (mutRateRange == 0)
    mutRateRange = 1.0;

  for (s = 0; s < S; s++)
    isleaf[s] = 1;

  for (s = 0; s < S; s++)
  {
    for (t = 0; t < S; t++)
      edges[s][t] = 0;
  }

  for (s = 0; s < S; s++)
  {
    for (t = s + 1; t < S; t++)
    {
      distance[s][t] = compute_distance(s, t, SINGLE_MUTATION_INSPENALTY);
      distance[t][s] = distance[s][t];
    }
  }

  if ( filename )
  {
    // create the viz file
    if ((fp = fopen(filename, "w")) == NULL)
    {
      printf("Could not open input file %s\n", filename);
      exit(1);
    }
    fprintf(fp, "digraph G {\n");
    fprintf(fp, "        size=\"8,10\";\n");
  }

  /*
     choose the oldest subfamily s0 
   */
  for (s = 0; s < S; s++)
    done[s] = 0;
  s0 = 0;
  s0age = 0.0;
  for (s = 0; s < S; s++)
  {
    sage = 0.0;
    for (x = 0; x < conLen; x++)
    {
      sage += 1.0 -
        ((double) count[s][x][pattern[s][x]]) / ((double) assigncount[s]);
    }
    if (sage > s0age)
    {
      s0age = sage;
      s0 = s;
    }
  }
  done[s0] = 1;
// RMH:
// was   parent[s0] = -1;
  for (s = 0; s < S; s++)
    parent[s] = -1;

  // Optional output data to graph file.
  if ( fp ) 
  {
    // Circle area is proportional to the size of the subfamily 
    diameter = (double) 2.0 *sqrt((double) assigncount[s0] / (double) maxCount);
  
    // Color heatmap ranges from 0-0.6 ( or 0-216 degrees )
    if (useOriginalGraphColors)
      hue = ((mutrate_CpGMod[s0] - minMutRate) * (double) 0.6) / mutRateRange;
    else
      hue =
        ((double) 1.0 -
         ((mutrate_CpGMod[s0] - minMutRate) / mutRateRange)) * (double) 0.6;

    if (s0 <= lastScaffoldIndex)
      fprintf(fp,
              "        %d [shape=circle, label=\"sub%d\\nc=%d\\npv=%s\\ndiv=%.3f\", style=filled,"
              " height=%.1f, width=%.1f, color=\"%.2f,0.9,0.6\"", s0, s0,
              assigncount[s0], pValueStr(mstLogPValues[s0]), mutrate_CpGMod[s0],
              diameter, diameter, hue);
    else
      fprintf(fp,
              "        %d [shape=circle, label=\"sub%d\\nc=%d\\npv=%s\\ndiv=%.3f\", style=solid,"
              " height=%.1f, width=%.1f, color=\"%.2f,0.15,0.95\"", s0, s0,
              assigncount[s0], pValueStr(sigmage_to_logpvalue(mutsigma[s0])),
              mutrate_CpGMod[s0], diameter, diameter, hue);

    fprintf(fp, "];\n");
  } // if ( fp )

  numdone = 1;
  while (numdone < S)
  {
    thiss = -1;
    thist = -1;
    mindist = 1000000000.0;
    for (s = 0; s < S; s++)     /* want s done */
    {
      if (done[s] == 0)
        continue;
      for (t = 0; t < S; t++)   /* want t not done */
      {
        if (done[t] == 1)
          continue;             /* guarantees t != s */
        if (distance[s][t] < mindist)
        {
          mindist = distance[s][t];
          thiss = s;
          thist = t;
        }
      }
    }
    if (thiss < 0)
    {
      printf("OOPS numdone=%d thiss<0\n", numdone);
      exit(1);
    }
    s = thiss;
    t = thist;

    done[t] = 1;
    numdone++;
    parent[t] = s;

    if ( fp ) 
    {
      // Circle area is proportional to the size of the subfamily 
      diameter =
        (double) 2.0 *sqrt((double) assigncount[t] / (double) maxCount);
  
      // Color heatmap ranges from 0-0.6 ( or 0-216 degrees )
      if (useOriginalGraphColors)
        hue = ((mutrate_CpGMod[t] - minMutRate) * (double) 0.6) / mutRateRange;
      else
        hue =
          ((double) 1.0 -
           ((mutrate_CpGMod[t] - minMutRate) / mutRateRange)) * (double) 0.6;

      if (t <= lastScaffoldIndex)
        fprintf(fp,
                "        %d [shape=circle, label=\"sub%d\\nc=%dpv=%s\\ndiv=%.3f\", style=filled,"
                " height=%.1f, width=%.1f, color=\"%.2f,0.9,0.6\"", t, t,
                assigncount[t], pValueStr(mstLogPValues[t]), mutrate_CpGMod[t],
                diameter, diameter, hue);
      else
        fprintf(fp,
                "        %d [shape=circle, label=\"sub%d\\nc=%dpv=%s\\ndiv=%.3f\", "
                "style=solid, height=%.1f, width=%.1f, color=\"%.2f,0.15,0.95\"",
                t, t, assigncount[t],
                pValueStr(sigmage_to_logpvalue(mutsigma[t])), mutrate_CpGMod[t],
                diameter, diameter, hue);
      fprintf(fp, "];\n");
      fprintf(fp, "        %d -> %d;\n", s, t);
    } // if ( fp )
    edges[s][t] = 1;
    isleaf[s] = 0;
  } // while ( ...

  // Complete file output if necessary
  if ( fp ) 
  {
    fprintf(fp, "}\n");
    fclose(fp);
  }

  nleaves = 0;
  for (s = 0; s < S; s++)
    nleaves += isleaf[s];
  fprintf(outFP, "%d subfamilies overall (%d leaves in tree)\n", S, nleaves);
  printf("%d subfamilies overall (%d leaves in tree)\n", S, nleaves);

  // RMH: Sanity Check
  int rootCnt = 0;
  for(s=0; s<S; s++) 
  {
    if ( parent[s] == -1 )
    {
      printf("Tree root = subfamily-%d\n", s );
      rootCnt++;
    }
  }
  if ( rootCnt > 1 )
  {
    printf("Error: tree has more than one root!!!\n");
    exit(1);
  }

} // build_MST_full(...


// from singlemut
double
sigmage_to_logpvalue(double sigmage)
{
  double GRANULARITY;
  double ans, x, xavg, pvalue;

  GRANULARITY = 0.00001;
  /*
     compute integral from 0 to SIGMA of exp(-x*x/2) 
   */
  ans = 0.0;
  for (x = sigmage; x < sigmage + 10.0; x += GRANULARITY)
  {
    xavg = x + 0.5 * GRANULARITY;
    ans += GRANULARITY * exp(-(xavg * xavg - sigmage * sigmage) / 2.0);
  }
  ans *= 0.39894228;            /* 1/sqrt(2*pi) */
  //ans *= 4.0 * ((double)conLen) * ((double)MAXS); // Bonferroni corrections 
  ans *= sngBonferroni;         /* Bonferroni corrections */
  ans = log(ans);
  ans -= sigmage * sigmage / 2.0;
  return ans;
} // sigmage_to_logpvalue(...

// from singlemut obviously 
// TODO: merger with mst?
void
build_singlemut_MST()
{
  int s, t, s0, done[S], numdone, thiss, thist, mindist, x, a;
  double sage, s0age, pvalue;
  double emutfrac, sigma;

  for (s = 0; s < S; s++)
  {
    for (t = 0; t < S; t++)
      edges[s][t] = 0;
  }

  for (s = 0; s < S; s++)
  {
    for (t = s + 1; t < S; t++)
    {
      distance[s][t] = compute_distance(s, t, SINGLE_MUTATION_INSPENALTY);
      distance[t][s] = distance[s][t];
    }
  }

  /*
     choose the oldest subfamily s0 
   */
  for (s = 0; s < S; s++)
    done[s] = 0;
  s0 = 0;
  s0age = 0.0;
  for (s = 0; s < S; s++)
  {
    if (assigncount[s] == 0)
    {
      printf("OOPS build_MST_full s=%d assigncount=0\n", s);
      exit(1);
    }
    sage = 0.0;
    for (x = 0; x < conLen; x++)
    {
      sage +=
        1.0 -
        ((double) count[s][x][pattern[s][x]]) / ((double) assigncount[s]);
    }
    if (sage > s0age)
    {
      s0age = sage;
      s0 = s;
    }
  }
  done[s0] = 1;
  numdone = 1;
  for (s = 0; s < S; s++)
    parent[s] = -1;

  while (numdone < S)
  {
    thiss = -1;
    thist = -1;
    mindist = 1000000000.0;
    for (s = 0; s < S; s++)     /* want s done */
    {
      if (done[s] == 0)
        continue;
      for (t = 0; t < S; t++)   /* want t not done */
      {
        if (done[t] == 1)
          continue;             /* guarantees t != s */
        if (distance[s][t] < mindist)
        {
          mindist = distance[s][t];
          thiss = s;
          thist = t;
        }
      }
    }
    if (thiss < 0)
    {
      printf("OOPS numdone=%d thiss<0\n", numdone);
      exit(1);
    }
    s = thiss;
    t = thist;

    /*
       about to add edge from s to t.  But first, check if shd fill 
     */
    if (S < MAXS)
    {
      bestmutsigma = sigmaThreshold;
      /*
         Look for mut from s to t 
       */
      for (x = 0; x < conLen; x++)
      {
        if (pattern[s][x] == pattern[t][x])
          continue;
        if (DISALLOW_INDEL)
        {
          if ((pattern[s][x] == 4) || (pattern[t][x] == 4))
            continue;
        }
        a = pattern[t][x];
        if (count[s][x][a] < minCount)
          continue;
        if (DISALLOW_CG)
        {
          if (NEW_CG)
          {
            // CG->TG
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 2))
              continue;
            // TG->CG
            if ((x < conLen - 1)
                && (pattern[s][x] == 3) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 1) && (pattern[t][x + 1] == 2))
              continue;
            // CG->CA
            if ((x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 0))
              continue;
            // CA->CG
            if ((x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 0)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 2))
              continue;
            // CA->TG both columns  THIS IS NEW
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 0)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 2)
                || (x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 0)
                && (pattern[t][x - 1] == 3) && (pattern[t][x] == 2))
              continue;
            // TG->CA both columns  THIS IS NEW
            if ((x < conLen - 1)
                && (pattern[s][x] == 3) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 1) && (pattern[t][x + 1] == 0)
                || (x > 0)
                && (pattern[s][x - 1] == 3) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 0))
              continue;
            // CA->TA ( more common transition )
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 0)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 0))
              continue;
            // TG->TA ( more common transition )
            if ((x > 0)
                && (pattern[s][x - 1] == 3) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 3) && (pattern[t][x] == 0))
              continue;
          }
          else
          {
            // CG->TG
            if ((pattern[s][x] == 1) && (a == 3) && (x < conLen - 1)
                && (pattern[s][x + 1] == 2))
              continue;
            // TG->CG
            if ((pattern[s][x] == 3) && (a == 1) && (x < conLen - 1)
                && (pattern[s][x + 1] == 2))
              continue;
            // CG->CA
            if ((pattern[s][x] == 2) && (a == 0) && (x > 0)
                && (pattern[s][x - 1] == 1))
              continue;
            // CA->CG
            if ((pattern[s][x] == 0) && (a == 2) && (x > 0)
                && (pattern[s][x - 1] == 1))
              continue;

            // CA->TA ( more common transition )
            if ((pattern[s][x] == 1) && (a == 3) && (x < conLen - 1)
                && (pattern[s][x + 1] == 0))
              continue;
            // TG->TA ( more common transition )
            if ((pattern[s][x] == 2) && (a == 0) && (x > 0)
                && (pattern[s][x - 1] == 3))
              continue;
          }
        }


        emutfrac = age[s] * mutperage[pattern[s][x]][x][a];

        if ((emutfrac < 0.0) || (emutfrac > 1.0))
        {
          printf
            ("OOPS emutfrac=%f age[%d]=%f mutperage[ pattern[%d][%d] = %d ][%d][%d] = %f\n",
             emutfrac, s, age[s], s, x, pattern[s][x], x, a,
             mutperage[pattern[s][x]][x][a]);
          //exit(1);
          if (emutfrac < 0)
            emutfrac = 0;
          if (emutfrac > 1.0)
            emutfrac = 1.0;
        }
        sigma = compute_sigma(assigncount[s], count[s][x][a], emutfrac);
        //printf("Pos=%d %c -> %c %lf\n", x, num_to_char( pattern[s][x] ), num_to_char( a ), sigma );
        if (sigma <= bestmutsigma)
          continue;
        //printf("Winner!\n" );
        /*
           WE HAVE A NEW WINNER! 
         */
        bestmutsigma = sigma;
        bestmutSS = s;
        bestmutx = x;
        bestmuta = a;
      }
      /*
         Look for mut from t to s 
       */
      for (x = 0; x < conLen; x++)
      {
        if (pattern[s][x] == pattern[t][x])
          continue;
        if (DISALLOW_INDEL)
        {
          if ((pattern[s][x] == 4) || (pattern[t][x] == 4))
            continue;
        }
        a = pattern[s][x];
        if (count[t][x][a] < minCount)
          continue;
        if (DISALLOW_CG)
        {
          if (NEW_CG)
          {
            // CG->TG
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 2))
              continue;
            // TG->CG
            if ((x < conLen - 1)
                && (pattern[s][x] == 3) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 1) && (pattern[t][x + 1] == 2))
              continue;
            // CG->CA
            if ((x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 0))
              continue;
            // CA->CG
            if ((x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 0)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 2))
              continue;
            // CA->TG both columns
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 0)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 2)
                || (x > 0)
                && (pattern[s][x - 1] == 1) && (pattern[s][x] == 0)
                && (pattern[t][x - 1] == 3) && (pattern[t][x] == 2))
              continue;
            // TG->CA
            if ((x < conLen - 1)
                && (pattern[s][x] == 3) && (pattern[s][x + 1] == 2)
                && (pattern[t][x] == 1) && (pattern[t][x + 1] == 0)
                || (x > 0)
                && (pattern[s][x - 1] == 3) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 1) && (pattern[t][x] == 0))
              continue;
            // CA->TA ( more common transition )
            if ((x < conLen - 1)
                && (pattern[s][x] == 1) && (pattern[s][x + 1] == 0)
                && (pattern[t][x] == 3) && (pattern[t][x + 1] == 0))
              continue;
            // TG->TA ( more common transition )
            if ((x > 0)
                && (pattern[s][x - 1] == 3) && (pattern[s][x] == 2)
                && (pattern[t][x - 1] == 3) && (pattern[t][x] == 0))
              continue;
          }
          else
          {
            if ((pattern[t][x] == 1) && (a == 3) && (x < conLen - 1)
                && (pattern[t][x + 1] == 2))
              continue;
            if ((pattern[t][x] == 3) && (a == 1) && (x < conLen - 1)
                && (pattern[t][x + 1] == 2))
              continue;
            if ((pattern[t][x] == 2) && (a == 0) && (x > 0)
                && (pattern[t][x - 1] == 1))
              continue;
            if ((pattern[t][x] == 0) && (a == 2) && (x > 0)
                && (pattern[t][x - 1] == 1))
              continue;
            if ((pattern[t][x] == 1) && (a == 3) && (x < conLen - 1)
                && (pattern[t][x + 1] == 0))
              continue;
            if ((pattern[t][x] == 2) && (a == 0) && (x > 0)
                && (pattern[t][x - 1] == 3))
              continue;
          }
        }

        emutfrac = age[t] * mutperage[pattern[t][x]][x][a];
        if ((emutfrac < 0.0) || (emutfrac > 1.0))
        {
          printf("OOPS emutfrac=%f\n", emutfrac);
          exit(1);
        }
        sigma = compute_sigma(assigncount[t], count[t][x][a], emutfrac);
        //printf("Pos=%d %c -> %c %lf\n", x, num_to_char( pattern[t][x] ), num_to_char( a ), sigma );
        if (sigma <= bestmutsigma)
          continue;
        //printf("WinnerHERE!\n" );
        /*
           WE HAVE A NEW WINNER! 
         */
        bestmutsigma = sigma;
        bestmutSS = t;
        bestmutx = x;
        bestmuta = a;
      }
      if (bestmutsigma > sigmaThreshold)        /* Found a fill vertex */
      {
        //printf( "Going to do it: [parent %d] %d %c -> %c %.2f\n", bestmutSS, bestmutx, num_to_char( pattern[bestmutSS][bestmutx] ), num_to_char( bestmuta ), bestmutsigma );
        mutsigma[S] = bestmutsigma;
        build_new_singlemut_subfamily();
        build_singlemut_MST();
        return;
      }
    }

    done[t] = 1;
    numdone++;
    parent[t] = s;
    for (s0 = 0; s0 < S; s0++)
      labels[s0] = 0;
    labels[s] = 1;
    labels[t] = 1;
    edges[s][t] = 1;
  }
}


// from singlemut obviously
void
build_new_singlemut_subfamily()
{
  int n, s, thisx, thisa, x, a, SS, besta, bestb, b;
  int ccount[5], ccounti[2];
  int oldcount, newcount;
  int totalcount, count1;
  double emutfrac;
  double bestscore, countratio;

  /*
     we are given bestmutSS, bestmutx, bestmuta 
   */
  SS = bestmutSS;
  thisx = bestmutx;
  thisa = bestmuta;

  totalcount = assigncount[SS];
  count1 = count[SS][thisx][thisa];
  emutfrac = age[SS] * mutperage[pattern[SS][thisx]][thisx][thisa];
  fprintf(outFP,
          "Building subfamily %d (parent %d, logpvalue %f, sigma %.2f): pos %d %c to %c (%d %d %f)\n",
          S, SS, sigmage_to_logpvalue(bestmutsigma), bestmutsigma, thisx,
          num_to_char(pattern[SS][thisx]), num_to_char(thisa), totalcount,
          count1, emutfrac);
  printf
    ("  Building subfamily %d (parent %d, logpvalue %f, sigma %.2f): pos %d %c to %c (%d %d %f)\n",
     S, SS, sigmage_to_logpvalue(bestmutsigma), bestmutsigma, thisx,
     num_to_char(pattern[SS][thisx]), num_to_char(thisa), totalcount,
     count1, emutfrac);

  for (x = 0; x < conLen; x++)
  {
    for (a = 0; a < 5; a++)
      count[S][x][a] = 0;
  }

  for (x = 0; x < conLen; x++)
  {
    pattern[S][x] = pattern[SS][x];
    patterni[S][x] = patterni[SS][x];
  }
  pattern[S][thisx] = thisa;
  S++;

  for (s = 0; s < S; s++)
    mark[s] = 0;
  mark[SS] = 1;
  mark[S - 1] = 1;

  pattern_to_assign_mark(thisx);
  assign_to_pattern_mark();

}


//
// computeMutationRatesKimura()
//
//   Calculate the Kimura substition divergence for each 
//   subfamily.
//
//   DNA Transitions/Transversions
//
//         A -V- C
//         |\   /|
//         I  V  I
//         |/   \|
//         G -V- T
//
//   CpG: Transition mutations are 10-15x more likely 
//        at CpG sites in many mammalian species.
//
//             CpG---deamination-->TG
//             CpG---deamination-->CA
//
//    Arian's scoring method:
//       CpG to:
//          AA = Vi   CA = i   GA = Vi   TA = I
//          AC = VV   CC = V   GC = VV   TC = iV
//          AG = V    CG =     GG = V    TG = i
//          AT = VV   CT = V   GT = VV   TT = iV
//
//     Where:  V = 1    - Transversion
//             i = 1/10 - Transition
//             I = 1    - Transition
//
//     **SPECIAL IMPLEMENTATION**
//     In addition the consnensus sequence is crudely 
//     generated ( ie. CpG sites are not inferred ).
//     Therefore it makes sense to also treat CA and
//     TG sites as possible CpG sites and treat 
//     TG->CA and CA->TG as a single transition event.
//
//     Then plug in the integer transversion count and real
//     transition counts into the 2-parameter Kimura divergence
//     formula.
//
//          K = - 1/2 ln(( 1 - 2p -q) * sqrt(1-2q))
//        Where:
//          p = transition proportion
//          q = transversion proportion
//          
//     Alternatively Arian reports the simple divergence throwing
//     out single transitions at CpG sites and collapsing two
//     transitions at a CpG site into one. 
//
//
//   Previous implementation ( and original ) calculated mutrate
//   as: 
//
//      - In build_scaffold.c if the consensus is CG/TG/CA do not 
//        record any mutation at either of the bases.
//      - In fill_scaffold.c if the consensus is CG/TG/CA do not 
//        record any mutation at the first base and the skip the second
//        entirely ( mutrate[] ) or skip both sites ( mutrate2[] ).
//        Reported in subfamilies.seq as "mutrate mutrate[]/mutrate2[]".
//      - Include deletions individually and insertion of any size
//        as a single penalty.
//
//   We no longer include insertions/deletions in the calculation except
//   to invalidate a CG/TG/CA in the consensus if there is an insertion
//   likely ( patterni[][] ) in between the two bases.
//   
// Updated 11/2016 RMH
//
// NOTE: Expects two arrays of MAXS size
void
computeMutationRatesKimura(double *mutrate, double *mutrate_CpGMod)
{
  int ccount[MAXS];
  int trans[MAXS], transv[MAXS];
  double trans_CpGMod[MAXS];
  double logOperand, kimura;
  int n, s, x, totcount;
  double p, q;

  // Initialize counters
  for (s = 0; s < S; s++)
  {
    ccount[s] = 0;
    trans[s] = 0;
    trans_CpGMod[s] = 0;
    transv[s] = 0;
  }
  totcount = 0;

  for (n = 0; n < N; n++)
  {
    s = assign[n];
    if (s < 0)
      continue;
    totcount++;
    ccount[s] += 1;
    for (x = 0; x < conLen; x++)
    {
      // Is this a CpG site? NOTE: Not considering sites with intervening insertions
      //     0=a,1=c,2=g,3=t,4=-
      if ( x < conLen - 1 && patterni[s][x] != 1 )
      {
        if ( (pattern[s][x] == 1) && (pattern[s][x + 1] == 2) )  // CpG
        {
          // Transversions = 0,1,2 Transitions = 0, 1/10, 1
          int tmpTrans = 0;

          // C Mutations
          if ( ele[n][x] == 3 )  // C->T Transition
             tmpTrans++;
          if ( ele[n][x] == 0 || ele[n][x] == 2 ) // C->A, C->G Transversions
             transv[s]++;

          // G Mutations
          if ( ele[n][x+1] == 0 )  // G->A Transition
             tmpTrans++;
          if ( ele[n][x+1] == 1 || ele[n][x+1] == 3 ) // G->C, G->T Transversions
             transv[s]++;

          if ( tmpTrans == 2 )
            trans_CpGMod[s] += 1;
          else if ( tmpTrans == 1 )
            trans_CpGMod[s] += 0.1;

          trans[s] += tmpTrans;

          x++;
          continue;
        }else if ( (pattern[s][x] == 3) && (pattern[s][x + 1] == 2) )  // TpG
        {
          // CA = 1 transition -- all other normal accounting
          if ( ele[n][x] == 1 || ele[n][x+1] == 0 )
          {
            trans_CpGMod[s] += 1;
            trans[s] += 2;
            x++;
            continue;
          }
        }else if ( (pattern[s][x] == 1) && (pattern[s][x + 1] == 0) )  // CpA
        {
          // TG = 1 transition -- all other normal accounting
          if ( ele[n][x] == 3 || ele[n][x+1] == 2 )
          {
            trans_CpGMod[s] += 1;
            trans[s] += 2;
            x++;
            continue;
          }
        }
      }
           
      if ( ( pattern[s][x] == 0 && ele[n][x] == 2 ) || // A->G trans
           ( pattern[s][x] == 1 && ele[n][x] == 3 ) || // C->T trans
           ( pattern[s][x] == 2 && ele[n][x] == 0 ) || // G->A trans
           ( pattern[s][x] == 3 && ele[n][x] == 1 ) )  // T->C trans
      {
             trans_CpGMod[s]++;
             trans[s]++;
      }
      else if ( pattern[s][x] != ele[n][x] )
      {
             transv[s]++;
      }
    }
  }
  for (s = 0; s < S; s++)
  {
    if (ccount[s] == 0)
    {
      mutrate[s] = 0.0;
      mutrate_CpGMod[s] = 0.0;
      continue;
    }

    p = (double)trans[s] / ( ccount[s] * conLen );
    q = (double)transv[s] / ( ccount[s] * conLen );
    kimura = 0.0;
    logOperand = ( ( 1 - ( 2 * p ) - q ) * 
                   pow(( 1 - ( 2 * q ) ),0.5 ) );
    if ( logOperand > 0.0 )
      kimura = fabs( ( -0.5 * log( logOperand ) ) );

    mutrate[s] = kimura;

    kimura = 0.0;
    q = (double)trans_CpGMod[s] / ( ccount[s] * conLen );
    logOperand = ( ( 1 - ( 2 * p ) - q ) * 
                   pow(( 1 - ( 2 * q ) ),0.5 ) );
    if ( logOperand > 0 )
      kimura = fabs( ( -0.5 * log( logOperand ) ) );

    mutrate_CpGMod[s] = kimura;
  }  
} // computeMutationRatesKimura()


//
// print_subfamilies()
//
//   Generate the final *.subfamilies output file.
//
void
print_subfamilies(int S)
{
  int n, s, x, a, b, freq, maxfreq, totcount, t, r;
  int rdist[MAXS], sdist[MAXS], hit[MAXS], mindist;
  int nmut[MAXS], ntot[MAXS], nmut2[MAXS], ntot2[MAXS];
  double mutrate[MAXS], mutrate2[MAXS];
  FILE *fp;


  // NOTE: Do not need to calculate mindist anymore

  if ((fp = fopen(subfamFile, "w")) == NULL)
  {
    printf("Could not open input file %s\n", subfamFile);
    exit(1);
  }

  computeMutationRatesKimura(mutrate, mutrate2);

  /*
     print subfamilies 
   */
  if (S == MAXS)
  {
    printf("OOPS S=MAXS=%d\n", S, MAXS);
    exit(1);
  }

  for (s = 0; s < S; s++)
  {
    if (s < numScaffolds)
      fprintf(fp, "Subfamily %d: count %d mutrate %.03f/%.03f "
              "mstLogPValue %f\n", s, assigncount[s], mutrate[s],
              mutrate2[s], mstLogPValues[s]);
    else
      fprintf(fp,
              "Subfamily %d: count %d mutrate %.03f/%.03f sigma %f, logpvalue %f\n",
              s, assigncount[s], mutrate[s], mutrate2[s], mutsigma[s],
              sigmage_to_logpvalue(mutsigma[s]));


    for (x = 0; x < conLen; x++)
    {
      if (pattern[s][x] != Sxsequence[x])
        fprintf(fp, "%d:%c ", x, num_to_char(pattern[s][x]));
      if (patterni[s][x] == 1)
        fprintf(fp, "%d:+ ", x);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

}
