/*
  PDF.H - Header file for use with PDF routines.
          Copyright 1999 by Michael Wall and Rice University.  

	  Distribution and modification of this package is governed
	  under the "Artistic License" (perl license) described in the
	  file LICENSE in the SOESA distribution.

  Author: Michael Wall
  Date: 11/19/99
  Version: 0.2
  Date: 2/11/1999
  Version: 0.1

*/

/*
 * These are standard header files that are required:
 */

#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<signal.h>
#include<unistd.h>
#include<ctype.h>
#include<fcntl.h>
#include<stropts.h>

/*
 * Commands recognized through the FIFO in servermode:
 */
#define TNT_LONG_CMD         "long"
#define TNT_SHORT_CMD      "short"
#define EVAL_CMD "eval"
#define XPLOR_CMD "xplor"
#define SHUTDOWN_CMD "shutdown"

#define FIFO_FILENAME "soesa.fifo" /* Default name of the FIFO */
#define ZERO 0.0
#define MSG_SIZE 256 /* Bytes allocated for the message passed through the FIFO */
#define INPUT_LINE_SIZE 120 /* # of bytes alloc'd for a single input line */
#define CMD_SIZE 80 /* # of bytes alloc'd for a single command word */
#define MAX_NATOM 10000 /* Maximum # of atoms in a structure */
#define FNAME_SIZE 120 /* # of bytes alloc'd for a filename */
#define TYPE_SIZE 5 /* # of bytes alloc'd for a single atom type label */
#define RESID_SIZE 5 /* # of bytes alloc'd for a single rsidue type label */
#define SEGID_SIZE 5 /* # of bytes alloc'd for a single segment name */
#define RESN_SIZE 5 /* # of bytes alloc'd for a residue number (#digits+1) */
#define MAX_AA 50 /* Maximum # amino acid types */
#define MAX_AA_NATOM 20 /* Maximum # of atoms in an amino acid */
#define MAX_M 6 /* Maximum value of m capable of being flagged via mexc, minc */
#define DEFAULT_RCUT 10.0 /* Dflt atom pair cutoff distance for PDF calculation */
#define IND_HEADER_LENGTH 20 /* # bytes in the pdf data structure header */
#define PDF_HEADER_LENGTH 36 /* #bytes in the pdf distn header */

/*
 * Types for client/server support:
 */

typedef enum {SHORT,LONG,XPLOR,EVAL,SHUTDOWN} CMD; /* Recognized commands */

/* Data structure for message passing: */
typedef struct {
  pid_t pid; /* PID of the calling process */
  CMD cmd; /* Command */
} MSG_STRUCT;

/*
 * The following data structures are derived from specifications of the PDF
 *   database described by Atipat Rojnuckarin:
 */

/*
 * Amino acid type:
 */

typedef struct {
  char AAcode1; /* Single-letter AA code */
  char AAcode3[RESID_SIZE]; /* 3-letter AA code */
  char *atname[MAX_AA_NATOM]; /* List of atom names in the amino acid */
  int natom; /* Number of atoms in the amino acid */
  int baseidx; /* Position in the amino acid index */
} AA_STRUCT;

/*
 * Data structure for a single PDF distn: 
 */

typedef struct {
  int dist; /* The 'm' value = number of peptide bonds separating the atoms */
  int res1; /* Residue type code for the first atom */
  int atm1; /* Atom type code for the first atom */
  int res2; /* Residue type code for the second atom */
  int atm2; /* Atom type code for the second atom */
  float p_ref; /* Reference value for the distribution values */
  float x_start; /* Smallest distance in the distn */
  float x_step; /* Distance bin size */
  int n_data; /* Number of data points */
  float *pdf; /* Array of data points */
} PDF_STRUCT;

/*
 * Data structure for the whole PDF database:
 */

typedef struct {
  int version; /* Version number of the database */
  int NLevels; /* Number of 'm' levels, including the tertiary level */
  int NPoints; /* Number of points in each distn */
  float Cutloc; /* Dist. cutoff for including local interactions */
  float Cutter; /* Dist. cutoff for including tertiary interactions */
  PDF_STRUCT *pdfpair; /* Array of pdf's */
} IND_STRUCT;

/*
 * Data structure for energies:
 */

typedef struct {
  float energy;
  float drvx; /* X-derivative */
  float drvy; /* Y-derivative */
  float drvz; /* Z-derivative */
  float ddrvx; /* X-curvature */
  float ddrvy; /* Y-curvature */
  float ddrvz; /* Z-curvature */
  float n; /* Number of data points used to calculate the energy */
} ENERGY_STRUCT;

/*
 * Data structure for PDB files:
 *   (not used as of 11/8/99 - MW)
 */

typedef struct {
  char **type;
  char **resid;
  char *chain;
  char **segid;
  int *resn;
  float *x;
  float *y;
  float *z;
  float *q;
  float *b;
  int natom;
} PDB_STRUCT;

/*
 * I/O subroutines:
 */

/* Get the size in bytes of a file from the file system: */
off_t getFileSize(const char *path);

/* Read a PDB file: */
int readPDB(char **type, char **resid, char *chain, int *resn, 
	    float *x, float *y, float *z, float *q, float *b, 
	    char **segid, int *natom, char *fname);

/* Read an amino acid descriptions file: */
int readAA(AA_STRUCT *aa, int *naa, int *natype, char *fname);

/* Read a file into a buffer: */
int readFile(const char *fname, char *buf, off_t fsize);

/* Write a file from a buffer: */
int writeFile(const char *fname, char *buf, off_t fsize);

/* Write the energy, derivs in TNT long format: */
int writeTNTlong(ENERGY_STRUCT *e, float *escale, int *natom,
		 char *pdb_fname, char *out_fname);

/* Write the energy in TNT short format: */
int writeTNTshort(ENERGY_STRUCT *e, float *escale, int *natom, char *out_fname);

/* Get the index associated with a 3-letter amino acid type: */
int getAA3idx(char *AAcode3,AA_STRUCT *aa,int *naa);

/* Get the index associated with an atom type plus a 3-letter amino acid type: */
int getAtomidx(char *AAcode3,char *atname,AA_STRUCT *aa,int *naa);

/* Calculate the hash key associated with a pair of atom types, */
/*   residue types plus 'm': */
long getPDFhash(int *m,char *resid1,char *type1,char *resid2,char *type2,
		AA_STRUCT *aa,int *naa,int *natype, int *TLevel);

/* Get the PDF score at a given distance: */
void getPDFVal(float *xlo, float *xhi,float *ya,float *y2a,
	       int *n,float *x,float *y);

/* Get the derivative of the PDF distn at a given distance: */
void getPDFDrv(float *xlo, float *xhi,float *ya,float *y2a,
	       int *n,float *x,float *y);

/* Get the 2d derivative of the PDF distn at a given distance: */
void getPDFCrv(float *xlo, float *xhi,float *ya,float *y2a,
	       int *n,float *x,float *y);

/* Reverse the byte order in a buffer (assume 4-byte values): */
int xformRevByte(char *buf, off_t bufsize);

/* Calculate the spline coefficients for all of the PDF data: */
int xformCalcSpline(IND_STRUCT *ind, off_t ind_fsize);

/* Output a requested PDF distn: */
int showPDF(char *fname, 
	    int *m, char *resid1, char *type1, char *resid2, char *type2,
	    char *zoom, IND_STRUCT *indspl, int *idx, 
	    AA_STRUCT *aa, int *naa, int *natype, int *natom);

/* Calculate the average PDF energy by residue: */
int calcAvgByRes(int *resn, ENERGY_STRUCT *energy, int *natom,
		 float *eval, int *neval);

/* Calculate the PDF energy of a structure for TNT, EVAL-type calcs: */
int calcEnergy(char **segid,char *chain,int *resn, char **resid, char **type,
	       float *x, float *y, float *z, 
	       IND_STRUCT *indspl, int *idx,
	       AA_STRUCT *aa, int *naa, int *natype, int *natom,
	       ENERGY_STRUCT *e, char *cmd, char *msel);

/* Calculate the PDF energy of a structure for XPLOR-type calcs: */
int calcEnergyXPLOR(char *chain,int *resn, char **resid, char **type,
	       float *x, float *y, float *z, 
	       IND_STRUCT *indspl, int *idx,
	       AA_STRUCT *aa, int *naa, int *natype, int *natom,
	       ENERGY_STRUCT *e, char *msel,
		    int *dselect, float *rcut);


