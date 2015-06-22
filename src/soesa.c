/*
  SOESA.C - Structure Optimization and Evaluation using Separations of Atoms
            Copyright 1999 by Michael Wall and Rice University

	    Distribution and modification of this package is governed
	    under the "Artistic License" (perl license) described in the
	    file LICENSE in the SOESA distribution.

  Author: Michael Wall
  Version: 0.2  SOESA  11/6/1999
     Based on OVID by Michael Wall
  Version: 0.1  OVID  5/19/1999
     Based on pdftst.c 2/11/1999 by Michael Wall
*/

#include<pdf.h>


int main(int argc, char *argv[])
{

  /*
   * File pointer used for i/o:
   */

  FILE *fptr;

  /*
   * Used in parsing command line arguments:
   */

  char *cmd; /* Points to a command from the argument list */
  int inpct=0; /* Keeps track of how many arguments are passed */
	     /* to a command line directive */

  /*
   * Miscellaneous file names:
   */

  char 
    *pdb_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* PDB file */
    *idx_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* PDF hash table */
    *ind_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* PDF data file */
    *spl_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* Spline coeffs. */
    *indspl_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /*Data + Spline*/
    *aacfg_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* Amino defs */
    *out_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* Output file */
    *fifo_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* FIFO name */
    *full_fname1, /* Used to hold file name */
    *full_fname2; /* Used to hold file name */

  size_t
    pdfsize,
    datasize,
    idx_offset,
    num_wrote;
  
  /*
   * File sizes reported by file system calls:
   */

  off_t
    idx_fsize, /* Size in bytes of PDF hash file */
    ind_fsize, /* Size in bytes of PDF data file */
    spl_fsize, /* Size in bytes of PDF spline coefficients file */
    indspl_fsize; /* Size in bytes of PDF data + spline coefficients file */


  /* 
   * Variables used to display the AA sequence of the PDB file:
   */

  int
    resn_ct, /* Total number of residues currently displayed */
    resn_old, /* Number of the last residue that was displayed */
    aaidx; /* Amino acid index used for displaying the sequence */


  /* 
   * Data structures for accessing the PDF data:
   */

  long PDFhash; /* Hash code to look up PDF data; */
                /* calculated from atom pair info. */

  int *idx; /* Array of PDF hash values */

  IND_STRUCT
    *ind, /* PDF data */
    *spl, /* PDF spline coefficient data */
    *indspl; /* Combined PDF data and spline coefficients */

  void
    *idx_buf, /* Buffer for PDF hash values */
    *ind_buf, /* Buffer for PDF data */
    *spl_buf, /* Buffer for PDF spline coefficients */
    *indspl_buf; /* Buffer for combined PDF data and spline coefficients */

  int TLevel; /* Indicates the lowest value of 'm' for which  */
              /* the interaction is to be deemed tertiary */


  /*
   * The following are used in function calls to retrieve a PDF distribution:
   */
 
  char
    *resid1=(char *)calloc(RESID_SIZE+1,sizeof(char)), /* Residue name */
    *resid2=(char *)calloc(RESID_SIZE+1,sizeof(char)), /* Residue name */
    *type1=(char *)calloc(TYPE_SIZE+1,sizeof(char)), /* Chemical type */
    *type2=(char *)calloc(TYPE_SIZE+1,sizeof(char)); /* Chemical type */

  int m; /* Parameter 'm' indicates # of peptide bonds between the atoms */

  /*
   * Used to manipulate amino acid info:
   */

  AA_STRUCT
    *aa; /* List of amino acid names, atom type composition, etc. */

  int
    naa, /* Total number of amino acid types */
    natype; /* Total number of atom types */ 

  /*
   * These variables and arrays store the PDB file information:
   */

  char 
    **type = (char **)calloc(MAX_NATOM,sizeof(char *)), /* Chemical types */
    **resid = (char **)calloc(MAX_NATOM,sizeof(char *)), /* Residue names */
    *chain = (char *)calloc(MAX_NATOM,sizeof(char)), /* Chain ID */
    **segid = (char **)calloc(MAX_NATOM,sizeof(char *)); /* Segment names */

  float
    *x = (float *)calloc(MAX_NATOM,sizeof(float)), /* X coordinates */
    *y = (float *)calloc(MAX_NATOM,sizeof(float)), /* Y coordinates */
    *z = (float *)calloc(MAX_NATOM,sizeof(float)), /* Z coordinates */
    *q = (float *)calloc(MAX_NATOM,sizeof(float)), /* Occupancies */
    *b = (float *)calloc(MAX_NATOM,sizeof(float)); /* B-factors */

  int
    *resn = (int *)calloc(MAX_NATOM,sizeof(int)); /* Residue numbers */

  int natom=0; /* Total number of atoms in the structure */

  /*
   * Variables for use in generating a PDF score by residue ('eval'):
   */

  int
    *neval; /* For each residue, the number of atom pairs that were used */
            /* to calculate the average score */

  float *eval; /* Array of PDF scores for each residue */

  /*
   * Generic counting variables:
   */

  int
    i,j,k;

  /*
   * The following are used in energy calculation:
   */

  ENERGY_STRUCT
    *e = (ENERGY_STRUCT *)calloc(MAX_NATOM,sizeof(ENERGY_STRUCT)); 
                       /* Calculated energies and derivatives */

  int *dselect; /* Double selection array for X-PLOR calculations */

  float
    escale=0.58, /* Overall energy scale */
    pdfe, /* Total energy */
    rcut=DEFAULT_RCUT; /* Distance cutoff beyond which no */
		       /* PDF energy is calculated */

  char
    *msel=(char *)calloc(MAX_M+1,sizeof(char)); 
             /* msel is an array of flags to indicate whether */
	     /* to include contributions from atom pairs separated by 'm' */
             /* peptide bonds */

  /*
   * Used for log information:
   */

  time_t
    t; /* Time for generating time stamps */

  /*
   * The following handle client requests received through the FIFO:
   */

  mode_t 
    mode;  /* File mode descriptor */

  int
    fifo_stream; /* The handle used to access the FIFO */

  struct strbuf
    raw_msg; /* Data structure used to read data from the FIFO into msg_buf */

  char
    *msg_buf; /* Buffer for the msg passed through the FIFO in server mode */

  MSG_STRUCT
    *parsed_msg; /* Data structure for parsing the msg read in from the FIFO */

  char servercmd[10]; /* Command string for the requested action */

  char badcmd=0; /* Flag to indicate that an invalid command was received */

  char serverquit=0; /* Flag to indicate server shutdown */

  int
    dummy=0; /* Passed as an unused parameter in FIFO i/o calls */

  char *calc_home; /* Path passed through the FIFO in server mode.  Look for */
                   /* PDB file here, and place output here. */
  size_t 
    path_len_bytes; /* Length of the path name in calc_home */

  /*
   * Command flags:
   */

  char
    debug=0, /* Display debugging info */
    readidx=0, /* Read the PDF hash values from a file */
    readind=0, /* Read the PDF data from the file */
    readspl=0, /* Read the spline coefficients from a file */
    readindspl=0, /* Read combined PDF data + spline coeffs from a file */
    readaa=0, /* Read the amino acid info from the amino acid confg. file */
    showpdf=0, /* Display a particular pdf (which one must be specified) */
    showdrv=0, /* Display the derivative of a particular pdf (ditto) */
    showpdb=0, /* Display the pdb file that has been read in */
    showaa=0, /* Display the information about amino acid definitions */
    showseq=0, /* Display the sequence of the PDB file */
    revidx=0, /* Reverse the byte order of the hash table */
    revind=0, /* Reverse the byte order of the PDF data */
    calcspline=0, /* Calculate spline coefficients for the PDF data */
    merge=0, /* Merge spline coefficients with the PDF data */
    zoompdf=1, /* Zoom factor for displaying a PDF */
    writeidx=0, /* Output the pdf hash table to a file */
    writeind=0, /* Output the pdf data to a file */
    servermode=0, /* Run in server mode, listening at the FIFO for commands */
    TNTlong=0, /* Calculate the energy in TNT long format */
    TNTshort=0, /* Calculate the energy in TNT short format */
    xplor=0, /* Calculate the energy in xplor format */
    evalpdb=0, /* Evaluate the pdf score by residue of the input pdb file */
    verbose=0; /* Run in verbose logging mode */


  /*
   * Allocate chemical type, residue name and segment name arrays for PDB:
   */

  for (i=0;i<MAX_NATOM;i++) {
    type[i] = (char *)calloc(TYPE_SIZE+1,sizeof(char));
    resid[i] = (char *)calloc(RESID_SIZE+1,sizeof(char));
    segid[i] = (char *)calloc(SEGID_SIZE+1,sizeof(char));
  }

  /*
   * By default include all 'm' values in energy calculation:
   */

  for (i=0;i<MAX_M;i++) {
    msel[i]=1;
  }

  /*
   * Parse command line arguments:
   *    These are in form: -command <arg1> <arg2> ...
   */

 --argc; 
  while (argc>0) {
    if (argv[argc][0] != '-') {
      ++inpct;
    }
    else {
      cmd = argv[argc]+1;
      /*
       * PDB file name:
       */
      if (!strcmp(cmd,"pdb")) {
	if (inpct == 1) {
	  strcpy(pdb_fname,argv[argc+1]);
	}
	else {
	  perror("SOESA: -pdb takes 1 argument only.");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Output file name:
       */
      else if (!strcmp(cmd,"out")) {
	if (inpct==1) {
	  strcpy(out_fname,argv[argc+1]);
	}
	else {
	  perror("SOESA: -out takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Hash vlaues file name ("index" file):
       */
      else if (!strcmp(cmd,"hashfile") || !strcmp(cmd,"idx")) {
	if (inpct==1) {
	  strcpy(idx_fname,argv[argc+1]);
	  readidx=1;
	}
	else {
	  perror("SOESA: -idx takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * PDF data file name ("induced" file):
       */
      else if (!strcmp(cmd,"ind")) {
	if (inpct==1) {
	  strcpy(ind_fname,argv[argc+1]);
	  readind=1;
	}
	else {
	  perror("SOESA: -ind takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Spline coefficients file name:
       */
      else if (!strcmp(cmd,"spl")) {
	if (inpct==1) {
	  strcpy(spl_fname,argv[argc+1]);
	  readspl=1;
	}
	else {
	  perror("SOESA: -spl takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * PDF data + spline coefficients file name:
       */
      else if (!strcmp(cmd,"indspl") || !strcmp(cmd,"datafile")) {
	if (inpct==1) {
	  strcpy(indspl_fname,argv[argc+1]);
	  readindspl=1;
	}
	else {
	  perror("SOESA: -indspl takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      else if (!strcmp(cmd,"merge")) {
	if (inpct == 1) {
	  strcpy(indspl_fname,argv[argc+1]);
	  merge=1;
	}
	else {
	  perror("SOESA: merge takes 1 argument");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * FIFO file name:
       */
      else if (!strcmp(cmd,"fifo")) {
	if (inpct==1) {
	  strcpy(fifo_fname,argv[argc+1]);
	}
	else {
	  perror("SOESA: -fifo takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Amino acid info. file name:
       */
      if (!strcmp(cmd,"aacfg")) {
	if (inpct == 1) {
	  strcpy(aacfg_fname,argv[argc+1]);
	  readaa=1;
	}
	else {
	  perror("SOESA: -aacfg takes 1 argument only.");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Display a PDF (selected by "<m> <resid1> <type1> <resid2> <type2>"):
       */
      else if (!strcmp(cmd,"showpdf")) {
	if (inpct==1) {
	  sscanf(argv[argc+1],"%d %s %s %s %s",&m,resid1,type1,resid2,type2);
	  showpdf=1;
	}
	else {
	  perror("SOESA: -showpdf takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Display the derivative of a PDF (as in showpdf):
       */
      else if (!strcmp(cmd,"showdrv")) {
	if (inpct==1) {
	  sscanf(argv[argc+1],"%d %s %s %s %s",&m,resid1,type1,resid2,type2);
	  showdrv=1;
	}
	else {
	  perror("SOESA: -showdrv takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Scale for the sampling frequency in diaplaying PDF, drv:
       */
      else if (!strcmp(cmd,"zoompdf")) {
	if (inpct==1) {
	  zoompdf=atoi(argv[argc+1]);
	}
	else {
	  perror("SOESA: zoompdf takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Cutoff distance beyond which the PDF energy is set to zero:
       */
      else if (!strcmp(cmd,"rcut")) {
	if (inpct==1) {
	  rcut=atof(argv[argc+1]);
	}
	else {
	  perror("SOESA: rcut takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Energy scale:
       */
      else if (!strcmp(cmd,"escale")) {
	if (inpct==1) {
	  escale=atof(argv[argc+1]);
	}
	else {
	  perror("SOESA: escale takes 1 argument only");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Output the PDB file specified by -pdb:
       */
      else if (!strcmp(cmd,"showpdb")) {
	if (inpct == 0) {
	  showpdb=1;
	}
	else {
	  perror("SOESA: showpdb takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Evaluate the PDF score of the PDB file by residue:
       */
      else if (!strcmp(cmd,"evalpdb") || !strcmp(cmd,"eval")) {
	if (inpct == 0) {
	  evalpdb=1;
	}
	else {
	  perror("SOESA: evalpdb takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Show the number of bytes of standard types for debugging:
       */
      else if (!strcmp(cmd,"sizeofs")) {
	if (inpct == 0) {
	  printf("sizeof(*) = *:\nchar = %ld\nshort = %ld\nint = %ld\nlong = %ld\nfloat = %ld\ndouble = %ld\n",
		 sizeof(char),
		 sizeof(short),
		 sizeof(int),
		 sizeof(long),
		 sizeof(float),
		 sizeof(double));
	}
	else {
	  perror("SOESA: sizeofs takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Show the amino acid information for debugging:
       */
      else if (!strcmp(cmd,"showaa")) {
	if (inpct == 0) {
	  showaa=1;
	}
	else {
	  perror("SOESA: showaa takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Display the amino acid sequence of the PDB file:
       */
      else if (!strcmp(cmd,"showseq")) {
	if (inpct == 0) {
	  showseq=1;
	}
	else {
	  perror("SOESA: showseq takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Reverse the byte order of the hash value file:
       */
      else if (!strcmp(cmd,"revidx")) {
	if (inpct == 0) {
	  revidx=1;
	}
	else {
	  perror("SOESA: revidx takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Reverse the byte order of the PDF data file:
       */
      else if (!strcmp(cmd,"revind")) {
	if (inpct == 0) {
	  revind=1;
	}
	else {
	  perror("SOESA: revind takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Output the hash table file:
       */
      else if (!strcmp(cmd,"writeidx")) {
	if (inpct == 0) {
	  writeidx=1;
	}
	else {
	  perror("SOESA: writeidx takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Output the PDF data file:
       */
      else if (!strcmp(cmd,"writeind")) {
	if (inpct == 0) {
	  writeind=1;
	}
	else {
	  perror("SOESA: writeind takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Calculate the spline coefficients for the PDF data:
       *   Results are left in 'ind'
       */
      else if (!strcmp(cmd,"calcspline")) {
	if (inpct == 0) {
	  calcspline=1;
	}
	else {
	  perror("SOESA: calcspline takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Calculate the spline coefficients for the PDF data:
       *   Results are left in 'ind'
       */
      /*
       * Calculate the energies and derivs and output in TNT format:
       */
      else if (!strcmp(cmd,"TNTlong")) {
	if (inpct == 0) {
	  TNTlong=1;
	}
	else {
	  perror("SOESA: TNTlong takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Calculate the total energy and output in TNT format:
       */
      else if (!strcmp(cmd,"TNTshort")) {
	if (inpct == 0) {
	  TNTshort=1;
	}
	else {
	  perror("SOESA: TNTshort takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Calculate the total energy and output in XPLOR format:
       */
      else if (!strcmp(cmd,XPLOR_CMD)) {
	if (inpct == 0) {
	  xplor=1;
	}
	else {
	  perror("SOESA: xplor takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Give verbose logging info:
       */
      else if (!strcmp(cmd,"verbose")) {
	if (inpct == 0) {
	  verbose=1;
	}
	else {
	  perror("SOESA: verbose takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Give debugging info (not much as of 11/7/99):
       */
      else if (!strcmp(cmd,"debug")) {
	if (inpct == 0) {
	  debug=1;
	}
	else {
	  perror("SOESA: debug takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * Operate as an energy "server," listening at the FIFO for commands:
       */
      else if (!strcmp(cmd,"servermode")) {
	if (inpct == 0) {
	  servermode=1;
	}
	else {
	  perror("SOESA: servermode takes 0 arguments");
	  exit(1);
	}
      }
      /*
       * A list of 'm' values to include in energy calculation:
       */
      else if (!strcmp(cmd,"minc")) {
	while (inpct > 0) {
	  msel[atoi(argv[argc+inpct])]=1;
	  inpct--;
	}
      }
      /*
       * A list of 'm' values to exclude in energy calculation:
       */
      else if (!strcmp(cmd,"mexc")) {
	while (inpct > 0) {
	  msel[atoi(argv[argc+inpct])]=0;
	  inpct--;
	}
      }
    }
    --argc;
  }

  fprintf(stderr,"\n\nSOESA - Structure Optimization and Evaluation using Separations of Atoms\n");
  fprintf(stderr,"  Copyright 1999 by Michael Wall and Rice University\n");
  t=time(NULL);
  fprintf(stderr,"SOESA: Started at %s",asctime(localtime(&t)));

  /*
   * Read the amino acid info. file:
   */
  if (readaa) {
    aa = (AA_STRUCT *)malloc(sizeof(AA_STRUCT)*MAX_AA);
    for (j=0;j<MAX_AA;j++) {
      for (k=0;k<MAX_AA_NATOM;k++) {
	aa[j].atname[k] = (char *)calloc(TYPE_SIZE,sizeof(char));
      }
    }
    readAA(aa,&naa,&natype,aacfg_fname);
    /*
     * Display the amino acid info that was read in:
     */
    if (showaa) {
      for (i=0;i<naa;i++) {
	fprintf(stderr,"%c %s",aa[i].AAcode1,aa[i].AAcode3);
	for (j=0;j<aa[i].natom;j++) {
	  fprintf(stderr," %s:%3d",aa[i].atname[j],
		 getAtomidx(aa[i].AAcode3,aa[i].atname[j],aa,&naa));
	}
	fprintf(stderr,"\n");
      }
      fprintf(stderr,"Number of amino acids = %d\n",naa);
      fprintf(stderr,"Number of atom types = %d\n",natype);
    }
    
  }

  /*
   * Read the hash value file ("index" file):
   */
  if (readidx) {
    idx_fsize = getFileSize(idx_fname);
    idx_buf = (void *)malloc(idx_fsize);
    readFile(idx_fname,idx_buf,idx_fsize);
    idx=idx_buf;
  }

  /*
   * Read the PDF data file ("induced" file):
   */
  if (readind) {
    if (verbose) {
      fprintf(stderr,"SOESA: Reading PDF data file...\n");
    }
    ind_fsize = getFileSize(ind_fname); /* Determine file size */
    /* Allocate the buffer: */
    if ((ind_buf = (void *)malloc(ind_fsize)) == NULL) {
      printf("SOESA: Can't allocate %ld bytes for ind.\n",(long)ind_fsize);
      exit(1);
    }
    readFile(ind_fname,ind_buf,ind_fsize); /* Read the data */
    ind=ind_buf; /* ind is the data structure used to parse the data */
    if (verbose) {
      fprintf(stderr,"SOESA: Header:\nSOESA:   version=%d\nSOESA:   NLevels=%d\nSOESA:   NPoints=%d\nSOESA:   Cutloc=%g\nSOESA:   Cutter=%g\n",
	     ind->version,
	     ind->NLevels,
	     ind->NPoints,
	     ind->Cutloc,
	     ind->Cutter);
    }
    TLevel=ind->NLevels-1; /* 'm' greater or equal this is tertiary */
  }

  /*
   * Read the spline coefficients:
   */
  if (readspl) {
    if (verbose) {
      fprintf(stderr,"SOESA: Reading spline coefficients...\n");
    }
    spl_fsize = getFileSize(spl_fname); /* Determine the file size */
    /* Allocate the buffer: */
    if ((spl_buf = (void *)malloc(spl_fsize)) == NULL) {
      printf("SOESA: Can't allocate %ld bytes for spl.\n",(long)spl_fsize);
      exit(1);
    }
    readFile(spl_fname,spl_buf,spl_fsize); /* Read the data */
    spl=spl_buf; /* spl is the data structure used to parse the buffer */
    if (verbose) {
      fprintf(stderr,"SOESA: Header:\nSOESA:   version=%d\nSOESA:   NLevels=%d\nSOESA:   NPoints=%d\nSOESA:   Cutloc=%g\nSOESA:   Cutter=%g\n",
	     spl->version,
	     spl->NLevels,
	     spl->NPoints,
	     spl->Cutloc,
	     spl->Cutter);
    }
  }

  if (readindspl) {
    if (verbose) {
      fprintf(stderr,"SOESA: Reading PDF data file...\n");
    }
    indspl_fsize = getFileSize(indspl_fname); /* Determine the file size */
    /* Allocate the buffer: */
    if ((indspl_buf = (void *)malloc(indspl_fsize)) == NULL) {
      printf("SOESA: Can't allocate %ld bytes for indspl.\n",(long)indspl_fsize);
      exit(1);
    }
    readFile(indspl_fname,indspl_buf,indspl_fsize); /* Read the data */
    indspl=indspl_buf; /* spl is the data structure used to parse the buffer */
    if (verbose) {
      fprintf(stderr,"SOESA: Header:\nSOESA:   version=%d\nSOESA:   NLevels=%d\nSOESA:   NPoints=%d\nSOESA:   Cutloc=%g\nSOESA:   Cutter=%g\n",
	     indspl->version,
	     indspl->NLevels,
	     indspl->NPoints,
	     indspl->Cutloc,
	     indspl->Cutter);
    }
  }

  /*
   * Output the PDB file specified by -pdb:
   */
  if (showpdb) {
    /* Read the file: */
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    if ((fptr = fopen(out_fname, "w")) == NULL) {
      perror("\nCan't open output file.\n");
      exit(1);
    }
    for (i=0;i<natom;++i) {
      fprintf(fptr,"%5d %s %s %c %5d %4.3f %4.3f %4.3f %2.2f %2.2f %s\n", 
  	     i,type[i],resid[i],chain[i],resn[i],x[i],y[i],z[i], 
  	     q[i],b[i],segid[i]); 
    }
    fclose(fptr);
  }

  /*
   * Output the single-letter AA sequence of the PDB file specified by -pdb:
   */
  if (showseq) {
    /* Read the file: */
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    if ((fptr = fopen(out_fname, "w")) == NULL) {
      perror("\nCan't open output file.\n");
      exit(1);
    }
    resn_old = -1; /* Residue number of the last residue */
    resn_ct = 0; /* Total number of residues currently displayed */
    for (i=0;i<natom;++i) {
      if (resn[i] != resn_old) {
	++resn_ct;
	aaidx=getAA3idx(resid[i],aa,&naa);
	fprintf(fptr,"%c",aa[aaidx].AAcode1);
	if (resn_ct%70 == 0) { /* Display 70 AA's per line */
	  fprintf(fptr,"\n");
	}
	resn_old=resn[i];
      }
    }
    fclose(fptr);
  }

  /*
   * Output the PDF for the specified atoms + m:
   */
  if (showpdf) {
    showPDF(out_fname, &m, resid1, type1, resid2, type2, &zoompdf, indspl,
	    idx, aa, &naa, &natype, &natom);
  }

  /*
   * Enter into server mode, waiting for commands at the FIFO:
   */
  serverquit=0;
  if (servermode && pdb_fname[0] && out_fname[0]) {
    /* Allocate the msg buffer */
    msg_buf = (char *)calloc(MSG_SIZE,sizeof(char)); 
    parsed_msg = (MSG_STRUCT *)msg_buf; /* Use parsed_msg to parse msg_buf */
    mode = S_IFIFO|S_IRUSR|S_IWUSR; /* File mode for the FIFO */
    /* If a name for the fifo not given, use default */
    if (!fifo_fname[0]) { 
      strcpy(fifo_fname,FIFO_FILENAME);
    }
    mknod(fifo_fname,mode,NULL); /* Create the FIFO */
    /*
     * raw_msg is used to read the FIFO.  It's a structure that consists of
     *   a buffer, an integer specifying the number of bytes read into the
     *   buffer, and an integer specifying the max number of bytes that
     *   can be read into the buffer.
     */
    raw_msg.buf = (char *)msg_buf;
    raw_msg.maxlen = (int)MSG_SIZE;
    /* 
     * Print some log info:
     */
    fprintf(stderr,"SOESA: PDB file name is '%s'\n",pdb_fname);
    fprintf(stderr,"SOESA: Output file name is '%s'\n",out_fname);
    t=time(NULL);
    fprintf(stderr,"SOESA: %s",asctime(localtime(&t)));
    fprintf(stderr,"SOESA: Listening at FIFO with name '%s'...\n",fifo_fname);
    /*
     * Open the FIFO.  fifo_stream is the handle:
     */
    fifo_stream = open(fifo_fname,O_RDONLY);
    /*
     * Enter server mode forever and ever...
     */
    while (!serverquit) {
      servercmd[0]=0; /* There is no requested command */
      dummy=0; /* This parameter is required for getmsg, and should be zero */
      /*
       * This will listen at the FIFO, waiting until a message comes:
       */
      getmsg(fifo_stream,NULL,&raw_msg,&dummy);
      /* 
       * A message has been received.
       */
      if (debug) printf("len=%d\n",raw_msg.len); 
      /*
       * Print a time stamp:
       */
      t=time(NULL);
      fprintf(stderr,"SOESA: %s",asctime(localtime(&t)));
      /*
       * Copy the command string from the message into 'servercmd':
       */
      badcmd=0;
      if (parsed_msg->cmd == SHORT) {
	strncpy(servercmd,TNT_SHORT_CMD,strlen(TNT_SHORT_CMD));
	servercmd[strlen(TNT_SHORT_CMD)]=0;
      } else if (parsed_msg->cmd == LONG) {
	strncpy(servercmd,TNT_LONG_CMD,strlen(TNT_LONG_CMD));
	servercmd[strlen(TNT_LONG_CMD)]=0;
      } else if (parsed_msg->cmd == EVAL) {
	strncpy(servercmd,EVAL_CMD,strlen(EVAL_CMD));
	servercmd[strlen(EVAL_CMD)]=0;
      } else if (parsed_msg->cmd == XPLOR) {
	strncpy(servercmd,XPLOR_CMD,strlen(XPLOR_CMD));
	servercmd[strlen(XPLOR_CMD)]=0;
      } else if (parsed_msg->cmd == SHUTDOWN) {
	serverquit = 1;
	strncpy(servercmd,SHUTDOWN_CMD,strlen(SHUTDOWN_CMD));
	servercmd[strlen(SHUTDOWN_CMD)]=0;
      } else {
	perror("SOESA: Received invalid command from client.\n");
	servercmd[0]=0;
	badcmd=1;
      }
      /*
       * Print some log info:
       */
      if (!badcmd) {
	fprintf(stderr,"SOESA: Client with PID %d sent '%s' request.\n",
		(int)parsed_msg->pid,servercmd);
      }
      else {
	fprintf(stderr,
		"SOESA: Client with PID %d sent an invalid request.\n",
		(int)parsed_msg->pid);
      }

      /*
       * Point calc_home to the position for the path name:
       */
      calc_home = (char *)(parsed_msg + sizeof(pid_t) + sizeof(CMD));

      /*
       * Print some log info:
       */
      fprintf(stderr,"SOESA: Path name is '%s'\n",calc_home);  

      /*
       * Generate full path names for the pdb file and output file:
       *  11/7/99:Note that the names are currently provided at command line,
       *  not through the FIFO. - MW
       */
      path_len_bytes = strlen(calc_home);
      full_fname1 = (char *)calloc(path_len_bytes + 
				   strlen(pdb_fname) + 1,sizeof(char));
      strcpy(full_fname1,calc_home);
      strcpy(&full_fname1[path_len_bytes],pdb_fname);

      full_fname2 = (char *)calloc(path_len_bytes + 
				   strlen(out_fname) + 1,sizeof(char));
      strcpy(full_fname2,calc_home);
      strcpy(&full_fname2[path_len_bytes],out_fname);

      /*
       * Read the pdb file:
       */
      if (!readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,full_fname1)) {
	if (verbose) {
	  fprintf(stderr,"SOESA: Read %d atoms from input PDB file.\n",natom);
	}
	/*
	 * Perform the action requested:
	 */

	/*
	 * Calculate energy and display in TNT long format:
	 */
      if (parsed_msg->cmd == LONG) {
	calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
		   &natype,&natom,e,servercmd,msel);
	if (writeTNTlong(e,&escale,&natom,full_fname1,full_fname2)) {
	  fprintf(stderr, "SOESA: Error in writeTNTlong(), filename '%s'",
		  full_fname2);
	}
      } 

      /*
       * Calculate energy and display in TNT short format (just the energy:
       */
      else if (parsed_msg->cmd == SHORT) {
	calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
		   &natype,&natom,e,servercmd,msel);
	if (writeTNTshort(e,&escale,&natom,full_fname2)) {
	  fprintf(stderr, "SOESA: Error in writeTNTshort(), filename '%s'",
		  full_fname2);
	}
      } 

      /*
       * Calculate energy using dselect, rcut info and display in XPLOR format:
       */
      else if (parsed_msg->cmd == XPLOR) {
	/*
	 * Capture double selection info in dselect:
	 */
	dselect = (int *)calloc(natom,sizeof(int));
	for (i=0;i<natom;i++) {
	  dselect[i] = (int)atoi(segid[i]);
	}
	/*
	 * Calculate the energy:
	 */
	calcEnergyXPLOR(chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
			&natype,&natom,e,msel,dselect,&rcut); 
	/*
	 * Calculate the total energy:
	 */
	pdfe=0;
	for (i=0;i<natom;i++) {
	  pdfe += e[i].energy;
	}
	/*
	 * Output the energy, derivs:
	 */
	if ((fptr = fopen(full_fname2, "w")) == NULL) {
	  perror("SOESA: Can't open output file.\n");
	  exit(1);
	}
	fprintf(fptr,"%e\n",-pdfe*escale);
	for (i=0;i<natom;i++) {
	  fprintf(fptr,"%e %e %e\n",
		  -e[i].drvx*escale,
		  -e[i].drvy*escale,
		  -e[i].drvz*escale);
	}
	fclose(fptr);
	free((int *)dselect);
      }

      /*
       * Calculate the PDF score (energy) by residue:
       */
      else if (parsed_msg->cmd == EVAL) {
	/*
	 * Allocate arrays for scores eval, 
	 *  number of atom pairs used to calculate scores neval
	 */
	eval=(float *)calloc(natom,sizeof(float));
	neval=(int *)calloc(natom,sizeof(int));
	/*
	 * Calculate the energy:
	 */
	calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
		   &natype,&natom,e,"short",msel);
	/*
	 * Average the energy by residue, store results in eval:
	 */
	calcAvgByRes(resn,e,&natom,eval,neval);
	if ((fptr = fopen(full_fname2, "w")) == NULL) {
	  fprintf(stderr,"SOESA: Can't open output file '%s'\n",full_fname2);
	} else {
	  for (i=0;i<natom;i++) { 
	    if (neval[i]>0) { 
	      fprintf(fptr,"%d %f\n",i,-escale*eval[i]); 
	    } 
	  } 
	}
	free((float *)eval);
	free((int *)neval);
	fclose(fptr);
      }
      free((char *)full_fname2);
      } else {
	fprintf(stderr,"SOESA: Error in readPDB(), filename '%s'\n",
		full_fname1);
      }
      /*
       * Send a signal to the client that we're done:
       */
      if (kill(parsed_msg->pid,SIGHUP) < 0) {
	perror("SOESA: Error sending SIGHUP to client.\n");
      }
      if (serverquit) {
	close(fifo_stream);
	unlink(fifo_fname);
	fprintf(stderr,"SOESA: Shutting down.\n");
      }
      else {
	fprintf(stderr,"SOESA: Request '%s' has been carried out.\n",servercmd);
	t=time(NULL);
	fprintf(stderr,"SOESA: %s",asctime(localtime(&t)));
	fprintf(stderr,"SOESA: Listening...\n");
      }
      free((char *)full_fname1);
    }
  }

  /*
   * The following are energy calculations are requested directly by command
   *   line directive, rather than in server mode.  Lines are just copied
   *   from relevant sections of servermode, with file name substitutions.
   */

  /*
   * Calculate energy in TNT long format (with derivs):
   */
  if (TNTlong) {
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
	       &natype,&natom,e,"long",msel);
    if (writeTNTlong(e,&escale,&natom,pdb_fname,out_fname)) {
      fprintf(stderr, "SOESA: Error in writeTNTlong(), filename '%s'",
	      out_fname);
    }
  }

  /*
   * Calculate energy in TNT short format (just energy):
   */
  if (TNTshort) {
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
	       &natype,&natom,e,"short",msel);
    if (writeTNTshort(e,&escale,&natom,out_fname)) {
      fprintf(stderr, "SOESA: Error in writeTNTshort(), filename '%s'",
	      full_fname2);
    }
  }

  /*
   * Calculate energy in xplor format:
   */
  if (xplor) {
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    /*
     * Capture double selection info in dselect:
     */
    dselect = (int *)calloc(natom,sizeof(int));
    for (i=0;i<natom;i++) {
      dselect[i] = (int)atoi(segid[i]);
    }
    /*
     * Calculate the energy:
     */
    calcEnergyXPLOR(chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
		    &natype,&natom,e,msel,dselect,&rcut); 
    /*
     * Calculate the total energy:
     */
    pdfe=0;
    for (i=0;i<natom;i++) {
      pdfe += e[i].energy;
    }
    /*
     * Output the energy, derivs:
     */
    if ((fptr = fopen(out_fname, "w")) == NULL) {
      perror("SOESA: Can't open output file.\n");
      exit(1);
    }
    fprintf(fptr,"%e\n",-pdfe*escale);
    for (i=0;i<natom;i++) {
      fprintf(fptr,"%e %e %e\n",
	      -e[i].drvx*escale,
	      -e[i].drvy*escale,
	      -e[i].drvz*escale);
    }
    fclose(fptr);
    free((int *)dselect);
}

  /*
   * Calculate the PDF score (energy) by residue):
   */
  if (evalpdb) {
    readPDB(type,resid,chain,resn,x,y,z,q,b,segid,&natom,pdb_fname);
    /*
     * Allocate arrays for scores eval, 
     *  number of atom pairs used to calculate scores neval
     */
    eval=(float *)calloc(natom,sizeof(float));
    neval=(int *)calloc(natom,sizeof(int));
    /*
     * Calculate the energy:
     */
    calcEnergy(segid,chain,resn,resid,type,x,y,z,indspl,idx,aa,&naa,
	       &natype,&natom,e,"short",msel);
    /*
     * Average the energy by residue, store results in eval:
     */
    calcAvgByRes(resn,e,&natom,eval,neval);
    if ((fptr = fopen(out_fname, "w")) == NULL) {
      fprintf(stderr,"SOESA: Can't open output file '%s'\n",out_fname);
    } else {
      for (i=0;i<natom;i++) { 
	if (neval[i]>0) { 
	  fprintf(fptr,"%d %f\n",i,-escale*eval[i]); 
	} 
      } 
    }
    free((float *)eval);
    free((int *)neval);
    fclose(fptr);
  }

  /*
   * Reverse the byte order of the pdf data buffer:
   */
  if (revind) {
    xformRevByte(ind_buf,ind_fsize);
    fprintf(stderr,"Byte order reversed for IDX buffer.\n");
  }
  
  /*
   * Reverse the byte order of the pdf hash value buffer:
   */
  if (revidx) {
    xformRevByte(idx_buf,idx_fsize);
    fprintf(stderr,"Byte order reversed for IND buffer.\n");
  }

  /*
   * Calculate the spline coefficients for the pdf data:
   */
  if (calcspline) {
    xformCalcSpline(ind,ind_fsize);
  }

  /*
   * Merge the PDF data plus spline coefficients into a single file:
   */
  if (merge) {
    if ((fptr = fopen(indspl_fname, "w")) == NULL) {
      perror("\nSOESA: Can't open file for merged data.\n");
      exit(1);
    }
    /*
     * Define the size of a single PDF distn:
     */
    datasize = ind->NPoints*sizeof(float);
    /*
     * Define the size of a single PDF record:
     */
    pdfsize = datasize + PDF_HEADER_LENGTH;

    /*
     * Write the ind header to the indspl file:
     */
    num_wrote = fwrite(ind_buf,1,IND_HEADER_LENGTH,fptr);

    /*
     * Point the ind buffer to the beginning of the data:
     */
    ind_buf = (void *)&((char *)ind_buf)[IND_HEADER_LENGTH];

    /*
     * Point the spl buffer to the first spline coefficient array:
     */
    spl_buf = (void *)&((char *)spl_buf)[IND_HEADER_LENGTH + datasize];

    if (verbose) fprintf(stderr,"SOESA: Writing a merged ind, spl file...");
    while (i<ind_fsize) {
      /* Write the entire ind pdf: */
      num_wrote = fwrite(ind_buf,1,pdfsize,fptr);
      /* Write the spline coeffs: */
      num_wrote = fwrite(spl_buf,1,datasize,fptr);
      /* Point the ind, spl buffers to the next pdf record: */
      ind_buf = (void *)&((char *)ind_buf)[pdfsize];
      spl_buf = (void *)&((char *)spl_buf)[pdfsize];
      i += pdfsize;
    }
    /*
     * Update the hash table:
     */
    i=j=0;
    idx_offset = 0;
    while (i<idx_fsize) {
      if (idx[j] != -1) {
	idx[j] += idx_offset;
	idx_offset += datasize;
      }
      i += sizeof(long);
      j++;
    }
    if (verbose) fprintf(stderr,"done.\n");
    fclose(fptr);
  }

  /*
   * Output the pdf data buffer in a file:
   */
  if (writeind) {
    if (!out_fname[0]) {
      perror("No file name specified for writing IND.\n");
      exit(1);
    }
    writeFile(out_fname,ind_buf,ind_fsize);
  }

  /*
   * Output the pdf hash value buffer in a file:
   */
  if (writeidx) {
    if (!out_fname[0]) {
      perror("No file name specified for writing IDX.\n");
      exit(1);
    }
    writeFile(out_fname,idx_buf,idx_fsize);
  }

  if (readidx) free(idx_buf);
  if (readind) free(ind_buf);
  if (readspl) free(spl_buf);
  if (readaa) free(aa);
  t=time(NULL);
  fprintf(stderr,"SOESA: Stopped at %s",asctime(localtime(&t)));
}

