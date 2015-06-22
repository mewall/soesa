/* 

   PDFTOOL.C - A library of useful routines for manipulating PDFs
               Copyright 1999 by Michael Wall and Rice University
   
               Distribution and modification of this package is governed
	       under the "Artistic License" (perl license) described in the
	       file LICENSE in the SOESA distribution.

   Author: Michael Wall
   Date: 11/7/1999
   Version: 0.2
      Added more comments for public release
   Date: 2/11/1999
   Version: 0.1

   
*/

#include<pdf.h>

/*
 * The following routines are based on ones in the Numerical Recipes library:
 */

/*
 * Error printing:
 */
void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

/*
 * Allocate a vector:
 */
float *vector(nl,nh)
int nl,nh;
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

/*
 * Free a vector:
 */
void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

/*
 * Calculate spline coefficients:
 */
void spline(x,y,n,yp1,ypn,y2)
float x[],y[],yp1,ypn,y2[];
int n;
{
	int i,k;
	float p,qn,sig,un,*u,*vector();
	void free_vector();

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

/*
 * Interpolate the PDF using spline coefficients:
 */
void getPDFVal(float *xlo, float *xhi,float *ya,float *y2a,int *n,float *x,float *y)
{
	float h,b,a;
	void nrerror();

	h=*xhi-*xlo;
	if (h == 0.0) perror("Bad XA input to routine getPDFVal");
	a=(*xhi-*x)/h;
	b=(*x-*xlo)/h;
	*y=a*ya[*n]+b*ya[*n+1]+((a*a*a-a)*y2a[*n]+(b*b*b-b)*y2a[*n+1])*(h*h)/6.0;
}

/*
 * Calculate the derivative of the PDF using the spline coefficients:
 */
void getPDFDrv(float *xlo, float *xhi,float *ya, float *y2a,int *n,float *x,float *y)
{
	float h,b,a;
	void nrerror();

	h=*xhi-*xlo;
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(*xhi-*x)/h;
	b=(*x-*xlo)/h;
	*y=(ya[*n+1]-ya[*n])/h+(y2a[*n]*(3.0*b*b-1.0)-y2a[*n+1]*(3.0*a*a-1.0))*h/6.0;
}

/*
 * Calculate the second derivative of the PDF using the spline coefficients:
 */
void getPDFCrv(float *xlo, float *xhi,float *ya,float *y2a,int *n,float *x,float *y)
{
	float h,b,a;
	void nrerror();

	h=*xhi-*xlo;
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(*xhi-*x)/h;
	b=(*x-*xlo)/h;
	*y=(y2a[*n+1]*(*x-*xlo) + y2a[*n]*(*xhi-*x))/h;
}

/*
 * The following were wholly written by MW:
 */

/*
 * Read a PDB file:
 */
int readPDB(char **type, char **resid, char *chain, int *resn, 
	    float *x, float *y, float *z, float *q, float *b, 
	    char **segid, int *natom, char *fname)
{

  FILE *fptr;
  int 
    iatn,
    i=0,
    idig,
    ict=0;
  char *input_line=(char *)calloc(INPUT_LINE_SIZE+1,sizeof(char));
  char *fil=(char *)calloc(INPUT_LINE_SIZE+1,sizeof(char));
  char *typestr=(char *)calloc(TYPE_SIZE+1,sizeof(char));
  char *resnstr=(char *)calloc(RESN_SIZE+1,sizeof(char));
  char *token;

  if ((fptr = fopen(fname, "r")) != NULL) {

  typestr[TYPE_SIZE-1]=0;

  while (fgets(input_line,INPUT_LINE_SIZE,fptr) != NULL) {
    if (!strncmp(input_line, "ATOM", 4)) {
      sscanf(input_line, 
	     "%6c%5d%1c%4c%1c%4c%1c%4c%5c%8g%8g%8g%1c%5g%1c%5g%6c%4c%4c",
	     fil,&iatn,fil,typestr,fil,resid[i],&chain[i],resnstr,fil,
	     &x[i],&y[i],&z[i],
	     fil,&q[i],fil,&b[i],fil,segid[i],fil);
      if ((token=strtok(typestr," ")) == NULL) {
	perror("READPDB: Can't parse atom type.\n");
      }
      strcpy(type[i],token);
      if ((token=strtok(resnstr," ")) == NULL) {
	perror("READPDB: Can't parse residue number.\n");
      }
      idig=0;
      while ((idig<strlen(token))&&isdigit(token[idig])) idig++; 
      if (idig==strlen(token)) {
	resn[i]=atoi(resnstr);
	i++;
      }
    }
    else {
      ict++;
    }
  }
  *natom=i;
  fclose(fptr);
  return(0);
  } else {
    return(1);
  }


}

/*
 * Read an amino acid configuration file:
 */
int readAA(AA_STRUCT *aa, int *naa, int *natype, char *fname)
{
  FILE *fptr;

  int 
    i,
    j,
    k,
    aaidx=0;

char 
  input_line[INPUT_LINE_SIZE],
  *token;

  if ((fptr = fopen(fname, "r")) == NULL) {
    perror("\nREADPDB: Can't open amino acid configuration file.\n");
    exit(1);
  }

  i=-1;  
  while (fgets(input_line,INPUT_LINE_SIZE-1,fptr) != NULL) {
    if (input_line[strlen(input_line)-1]=='\n') {
      input_line[strlen(input_line)-1]=0;
    } 
    if (input_line[0] != '#') {
      if (i==-1) {
	sscanf(input_line,"%d %d",naa,natype);
	++i;
      }
      else {
	aa[i].baseidx=aaidx;
	if ((token = strtok(input_line, " ")) == NULL) {
	  perror("READAA: Couldn't read single-letter AA code.");
	  exit(1);
	}
	else {
	  aa[i].AAcode1=token[0];
	}
	if ((token = strtok(NULL, " ")) == NULL) {
	  perror("READAA: Couldn't read three-letter AA code.");
	  exit(1);
	}
	else {
	  strncpy(aa[i].AAcode3,token,RESID_SIZE);
	}
	aa[i].natom=0;
	while ((token = strtok(NULL, " ")) != NULL) {
	  strcpy(aa[i].atname[aa[i].natom], token);
	  ++aa[i].natom;
	}
	aaidx += aa[i].natom;
	++i;
      }
    }
  }
  if (i != *naa) {
    perror("READAA: AA count doesn't match that in file header.\n");
    exit(1);
  }
  if (aaidx != *natype) {
    perror("READAA: Number of atom types doesn't match that in file header.\n");
    exit(1);
  }
  return(0);
}

/*
 * Get the amino acid index based on a 3-letter AA type:
 */
int getAA3idx(char *AAcode3,AA_STRUCT *aa,int *naa) 
{
  int i;

  i=0;
  while ((i<*naa) && strncmp(AAcode3,aa[i].AAcode3,strlen(aa[i].AAcode3))) ++i;

  if (i==*naa) {
    return(-1);
  }
  else {
    return(i);
  }
}

/*
 * Get the atom index, accounting for both AA type and chemical type:
 */
int getAtomidx(char *AAcode3,char *atname,AA_STRUCT *aa,int *naa) 
{
  int 
    i,
    aaidx;

  i=0;
  
  aaidx=getAA3idx(AAcode3,aa,naa);

  if (aaidx >= 0) {
    while ((i<aa[aaidx].natom) && strcmp(atname,aa[aaidx].atname[i])) ++i;
  }

  if ((aaidx < 0)) {
    return(-1);
  }
  else if (i==aa[aaidx].natom) {
    return(-1);
  }
  else {
    return(aa[aaidx].baseidx+i);
  }
}

/*
 * Get the hash value associated with a PDF distribution:
 */
long getPDFhash(int *m,char *resid1,char *type1,char *resid2,char *type2,
		AA_STRUCT *aa,int *naa,int *natype, int *TLevel)
{

  int
    tmp,
    atomidx1,
    atomidx2;

  long 
    PDFhash;

  atomidx1 = getAtomidx(resid1,type1,aa,naa);
  atomidx2 = getAtomidx(resid2,type2,aa,naa);

  if ((*m<0) || (atomidx1<0) || (atomidx2<0)) {
    PDFhash=-1;
  }
  else if (*m==0) {
    if (atomidx1>atomidx2) {
      tmp=atomidx1;
      atomidx1=atomidx2;
      atomidx2=tmp;
    }
    if (atomidx1==atomidx2) {
      PDFhash=-1;
    }
    else {
      PDFhash=(atomidx1*(*natype) + atomidx2);
    }
  }
  else if (*m>0) {
    if (*m<*TLevel) {
      PDFhash=(*m*(*natype)*(*natype) + atomidx1*(*natype) + atomidx2);
    }
    else {
      PDFhash=(*TLevel*(*natype)*(*natype) + atomidx1*(*natype) + atomidx2);
    }
  }
  return(PDFhash);
}

/*
 * Get the file size of a file from the file system:
 */
off_t getFileSize(const char *path)
{

  struct stat buf;

  stat(path,&buf);

  return(buf.st_size);

}

/*
 * Read a file into a buffer:
 */
int readFile(const char *fname, char *buf, off_t fsize)
{

  FILE *fptr;
  off_t num_read;

  if ((fptr = fopen(fname, "r")) == NULL) {
    perror("\nGETFILE: Can't open file.\n");
    exit(1);
  }

  num_read = fread(buf,1,fsize,fptr);

  if (num_read != fsize) {
    perror("\nGETFILE: Number of bytes read different from length.\n");
    exit(1);
  }
  fclose(fptr);
  return(0);
}

/*
 * Write a buffer to a file:
 */
int writeFile(const char *fname, char *buf, off_t fsize)
{

  FILE *fptr;
  off_t num_wrote;

  if ((fptr = fopen(fname, "w")) == NULL) {
    perror("\nWRITEFILE: Can't open file.\n");
    exit(1);
  }

  num_wrote = fwrite(buf,1,fsize,fptr);

  if (num_wrote != fsize) {
    perror("\nPUTFILE: Number of bytes wrote different from buffer length.\n");
    exit(1);
  }
  fclose(fptr);
  return(0);
}

/*
 * Write an energy, deriv file in TNT format:
 */
int writeTNTlong(ENERGY_STRUCT *e, float *escale, int *natom,
		 char *pdb_fname, char *out_fname)
{
  FILE *pdbfptr, *outfptr;
  int 
    iatn,
    i=0,
    j=0,
    idig,
    ict=0;
  char input_line[INPUT_LINE_SIZE];
  char fil[INPUT_LINE_SIZE];
  char typestr[TYPE_SIZE];
  char resnstr[RESN_SIZE];
  char chain;
  char *token;
  float pdfe;
  float x,y,z,q,b;

  if ((pdbfptr = fopen(pdb_fname, "r")) == NULL) {
    perror("\nWRITETNTLONG: Can't open PDB file.\n");
    exit(1);
  }
  if ((outfptr = fopen(out_fname, "w")) == NULL) {
    perror("\nWRITETNTLONG: Can't open output file.\n");
    exit(1);
  }

  pdfe=0;

  for (i=0;i<*natom;i++) {
    pdfe += e[i].energy;
  }

  fprintf(outfptr,"FUNCVAL PDFE %e\n",-pdfe*(*escale));

  i=0;

  while (fgets(input_line,INPUT_LINE_SIZE,pdbfptr) != NULL) {
    if (!strncmp(input_line, "ATOM", 4)) {
      for (j=0;j<TYPE_SIZE;j++) {
	typestr[j]=resnstr[j]=0;
      }
      sscanf(input_line, 
	     "%6c%5c%1c%4c%1c%4c%1c%4c%5c%8g%8g%8g%1c%5g%1c%5g%6c%4c%4c",
	     fil,fil,fil,typestr,fil,fil,&chain,resnstr,fil,
	     &x,&y,&z,
	     fil,&q,fil,&b,fil,fil,fil);
      if ((token=strtok(typestr," ")) == NULL) {
	perror("READPDB: Can't parse atom type.\n");
      }
      if ((token=strtok(resnstr," ")) == NULL) {
	perror("WRITETNTLONG: Can't parse residue number.\n");
      }
      idig=0;
      while ((idig<strlen(token))&&isdigit(token[idig])) idig++; 
      if (idig==strlen(token)) {
	fprintf(outfptr,"DRVC PDFE %+e %+e %+e %e %e %s %s %c\n",
		-e[i].drvx*(*escale),
		-e[i].drvy*(*escale),
		-e[i].drvz*(*escale),
		ZERO,ZERO,
		typestr,token,chain);
/*  	fprintf(outfptr,"CRVC PDFE %+e %+e %+e %e %e %s %s %c\n", */
/*  		-e[i].ddrvx*(*escale), */
/*  		-e[i].ddrvy*(*escale), */
/*  		-e[i].ddrvz*(*escale), */
/*  		ZERO,ZERO, */
/*  		typestr,token,chain); */
	i++;
      }
      else {
	fprintf(outfptr,"DRVC PDFE %+e %+e %+e %e %e %s %s %c\n",
		ZERO, ZERO, ZERO, ZERO, ZERO,
		typestr,token,chain);
/*  	fprintf(outfptr,"CRVC PDFE %+e %+e %+e %e %e %s %s %c\n", */
/*  		ZERO, ZERO, ZERO, ZERO, ZERO, */
/*  		typestr,token,chain); */
      }
    }
    else {
      ict++;
    }
  }
  if (*natom!=i) {
    perror("WRITETNTLONG: Atom count has changed.\n");
  }
  fclose(pdbfptr);
  fclose(outfptr);
  return(0);
}

/*
 * Write the energy in TNT short format:
 */
int writeTNTshort(ENERGY_STRUCT *e, float *escale, int *natom, char *out_fname)
{

  FILE *outfptr;

  int 
    i=0;

  float pdfe;

  if ((outfptr = fopen(out_fname, "w")) == NULL) {
    perror("\nWRITETNTSHORT: Can't open output file.\n");
    exit(1);
  }

  pdfe=0;
  for (i=0;i<*natom;i++) {
    pdfe += e[i].energy;
  }

  fprintf(outfptr,"FUNCVAL PDFE %e\n",-pdfe*(*escale));

  fclose(outfptr);
  return(0);
}

/*
 * Reverse the byte order of a buffer, assuming 4-byte values:
 */
int xformRevByte(char *buf, off_t bufsize) 
{
  size_t i;
  char tmp;

  if ((bufsize % 4) != 0) {
    perror("REVBYTE: Size of buffer not divisible by 4.\n");
    exit(1);
  }

  for (i=0;i<bufsize;i=i+4) {
    tmp=buf[i];
    buf[i]=buf[i+3];
    buf[i+3]=tmp;
    tmp=buf[i+1];
    buf[i+1]=buf[i+2];
    buf[i+2]=tmp;
  }
  return(0);
}

/*
 * Calculate the spline coefficients for the PDF data:
 */
int xformCalcSpline(IND_STRUCT *ind, off_t ind_fsize)
{

  off_t 
    pdfsize,
    pdfskip,
    i,
    j;

  PDF_STRUCT
    *pdf;

  float
    *x,
    *y2,
    *pdfpdf;

  pdfsize=ind->NPoints*sizeof(float);
  pdfskip=PDF_HEADER_LENGTH+pdfsize;

  x = (float *)malloc(ind->NPoints*sizeof(float));
  y2 = (float *)malloc(ind->NPoints*sizeof(float));

  i=IND_HEADER_LENGTH;

  while (i<ind_fsize) {
    pdf=(PDF_STRUCT *)&((char *)ind)[i];
    pdfpdf=(float *)&((char *)pdf)[PDF_HEADER_LENGTH];
    for(j=0;j<ind->NPoints;j++) {
      x[j]=pdf->x_start+j*pdf->x_step;
    }
    spline(x-1,pdfpdf-1,ind->NPoints,0.0,0.0,y2-1);
    bcopy(y2,pdfpdf,pdfsize);
    i += pdfskip;    
  }
  printf("\n");
  free((float *)x);
  free((float *)y2);
  return(0);
}

/*
 * Output a PDF distn:
 */
int showPDF(char *fname, 
	    int *m, char *resid1, char *type1, char *resid2, char *type2,
	    char *zoom, IND_STRUCT *indspl, int *idx, 
	    AA_STRUCT *aa, int *naa, int *natype, int *natom)
{

  FILE *fptr;

  long 
    PDFhash,
    datasize;

  float
    xlo,
    xhi,
    d,
    y,
    showpdf_step,
    *pdf,
    *spl;

  PDF_STRUCT
    *pdfspl;

  int
    nshow,
    i,
    didx,
    TLevel;

  TLevel=indspl->NLevels-1;

  datasize = indspl->NPoints*sizeof(float);

  PDFhash = getPDFhash(m,resid1,type1,resid2,type2,aa,naa,natype,&TLevel);
  if (PDFhash == -1) {
    perror("SHOWPDF: No distribution found for specified atom pairs.\n");
    exit(1);
  }
  pdfspl=(PDF_STRUCT *)&((char *)indspl)[idx[PDFhash]];
  pdf=(float *)&(((char *)pdfspl)[PDF_HEADER_LENGTH]);
  spl=(float *)&(((char *)pdfspl)[PDF_HEADER_LENGTH+datasize]);
  nshow=(indspl->NPoints-1)*(*zoom);
  showpdf_step=pdfspl->x_step/(float)*zoom;
  if ((fptr = fopen(fname, "w")) == NULL) {
    perror("\nCan't open output file.\n");
    exit(1);
  }
  for (i=0;i<nshow;i++) {
    d=pdfspl->x_start+showpdf_step*i;
    didx=(d-pdfspl->x_start)/pdfspl->x_step;
    xlo=pdfspl->x_start+pdfspl->x_step*didx;
    xhi=xlo+pdfspl->x_step;
    getPDFVal(&xlo,&xhi,pdf,spl,&didx,&d,&y);
    fprintf(fptr,"%g %g",d,y);
    getPDFDrv(&xlo,&xhi,pdf,spl,&didx,&d,&y);
    fprintf(fptr," %g",y);
    getPDFCrv(&xlo,&xhi,pdf,spl,&didx,&d,&y);
    fprintf(fptr," %g\n",y);
  }
  fclose(fptr);
  return(0);
}

/*
 * Calculate the PDF energy for TNT, eval:
 */
int calcEnergy(char **segid,char *chain,int *resn, char **resid, char **type,
	       float *x, float *y, float *z, 
	       IND_STRUCT *indspl, int *idx,
	       AA_STRUCT *aa, int *naa, int *natype, int *natom,
	       ENERGY_STRUCT *e, char *cmd, char *msel)
{

  char
    long_calc=0;

  int 
    didx,
    i,
    j,
    m;

  int
    ishow=0,
    TLevel,
    *nenergy;

  float 
    d,
    dx,
    dy,
    dz,
    xlo,
    xhi,
    PDFVal,
    PDFDrv,
    PDFCrv,
    *pdf,
    *spl;

  long
    PDFhash,
    datasize;

  PDF_STRUCT
    *pdfspl;

  if (!strncmp("long",cmd,4)) {
    long_calc=1;
  }
  else if (!strncmp("short",cmd,5)) {
    long_calc=0;
  }
  else {
    perror("CALCENERGY: Must supply command 'long' or 'short'.\n");
    exit(1);
  }

  TLevel=indspl->NLevels-1;
  datasize = indspl->NPoints*sizeof(float);

  for (i=0;i<*natom;i++) {
    e[i].energy=0;
    e[i].drvx=0;
    e[i].drvy=0;
    e[i].drvz=0;
    e[i].ddrvx=0;
    e[i].ddrvy=0;
    e[i].ddrvz=0;
    e[i].n=0;
  }

  for (i=0;i<*natom-1;i++) {
    for (j=i+1;j<*natom;j++) {
      if (resn[j] >= resn[i]) {
	m=resn[j]-resn[i];
	if (m>=TLevel || strcmp(segid[i],segid[j]) || (chain[i] != chain[j])) {
	  m=TLevel;
	}
	if (msel[m]) {
	  PDFhash=getPDFhash(&m,resid[i],type[i],resid[j],type[j],
			     aa,naa,natype,&TLevel);
	  if (PDFhash != -1) {
	    pdfspl=(PDF_STRUCT *)&((char *)indspl)[idx[PDFhash]];
	    pdf=(float *)&(((char *)pdfspl)[PDF_HEADER_LENGTH]);
	    spl=(float *)&(((char *)pdfspl)[PDF_HEADER_LENGTH+datasize]);
	    dx=x[i]-x[j];
	    dy=y[i]-y[j];
	    dz=z[i]-z[j];
	    d=sqrtf(dx*dx+dy*dy+dz*dz);
	    didx=(d-pdfspl->x_start)/pdfspl->x_step;
	    if ((didx>=0) && (didx<(indspl->NPoints-1))) {
	      xlo=pdfspl->x_start+pdfspl->x_step*didx;
	      xhi=xlo+pdfspl->x_step;
	      getPDFVal(&xlo,&xhi,pdf,spl,&didx,&d,&PDFVal);
	      e[i].energy += PDFVal;
	      e[j].energy += PDFVal;
	      e[i].n++;
	      e[j].n++;
	      if (long_calc) {
		getPDFDrv(&xlo,&xhi,pdf,spl,&didx,&d,&PDFDrv);
		e[i].drvx += PDFDrv*dx/d;
		e[j].drvx -= PDFDrv*dx/d;
		e[i].drvy += PDFDrv*dy/d;
		e[j].drvy -= PDFDrv*dy/d;
		e[i].drvz += PDFDrv*dz/d;
		e[j].drvz -= PDFDrv*dz/d;
	      }
	    }
	  }
	}
      }
    }
  }
  return(0);
}

/*
 * Calculate the PDF energy for XPLOR:
 */
int calcEnergyXPLOR(char *chain,int *resn, char **resid, char **type,
	       float *x, float *y, float *z, 
	       IND_STRUCT *indspl, int *idx,
	       AA_STRUCT *aa, int *naa, int *natype, int *natom,
	       ENERGY_STRUCT *e, char *msel,
		     int *dselect, float *rcut)
{

  int 
    didx,
    i,
    j,
    m;

  int
    ishow=0,
    TLevel,
    *nenergy;

  float 
    d,
    dx,
    dy,
    dz,
    xlo,
    xhi,
    PDFVal,
    PDFDrv,
    *pdf,
    *spl;

  size_t 
    PDFhash,
    datasize;

  PDF_STRUCT
    *pdfspl;


  TLevel=indspl->NLevels-1;
  datasize = indspl->NPoints*sizeof(float);

  for (i=0;i<*natom;i++) {
    e[i].energy=0;
    e[i].drvx=0;
    e[i].drvy=0;
    e[i].drvz=0;
    e[i].n=0;
  }
  for (i=0;i<*natom-1;i++) {
    for (j=i+1;j<*natom;j++) {
      if (((dselect[i]+dselect[j])<=1) && (resn[j] >= resn[i])) {
	m=resn[j]-resn[i];
	if (m>=TLevel || (chain[i] != chain[j])) {
	  m=TLevel;
	}
	if (msel[m]) {
	  PDFhash=getPDFhash(&m,resid[i],type[i],resid[j],type[j],
			     aa,naa,natype,&TLevel);
	  if (PDFhash != -1) {
	    pdfspl=(PDF_STRUCT *)&((char *)indspl)[idx[PDFhash]];
	    pdf=(float *)&((char *)pdfspl)[PDF_HEADER_LENGTH];
	    spl=(float *)&((char *)pdfspl)[PDF_HEADER_LENGTH+datasize];
	    dx=x[i]-x[j];
	    dy=y[i]-y[j];
	    dz=z[i]-z[j];
	    d=sqrtf(dx*dx+dy*dy+dz*dz);
	    didx=(d-pdfspl->x_start)/pdfspl->x_step;
	    if ((didx>=0) && (didx<(indspl->NPoints-1)) && (d<*rcut)) {
	      xlo=pdfspl->x_start+pdfspl->x_step*didx;
	      xhi=xlo+pdfspl->x_step;
	      getPDFVal(&xlo,&xhi,pdf,spl,&didx,&d,&PDFVal);
	      e[i].energy += PDFVal;
	      e[j].energy += PDFVal;
	      e[i].n++;
	      e[j].n++;
	      getPDFDrv(&xlo,&xhi,pdf,spl,&didx,&d,&PDFDrv);
	      e[i].drvx += PDFDrv*dx/d;
	      e[j].drvx -= PDFDrv*dx/d;
	      e[i].drvy += PDFDrv*dy/d;
	      e[j].drvy -= PDFDrv*dy/d;
	      e[i].drvz += PDFDrv*dz/d;
	      e[j].drvz -= PDFDrv*dz/d;
	    }
	  }
	}
      }
    }
  }
  return(0);
}
	
/*
 * Calculate the average energy by residue:
 */
int calcAvgByRes(int *resn, ENERGY_STRUCT *e, int *natom,
		 float *eval, int *neval)
{

  int i;

  for (i=0;i<*natom;i++) {
    if (e[i].n>0) {
      eval[resn[i]] += e[i].energy;
    }
    neval[resn[i]]+=e[i].n;
  }

  for (i=0;i<*natom;i++) {
    if (neval[i]>0) {
      eval[i] /= neval[i];
    }
  }
  return(0);
}
