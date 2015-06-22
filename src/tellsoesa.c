/*
  TELLSOESA.C - A sample client for the SOESA server mode
                Copyright 1999 by Michael Wall and Rice University

		Distribution and modification of this package is governed
		under the "Artistic License" (perl license) described in the
		file LICENSE in the SOESA distribution.

  Author: Michael Wall
  Version: 0.2  SOESA  11/8/1999
  Version: 0.1  OVID  5/19/1999
*/

#include <pdf.h>

/*
 * Do the following this client receives a signal:
 */
void catchSIGHUP()
{

  signal(SIGHUP,SIG_DFL);
  
  exit(0);
  
}

int main( int argc, char *argv[] )
{

  /*
   * Used in parsing command line arguments:
   */

  char *cmd; /* Points to a command from the argument list */
  int inpct; /* Keeps track of how many arguments are passed */
	     /* to a command line directive */

  MSG_STRUCT
    *parsed_msg;
  
  char
    *fifo_fname = (char *)calloc(FNAME_SIZE,sizeof(char)), /* FIFO name */
    *pdf_home,
    *servercmd = (char *)calloc(INPUT_LINE_SIZE,sizeof(char)),
    *calc_home,
    *full_fname;

  int
    dummy=0,
    fifo_stream,
    len;
  struct strbuf
    raw_msg;

  char
    *msg_buf = (char *)calloc(MSG_SIZE,sizeof(char));


  /*
   * Point the message data structure to the buffer:
   */
  parsed_msg = (MSG_STRUCT *)msg_buf;

  /*
   * The client file name comes after the fixed length data in the message:
   */
  calc_home = (char *)(parsed_msg + sizeof(pid_t) + sizeof(CMD));

 inpct=0;
 --argc; 
  while (argc>0) {
    if (argv[argc][0] != '-') {
      ++inpct;
    }
    else {
      cmd = argv[argc]+1;
      /*
       * FIFO name:
       */
      if (!strcmp(cmd,"fifo")) {
	if (inpct == 1) {
	  strcpy(fifo_fname,argv[argc+1]);
	}
	else {
	  perror("TELLSOESA: -fifo takes 1 argument only.");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Path where the PDB file is and where the output file is written:
       */
      if (!strcmp(cmd,"home")) {
	if (inpct == 1) {
	  strcpy(calc_home,argv[argc+1]);
	}
	else {
	  perror("TELLSOESA: -home takes 1 argument only.");
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Command request:
       */
      else if (!strcmp(cmd,"cmd")) {
	if (inpct == 1) {
	  strcpy(servercmd,argv[argc+1]);
	}
	else {
	  perror("TELLSOESA: -cmd takes 1 argument only.");
	  exit(1);
	}
	inpct=0;
      }
    }
    --argc;
  }

  /*
   * Build the message to be passed to the FIFO:
   */

  /*
   * If the command is defined, place it in the message to the FIFO:
   */
  if( !strcmp( servercmd, "short" ) ) {
    parsed_msg->cmd = SHORT;
  } else if ( !strcmp( servercmd, "long") ) {
    parsed_msg->cmd = LONG;
  } else if ( !strcmp(servercmd, "xplor")) {
    parsed_msg->cmd = XPLOR;
  } else if ( !strcmp(servercmd,"eval")) {
    parsed_msg->cmd = EVAL;
  } else if ( !strcmp(servercmd,"shutdown")) {
    parsed_msg->cmd = SHUTDOWN;
  } else {
    fprintf( stderr, "TELLSOESA: Command %s is not defined.\n", servercmd );
    return( -1 );
  }

  /*
   * Add in the PID for this process, so that we can receive a signal
   *    telling us that the request is done:
   */

  parsed_msg->pid = getpid();

  /*
   * Set up the data structure for message passing:
   */
  raw_msg.buf = (char *)msg_buf;
  raw_msg.len = MSG_SIZE;
  msg_buf[raw_msg.len]=0;

  /*
   * Build the full file name of the FIFO, starting with any path info
   *   in the shell environment:
   */
  if ((pdf_home = getenv("PDF_HOME")) != NULL) {
    len = strlen(pdf_home);
  }
  else len=0;

  /*
   * Use the default fifo name if it wasn't given:
   */
  if (!fifo_fname[0]) strcpy(fifo_fname,FIFO_FILENAME);

  /*
   * Allocate the space for the full filename string:
   */
  full_fname = (char *)calloc(len + strlen(fifo_fname) + 1,sizeof(char));

  /*
   * If the path is supplied, make sure there's a trailing slash:
   */
  if (pdf_home != NULL) {
    strcpy(full_fname,pdf_home);
    if (full_fname[len-1] != '/') {
      full_fname[len++]='/';
    }
  }

  strcpy(&full_fname[len],fifo_fname);

/*  printf("PID file = %s\n",full_fname); */

  /*
   * If the client path isn't defined, get it from the working dir:
   */
  if (!calc_home[0]) {
    if (!strcpy(calc_home, getenv("PWD"))) {
      perror("TELLSOESA: Unable to get client path via PWD\n");
      perror("               Please supply via -home directive.");
      exit(1);
    }
  }
  
  /*
   * Make sure there's a trailing slash in the client path:
   */
  len = strlen(calc_home);
  if (calc_home[len-1] != '/') {
    calc_home[len++]='/';
  }

  /*
   * Open the FIFO for writing:
   */
  if ((fifo_stream = open(full_fname,O_WRONLY))==-1) {
    fprintf(stderr,"TELLSOESA: Error opening '%s'\n",full_fname);
    perror("    Note: PDF_HOME may be defined incorrectly.\n");
    exit(1);
  }

/*  printf("About to get PWD:\n"); */

  /*
   * Send the message and wait for a completion signal:
   */
  if (!putmsg(fifo_stream,NULL,&raw_msg,dummy)) {

    /*
     * Wait until SIGHUP received, then goto catchSIGHUP()
     */
    signal(SIGHUP,catchSIGHUP);

    pause();

  } else {
    perror("TELLSOESA: Error writing to the FIFO.\n");
    exit(1);
  }
}



