/*
  CALCUSER.C - Signal handling routine to perform calculations in SOESA.

            Copyright 1999 by Michael Wall and Rice University

	    Distribution and modification of this package is governed
	    under the "Artistic License" (perl license) described in the
	    file LICENSE in the SOESA distribution.

  Author: Michael Wall
  Version: 0.2  12/23/99
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#define OK                0
#define ERROR             (-1)

int calcuser_(void)
{
	int     dummy;

	if (!pcreatel("user.csh",0)) {
	  printf(" CALCUSER: Couldn't execute energy calculator script.\n");
	  exit(ERROR);
	} else {

	  printf(" CALCUSER: Executing energy calculator...\n");

	  /*
	   * Wait for completion of child process:
	   */

	  wait(&dummy);

	  printf(" CALCUSER: ...done\n");

	}

        return( OK );
}




