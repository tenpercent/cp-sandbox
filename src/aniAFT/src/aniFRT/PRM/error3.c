#include <stdio.h>
#include <stdlib.h>
#include "error3.h"

void errorExit3 (int group, char *number) {
	switch (group) {
		case  GE_Memory:
			printf("\nOut  of  memory!\n");
			break;
		case  GE_Critical:
			printf("\nCritical  error  in  program!\n");
			break;
		case  GE_User:
			printf("\nError  in  function  of  user!\n");
			break;
		case  GE_Temp:
			printf("\nTempoparal  error!\n");
			break;
	}
	printf("\n\n\n%s\n",number);

	exit(1);
	return;
} /* errorExit3 */

