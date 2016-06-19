#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"

PetscBool firstTime = PETSC_TRUE;
char externalSignalScriptName[PETSC_MAX_PATH_LEN];  
PetscBool useExternalSignalScript=PETSC_FALSE;
PetscBool useSignalFiles=PETSC_FALSE;
char *signalfiles[2];
PetscInt signalWaitTime = -1;

#undef __FUNCT__
#define __FUNCT__ "waitForSignal"
PetscErrorCode waitForSignal(PetscInt waitTime)
{
  PetscErrorCode ierr;
  
  FILE *fd1;  
  PetscViewer fd;
  PetscInt iDone=0;
  PetscInt zero=0;
  PetscInt one=1;
  PetscInt fp;
  PetscInt numsignalfiles = 2;
  PetscBool flg;  
  
  if (firstTime) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing signaling module ...\n");CHKERRQ(ierr);                  
    ierr = PetscOptionsGetString(PETSC_NULL,"-signalscript",externalSignalScriptName,PETSC_MAX_PATH_LEN-1,&useExternalSignalScript);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL,"-signalfiles",&useSignalFiles);CHKERRQ(ierr);
    if ((useExternalSignalScript) && (useSignalFiles)) {
      SETERRQ(PETSC_COMM_WORLD,1,"Cannot specify both an external signal script and signal files!");
    }  
    if ((!useExternalSignalScript) && (!useSignalFiles)) {
      SETERRQ(PETSC_COMM_WORLD,1,"Must specify an external signal script OR signal files!");
    }  
    if (useExternalSignalScript) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"External signaling script has been specified: %s\n",externalSignalScriptName);CHKERRQ(ierr);      
    }
    if (useSignalFiles) {
      ierr = PetscOptionsGetStringArray(PETSC_NULL,"-signalfiles",signalfiles,&numsignalfiles,&flg);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Signal files have been specified: %s, %s\n",signalfiles[0],signalfiles[1]);CHKERRQ(ierr);                  
      ierr = PetscOptionsGetInt(PETSC_NULL,"-signalwaittime",&signalWaitTime,&flg);CHKERRQ(ierr);
      if (flg) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Signal wait time of %d seconds has been specified\n",signalWaitTime);CHKERRQ(ierr);                        
      } else {
        signalWaitTime = -1;
      }
    }    
    
    firstTime = PETSC_FALSE;
  } else {  
    if (useExternalSignalScript) {
      ierr = PetscPOpen(PETSC_COMM_WORLD,PETSC_NULL,externalSignalScriptName,"r",&fd1);CHKERRQ(ierr);  
      ierr = PetscPClose(PETSC_COMM_WORLD,fd1,PETSC_NULL);CHKERRQ(ierr);  
    } else {  
      
      if (signalWaitTime>0) waitTime = signalWaitTime; /* overwrite with runtime option */

/*    overwrite external signal file       */
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,signalfiles[1],FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinarySynchronizedWrite(PETSC_COMM_WORLD,fp,&zero,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/*    send "ready" signal   */
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,signalfiles[0],&fd);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(fd,"%d\n",one);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/*    wait for external signal   */
      while (iDone==0) {
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,signalfiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);    
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinarySynchronizedRead(PETSC_COMM_WORLD,fp,&iDone,1,PETSC_INT);CHKERRQ(ierr);    
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        PetscSleep(waitTime);      
      }

/*    send "busy" signal   */
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,signalfiles[0],&fd);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(fd,"%d\n",zero);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    }    
  }

  return 0;    
}
