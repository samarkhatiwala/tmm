#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// #include "mpi.h"
#include "petsc.h"
#include "tmm_signal_utils.h"

PetscBool firstTime = PETSC_TRUE;
MPI_Comm client;
MPI_Status status;
char port_name[MPI_MAX_PORT_NAME];
PetscMPIInt activeRank = 0;

#undef __FUNCT__
#define __FUNCT__ "getExternalSignal"
PetscInt getExternalSignal(const char cmd[])
{
// Run the external program 'cmd' and return the integer value output by cmd.
// NOTE: It is assumed that cmd returns a single integer number.

  PetscErrorCode ierr;
  FILE *fd;
  PetscInt i = -1;
  char output[PETSC_MAX_PATH_LEN];
  int status;
  PetscMPIInt rank;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
//   ierr = PetscPOpen(PETSC_COMM_WORLD,PETSC_NULL,cmd,"r",&fd);CHKERRQ(ierr);  
  if (rank == 0) {
    fd = popen(cmd,"r");
    if (fd == NULL) {
      SETERRQ(PETSC_COMM_WORLD,1,"Problem with popen");
    }
	if (fgets(output, PETSC_MAX_PATH_LEN, fd) != NULL) {
	  i=atoi(output);
	}
    status = pclose(fd);
	if (status == -1) {
      SETERRQ(PETSC_COMM_WORLD,1,"Problem with pclose");	
    }
  }
//   ierr = PetscPClose(PETSC_COMM_WORLD,fd);CHKERRQ(ierr);
  MPI_Barrier(PETSC_COMM_WORLD);
  ierr = MPI_Bcast(&i,1,MPI_INT,0,PETSC_COMM_WORLD);CHKERRMPI(ierr);
  MPI_Barrier(PETSC_COMM_WORLD);
  
  return i;
}

#undef __FUNCT__
#define __FUNCT__ "getExternalSignalMPI"
PetscInt getExternalSignalMPI(PetscInt isend)
{
// Run the external program 'cmd' and return the integer value output by cmd.
// NOTE: It is assumed that cmd returns a single integer number.

  PetscErrorCode ierr;
  FILE *fd;
  PetscInt irecv = -1;
  char output[PETSC_MAX_PATH_LEN];
  MPI_Status status;
  PetscMPIInt rank;

//   MPI_Init(&argc, &argv);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);
//   if (size != 1) {
//     fprintf(stderr, "Server too big");
//     exit(EXIT_FAILURE);
//   }

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank==0) {
	if (firstTime) {
	  ierr = PetscPrintf(PETSC_COMM_SELF,"Initializing signaling module ...\n");CHKERRQ(ierr);                    
	  MPI_Open_port(MPI_INFO_NULL, port_name);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Server available at port: %s\n",port_name);CHKERRQ(ierr);                  
	  MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, PETSC_COMM_SELF, &client);
	  ierr = PetscPrintf(PETSC_COMM_SELF,"Connection accepted\n");CHKERRQ(ierr);                    
	  
	  firstTime = PETSC_FALSE;
	}
    MPI_Send(&isend, 1, MPI_INT, 0, MPI_ANY_TAG, client);
    
	MPI_Recv(&irecv, 1, MPI_INT, 0, MPI_ANY_TAG, client, &status);
  
	if (irecv==0) {
	  MPI_Comm_free(&client);
	  MPI_Close_port(port_name); 
	}
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  ierr = MPI_Bcast(&irecv,1,MPI_INT,0,PETSC_COMM_WORLD);CHKERRMPI(ierr);
    
  return irecv;
}

// #undef __FUNCT__
// #define __FUNCT__ "getExternalSignalFile"
// PetscInt getExternalSignalFile(PetscInt isend)
// {
// 
//   PetscErrorCode ierr;
//   PetscInt irecv = -1;
//   PetscMPIInt rank;
//   FILE * fp;
// 
//   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
// 
//   if (rank==0) {
// 	if (isend==1) {
//   // ready to run
// 	  fp = fopen("readytorun","w");
// 	  fputs("1",fp);
//   //   echo $signal > readytorun
// 	  fclose(fp);
// 	} else if (isend==2) {
//   // finished running
// 	  fp = fopen("donerunning","w");
// 	  fputs("2",fp);
//   //   echo $signal > donerunning
// 	  fclose(fp);
// 	} else {
// 	  fp = fopen("problemwithrun","w");
// 	  fprintf(fp,"%d",isend);
//   //   echo $signal > problemwithrun
// 	  fclose(fp);
// 	}
// 
// 	while (access("runmodel", F_OK) != 0) {
// 	  sleep(30);
// 	}
// 
// 	fp = fopen("runmodel","r");
// 	fscanf(fp,"%d",&irecv);
// 	fclose(fp);
// 
// 	system("rm -f runmodel");
// 
//   }
// 
//   MPI_Barrier(PETSC_COMM_WORLD);
//   ierr = MPI_Bcast(&irecv,1,MPI_INT,0,PETSC_COMM_WORLD);CHKERRMPI(ierr);
//   
//   return irecv;
// }

#undef __FUNCT__
#define __FUNCT__ "getExternalSignalFile"
PetscInt getExternalSignalFile(PetscInt isend, PetscInt waitTime)
{
  PetscErrorCode ierr;
  PetscInt irecv = -1;
  PetscMPIInt rank, numProcessors;
  FILE * fp;

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank==activeRank) {
	if (isend==1) {
  // ready to run
	  fp = fopen("readytorun","w");
	  fputs("1",fp);
  //   echo $signal > readytorun
	  fclose(fp);
	while (access("runmodel", F_OK) != 0) {
	  sleep(waitTime);
	}

	fp = fopen("runmodel","r");
	fscanf(fp,"%d",&irecv);
	fclose(fp);

	system("rm -f runmodel");
	  
	} else if (isend==2) {
  // finished running
	  fp = fopen("donerunning","w");
	  fputs("2",fp);
  //   echo $signal > donerunning
	  fclose(fp);
	  irecv=2; // This is not used and only here to give irecv a valid value
	} else {
	  fp = fopen("problemwithrun","w");
	  fprintf(fp,"%d",isend);
  //   echo $signal > problemwithrun
	  fclose(fp);
	}

// 	while (access("runmodel", F_OK) != 0) {
// 	  sleep(waitTime);
// 	}

// 	fp = fopen("runmodel","r");
// 	fscanf(fp,"%d",&irecv);
// 	fclose(fp);

// 	system("rm -f runmodel");

  }

  MPI_Barrier(PETSC_COMM_WORLD);
  ierr = MPI_Bcast(&irecv,1,MPI_INT,activeRank,PETSC_COMM_WORLD);CHKERRMPI(ierr);

  activeRank++;
  activeRank = MIN(activeRank,numProcessors-1);
  
  return irecv;
}
