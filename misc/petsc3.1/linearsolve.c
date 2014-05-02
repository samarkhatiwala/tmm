static char help[] = "\n";

#include "petscksp.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Mat A;                /* matrix */
  Vec x,b;              /* vectors */
  KSP ksp;              /* linear solver context */
  PC pc;
  char matFile[PETSC_MAX_PATH_LEN];     /* input file name */
  char rhsFile[PETSC_MAX_PATH_LEN];
  char outFile[PETSC_MAX_PATH_LEN];  
  char iniFile[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  PetscInt n, its;
  PetscTruth flg;
  PetscTruth useMUMPS = PETSC_FALSE;
  PetscTruth useSuperLU = PETSC_FALSE;
  PetscViewer fd, fdrhs;
  PetscScalar  rnorm;
  PetscScalar t1, t2;
  PetscInt nrhs = 1;
  PetscInt itrhs;
  
  PetscInitialize(&argc,&args,(char *)0,help);

/* Some defaults */
  
/* Process options and load files */  

  ierr = PetscOptionsGetInt(PETSC_NULL,"-nrhs",&nrhs,&flg);CHKERRQ(ierr);

/* RHS     */
  ierr = PetscOptionsGetString(PETSC_NULL,"-r",rhsFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate binary RHS file with the -r option");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,rhsFile,FILE_MODE_READ,&fdrhs);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading RHS from file %s\n", rhsFile);CHKERRQ(ierr);
  ierr = VecLoad(fdrhs,PETSC_NULL,&b);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(fdrhs);CHKERRQ(ierr);  

/* Read matrix name */
  ierr = PetscOptionsGetString(PETSC_NULL,"-m",matFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate binary matrix file with the -m option");

/* Flag for direct solver */
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_mumps",&useMUMPS);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_superlu",&useSuperLU);CHKERRQ(ierr);

  if ((useMUMPS) && (useSuperLU)) {
    SETERRQ(1,"Cannot use both MUMPS and SuperLU!");
  }
  if (useMUMPS) ierr = PetscPrintf(PETSC_COMM_WORLD,"Direct solver (MUMPS) has been specified\n");CHKERRQ(ierr);
  if (useSuperLU) ierr = PetscPrintf(PETSC_COMM_WORLD,"Direct solver (SUPERLU) has been specified\n");CHKERRQ(ierr);

/* Read matrix */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading matrix from file %s\n", matFile);CHKERRQ(ierr);  
  if (useMUMPS) {
    ierr = MatLoad(fd,MATAIJ,&A);CHKERRQ(ierr);
  } else if (useSuperLU) {
    ierr = MatLoad(fd,MATAIJ,&A);CHKERRQ(ierr);  
  } else {
    ierr = MatLoad(fd,MATMPIAIJ,&A);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

  ierr = MatGetSize(A,0,&n);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix size: %d  x  %d\n", n,n);CHKERRQ(ierr);  

/* Output file */
  ierr = PetscOptionsGetString(PETSC_NULL,"-o",outFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate binary output file with the -o option");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be written to %s\n", outFile);CHKERRQ(ierr);  

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

/* -ksp_type preonly-pc_type lu-pc_ factor_mat_solver_package superlu_dist. */

  if (useMUMPS) {
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(pc,MAT_SOLVER_MUMPS);CHKERRQ(ierr);
  }  
  if (useSuperLU) {
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(pc,MAT_SOLVER_SUPERLU_DIST);CHKERRQ(ierr);
  }    
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

/* Solution vector*/
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

/* Optional initial guess (overwrite x) */
  ierr = PetscOptionsGetString(PETSC_NULL,"-i",iniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading initial guess from file %s\n", iniFile);CHKERRQ(ierr);	
	ierr = VecLoadIntoVector(fd,x);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);    
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,rhsFile,FILE_MODE_READ,&fdrhs);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);

  for (itrhs=1; itrhs<=nrhs; itrhs++) {
    ierr = VecLoad(fdrhs,PETSC_NULL,&b);CHKERRQ(ierr);  
    ierr = PetscGetTime(&t1);CHKERRQ(ierr); /* start counting wall clock time */        
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    ierr = PetscGetTime(&t2); CHKERRQ(ierr); /* stop counting wall clock time */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Wall clock time: %10.5f\n", t2-t1);CHKERRQ(ierr); 
  
    ierr = KSPGetResidualNorm(ksp,&rnorm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual norm=%15.10g\n",rnorm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of iterations=%d\n",its);CHKERRQ(ierr);
  
    ierr = VecView(x,fd);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(fdrhs);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);
  
  /* Free data structures */
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

