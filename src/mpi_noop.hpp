/**
 * @file mpi_noop.hpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Noop definitions for MPI for when compiling one of the non-MPI
 * versions to run on a machine without MPI installed
 */

#pragma once

#define MPI_INT 0
#define MPI_COMM_WORLD 0
#define MPI_BYTE 0
#define MPI_C_BOOL 0

void MPI_Init(...) {return;}
void MPI_Finalize(...) {return;}
void MPI_Comm_rank(...) {return;}
void MPI_Comm_size(...) {return;}
void MPI_Bcast(...) {return;}
void MPI_Allgather(...) {return;}
void MPI_Allgatherv(...) {return;}
