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

inline void MPI_Init(...) {return;}
inline void MPI_Finalize(...) {return;}
inline void MPI_Comm_rank(...) {return;}
inline void MPI_Comm_size(...) {return;}
inline void MPI_Bcast(...) {return;}
inline void MPI_Allgather(...) {return;}
inline void MPI_Allgatherv(...) {return;}
