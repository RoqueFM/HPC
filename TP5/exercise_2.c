/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement odd/even 1D blocking communication scheme 
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
// NEW:
//     - >>> Odd/even communication ordering <<<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex2(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex2(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions using
	//       odd/even communications.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
	const int TAG_LR = 1; // Left to Right
	const int TAG_RL = 2; // Right to Left

	int rank,size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	MPI_Status status;
	int isOdd = rank % 2 == 1;

	/*
	if(you have a neighbour at the left)
		left_rank = rank - 1
	else
		left_rank = MPI_PROC_NULL (you don't have any neighbours)

	(the same for the right neighbour)
	*/
	int left_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
	int right_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

	// Data buffers for Ghost Cell data
	const int buffer_size = comm->height * DIRECTIONS;
	double* send_buffer;
	double* recv_buffer;
	// First phase odd ranks send, evens receive
	if(isOdd){
		// Send to the right recieve from left
		// SEND
		if(right_rank != MPI_PROC_NULL){
			send_buffer = lbm_mesh_get_cell(mesh, comm->width - 2,1);
			MPI_Send(send_buffer,buffer_size,MPI_DOUBLE,right_rank,TAG_LR,comm->communicator);
		}

		// SEND
		if(left_rank != MPI_PROC_NULL){
			send_buffer = lbm_mesh_get_cell(mesh, 1,1);
			MPI_Send(send_buffer,buffer_size,MPI_DOUBLE,left_rank,TAG_RL,comm->communicator);
		}
	}else{
		// RECIEVE
		if(left_rank != MPI_PROC_NULL){
			recv_buffer = lbm_mesh_get_cell(mesh, 0, 1);
			MPI_Recv(recv_buffer,buffer_size,MPI_DOUBLE,left_rank,TAG_LR,comm->communicator, &status);
		}
		// RECIEVE
		if(right_rank != MPI_PROC_NULL){
			recv_buffer = lbm_mesh_get_cell(mesh, comm->width - 1, 1);
			MPI_Recv(recv_buffer,buffer_size,MPI_DOUBLE,right_rank,TAG_RL,comm->communicator, &status);
		}	
	}

	// Second phase, odds receive, evens send
	if(isOdd){
		// RECIEVE
		if(left_rank != MPI_PROC_NULL){
			recv_buffer = lbm_mesh_get_cell(mesh,0,1);
			MPI_Recv(recv_buffer,buffer_size,MPI_DOUBLE,left_rank,TAG_LR,comm->communicator, &status);
		}
		// RECIEVE
		if(right_rank != MPI_PROC_NULL){
			recv_buffer = lbm_mesh_get_cell(mesh, comm->width - 1, 1);
			MPI_Recv(recv_buffer,buffer_size,MPI_DOUBLE,right_rank,TAG_RL,comm->communicator, &status);
		}
	}else{
		// SEND
		if(right_rank != MPI_PROC_NULL){
			send_buffer = lbm_mesh_get_cell(mesh, comm->width - 2,1);
			MPI_Send(send_buffer,buffer_size,MPI_DOUBLE,right_rank,TAG_LR,comm->communicator);
		}

		// SEND
		if(left_rank != MPI_PROC_NULL){
			send_buffer = lbm_mesh_get_cell(mesh, 1,1);
			MPI_Send(send_buffer,buffer_size,MPI_DOUBLE,left_rank,TAG_RL,comm->communicator);
		}
	}
}
