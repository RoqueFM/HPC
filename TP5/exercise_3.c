/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with non-blocking MPI functions.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	MPI_Request request[4]; // Array of requests (2 for each direction)
	MPI_Status status[4]; // Array of statuses (2 for each direction)
	int request_count = 0; // Counter for number of requests

	const int TAG_LR = 1; // Left to Right
	const int TAG_RL = 2; // Right to Left

	int rank,size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	/*
	if(you have a neighbour at the left)
		left_rank = rank - 1
	else
		left_rank = MPI_PROC_NULL (you don't have any neighbours)

	(the same for the right neighbour)
	*/
	int left_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
	int right_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

	const int buffer_size = comm->height * DIRECTIONS; // Size of the buffer for ghost cell data

	// Data buffers for Ghost Cell data
	// double send_buffer[buffer_size];
	// double recv_buffer[buffer_size];

	if(right_rank != MPI_PROC_NULL){
		// Send to the right
		double* send_buffer = lbm_mesh_get_cell(mesh, comm->width - 2, 1); // Get the cell data
		MPI_Isend(send_buffer, buffer_size, MPI_DOUBLE, right_rank, TAG_LR, comm->communicator, request + request_count++);
	}
	if(left_rank != MPI_PROC_NULL){
		// Recieve from the left
		double* recv_buffer = lbm_mesh_get_cell(mesh, 0, 1); // Get the cell data
		MPI_Irecv(recv_buffer, buffer_size, MPI_DOUBLE, left_rank, TAG_LR, comm->communicator, request + request_count++);
	}
	if(left_rank != MPI_PROC_NULL){
		// Send to the left
		double* send_buffer = lbm_mesh_get_cell(mesh, 1, 1); // Get the cell data
		MPI_Isend(send_buffer, buffer_size, MPI_DOUBLE, left_rank, TAG_RL, comm->communicator, request + request_count++);
	}
	if(right_rank != MPI_PROC_NULL){
		// Recieve from the right
		double* recv_buffer = lbm_mesh_get_cell(mesh, comm->width - 1, 1); // Get the cell data
		MPI_Irecv(recv_buffer, buffer_size, MPI_DOUBLE, right_rank, TAG_RL, comm->communicator, request + request_count++);
	}
	// Wait for all requests to complete
	MPI_Waitall(request_count, request, status);	
}
