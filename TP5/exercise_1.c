/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
//
// GOAL: Implement a 1D communication scheme along
//       X axis with blocking communications.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex1(lbm_comm_t * comm, int total_width, int total_height)
{
	//
	// TODO: calculate the splitting parameters for the current task.
	//
	// HINT: You can look in exercise_0.c to get an example for the sequential case.
	//
	int rank,size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	// size = size;
	// TODO: calculate the number of tasks along X axis and Y axis.
	comm->nb_x = size;
	comm->nb_y = 1;

	if(total_width % (comm->nb_x) != 0){
		// std::cout << "The total width must be a multiple of " << comm->nb_x << endl;
		exit(1);
	}

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank;
	comm->rank_y = 0;

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = total_width/(comm->nb_x) + 2;
	// comm->height = total_height/(comm->nb_y) + 2
	comm->height = total_height + 2;

	// TODO : calculate the absolute position in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = (comm->rank_x)*total_width/(comm->nb_x);
	comm->y = 0;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions (MPI_Send & MPI_Recv)
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[DIRECTIONS] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
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
	// double send_buffer[comm->height * DIRECTIONS];
	// double recv_buffer[comm->height * DIRECTIONS];
	double* send_buffer;
	double* recv_buffer;

	// Send to the right recieve from left
	// SEND
	if(right_rank != MPI_PROC_NULL){
		// for(int y = 0; y<comm->height; y++){
		// 	// Get the cell data (9 directions)
		// 	double* cell = lbm_mesh_get_cell(mesh, comm->width - 2,y);
		// 	// Store it in the send_buffer
		// 	for(int direction = 0; direction < DIRECTIONS; direction++)
		// 		send_buffer[y*DIRECTIONS + direction] = cell[direction];
		// }
		send_buffer = lbm_mesh_get_cell(mesh, comm->width - 2, 1); // Get the cell data
		/*
			The cell that we are sending is the one at (comm->width - 2, 1)
			That means, the cell which is taken account for the calculation of the rank
			We're not taking the ghost cell which is located at (comm->width - 1, 0)
			The send buffer is the cell at (comm->width - 2, 1)
			(the pointer of the buffer is pointing to it because the colum-major order)
			The y direction is the row-major order in the memory (not physically)
			As the MPI_Send function reads the memory location, it gets all the next elements
			which are the next cells in the same column
			That means, the send_buffer is pointing to the cell at (comm->width - 2, 1)
			And the next elements are the cells at (comm->width - 2, 2), (comm->width - 2, 3), ...
			And so on...
			And so on...

			To clarify this, see figure 7 in the lab worksheet
		*/
		// Send the data
		MPI_Send(send_buffer,comm->height*DIRECTIONS,MPI_DOUBLE,right_rank,TAG_LR,comm->communicator);
	}

	// RECIEVE
	if(left_rank != MPI_PROC_NULL){
		// Recieve the data
		recv_buffer = lbm_mesh_get_cell(mesh, 0, 1); // Get the cell data (9 directions)
		MPI_Recv(recv_buffer,comm->height*DIRECTIONS,MPI_DOUBLE,left_rank,TAG_LR,comm->communicator, &status);
		// for(int y = 0; y<comm->height; y++){
		// 	// Get the cell location
		// 	double* cell = lbm_mesh_get_cell(mesh, 0, y); // First to the left
		// 	// Store it in the cell
		// 	for(int direction = 0; direction < DIRECTIONS; direction++)
		// 		cell[direction] = recv_buffer[y*DIRECTIONS + direction];
		// }

	}

	// Send to the left recieve from right
	// SEND
	if(left_rank != MPI_PROC_NULL){
		// for(int y = 0; y<comm->height; y++){
		// 	// Get the cell data (9 directions)
		// 	double* cell = lbm_mesh_get_cell(mesh, 1,y); // Sending from the first real cell on the left (not the ghost one)
		// 	// Store it in the send_buffer
		// 	for(int direction = 0; direction < DIRECTIONS; direction++)
		// 		send_buffer[y*DIRECTIONS + direction] = cell[direction];
		// }
		send_buffer = lbm_mesh_get_cell(mesh, 1, 1); // Get the cell data (9 directions)
		// Send the data
		MPI_Send(send_buffer,comm->height*DIRECTIONS,MPI_DOUBLE,left_rank,TAG_RL,comm->communicator);
	}

	// RECIEVE
	if(right_rank != MPI_PROC_NULL){
		// Recieve the data
		recv_buffer = lbm_mesh_get_cell(mesh, comm->width - 1, 1); // Get the cell data (9 directions)
		MPI_Recv(recv_buffer,comm->height*DIRECTIONS,MPI_DOUBLE,right_rank,TAG_RL,comm->communicator, &status);
		// for(int y = 0; y<comm->height; y++){
		// 	// Get the cell location
		// 	double* cell = lbm_mesh_get_cell(mesh, comm->width - 1, y); // comm->width - 2 + 1
		// 	// Store it in the cell
		// 	for(int direction = 0; direction < DIRECTIONS; direction++)
		// 		cell[direction] = recv_buffer[y*DIRECTIONS + direction];
		// }
	}
}
