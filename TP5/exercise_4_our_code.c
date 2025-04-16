/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"
/*
***************************************
	PAY ATTENTION AT THE NEXT LINE!!
***************************************
*/
#include <math.h>

void copyCell(double* origin, double* destination){
	for(int i = 0; i<DIRECTIONS; i++) destination[i] = origin[i];
}

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	//
	// TODO: calculate the splitting parameters for the current task.
	//
	int rank,size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	// We can change this to be more performancing
	int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    comm->nb_x = dims[0];
    comm->nb_y = dims[1];
	const int factor = comm->nb_x /( (float) comm->nb_y);
	if(!(size % factor && total_width % total_height)){
		printf("Error: The total width must be a multiple of %d and the total height must be a multiple of %d\n", comm->nb_x, comm->nb_y);
		printf("\nSize = %d, Factor = %d\n",size,factor);
		exit(1);
	}

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank % comm->nb_x;
	comm->rank_y = rank / comm->nb_x;

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = total_width/(comm->nb_x) + 2;
	comm->height = total_height/(comm->nb_y) + 2;

	// TODO : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = (comm->rank_x)*(comm->width);
	comm->y = (comm->rank_y)*(comm->height);

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	//comm->buffer_recv_down, comm->buffer_recv_up, comm->buffer_send_down, comm->buffer_send_up
	comm->buffer_recv_down = (double*) malloc(sizeof(double)*DIRECTIONS*(comm->width));
	comm->buffer_send_down = (double*) malloc(sizeof(double)*DIRECTIONS*(comm->width));
	comm->buffer_recv_up = (double*) malloc(sizeof(double)*DIRECTIONS*(comm->width));
	comm->buffer_send_up = (double*) malloc(sizeof(double)*DIRECTIONS*(comm->width));
	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	free(comm->buffer_recv_down);
	free(comm->buffer_send_down);
	free(comm->buffer_recv_up);
	free(comm->buffer_send_up);
	//free allocated ressources
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate. 
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	/*
	************************************
	MAYBE WE CAN DO 
	enum{
		TAG_LR = 1, // ->
		TAG_RL,		// <-
		TAG_UD = 10,// \./
		TAG_DU,		// /'\
		TAG_45 = 20,// '/'
		TAG_135,	// '\'
		TAG_225,	// ./.
		TAG_315		// .\.
	};
	IN THE MAIN FILE DECLARATIONS (DEFINITIONS)
	*************************************
	*/
	const int TAG_LR = 1;
	const int TAG_RL = 2;
	const int TAG_UD = 10;
	const int TAG_DU = 11;
	const int TAG_45 = 20;
	const int TAG_135 = 21;
	const int TAG_225 = 22;
	const int TAG_315 = 23;
	
	int rank,size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	// Do we have to do this or is it already done?
	comm->rank_x = rank % comm->nb_x;
	comm->rank_y = rank % comm->nb_y;
	MPI_Status status;
	int left_rank = (comm->rank_x == 0) ? MPI_PROC_NULL : (comm->rank_x - 1) + (comm->rank_y)*(comm->nb_y);
	int right_rank = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : (comm->rank_x + 1) + (comm->rank_y)*(comm->nb_y);
	int up_rank = (comm->rank_y == 0) ? MPI_PROC_NULL : (comm->rank_y - 1)*(comm->nb_y) + (comm->rank_x);
	int down_rank = (comm->rank_y == comm->nb_y - 1) ? MPI_PROC_NULL : (comm->rank_y + 1)*(comm->nb_y) + (comm->rank_x);

	// Corners
	int c_45 = (right_rank != MPI_PROC_NULL) && (up_rank != MPI_PROC_NULL);
	int c_135 = (left_rank != MPI_PROC_NULL) && (up_rank != MPI_PROC_NULL);
	int c_225 = (left_rank != MPI_PROC_NULL) && (down_rank != MPI_PROC_NULL);
	int c_315 = (right_rank != MPI_PROC_NULL) && (down_rank != MPI_PROC_NULL);
	//====================================================
	// HORIZONTAL COMMUNICATIONS
	//====================================================
	//****************************************************
	// Communication (-->)
	//****************************************************
	if(c_45 || c_315) // send to right_rank
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2,c_45),(comm->height - 2 + !c_315 + !c_45)*DIRECTIONS,MPI_DOUBLE,right_rank,TAG_LR,comm->communicator);
	if(c_135 || c_225) // recieve from left_rank
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, c_135),(comm->height - 2 + !c_225 + !c_135)*DIRECTIONS,MPI_DOUBLE,left_rank,TAG_LR,comm->communicator,&status);
	//****************************************************
	// Communication (<--)
	//****************************************************
	if(c_135 || c_225) // send to left_rank
		MPI_Send(lbm_mesh_get_cell(mesh, 1, c_135),(comm->height - 2 + !c_225 + !c_135)*DIRECTIONS,MPI_DOUBLE,right_rank,TAG_RL,comm->communicator);
	if(c_45 || c_315) // recieve from right_rank
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, c_45),(comm->height - 2 + !c_315 + !c_45)*DIRECTIONS,MPI_DOUBLE,left_rank,TAG_RL,comm->communicator,&status);
	
	
	//====================================================
	// VERTICAL COMMUNICATIONS
	//====================================================
	//****************************************************
	// Communication (\./)
	//****************************************************
	if(c_225 || c_315){//send to down_rank
		int l = comm->height - 2 + !c_315 + !c_225;
		for(int x = 0; x < l; x++){
			copyCell(
				lbm_mesh_get_cell(mesh,x + c_225,comm->height - 2),
				comm->buffer_send_down + x
			);
		}
		MPI_Send(comm->buffer_send_down,l*DIRECTIONS,MPI_DOUBLE,down_rank,TAG_UD,comm->communicator);
	}
	if(c_45 || c_135){//recieve from up_rank
		int l = comm->height - 2 + !c_135 + !c_45;
		MPI_Recv(comm->buffer_recv_up,l*DIRECTIONS,MPI_DOUBLE,up_rank,TAG_UD,comm->communicator,&status);
		for(int x = 0; x < l; x++){
			copyCell(
				comm->buffer_recv_up + x,
				lbm_mesh_get_cell(mesh,x + c_135,1)
			);
		}
	}
	//****************************************************
	// Communication (/.\)
	//****************************************************
	if(c_45 || c_135){//send to up_rank
		int l = comm->height - 2 + !c_315 + !c_225;
		for(int x = 0; x < l; x++){
			copyCell(
				lbm_mesh_get_cell(mesh,x + c_135,1),
				comm->buffer_send_up + x
			);
		}
		MPI_Send(comm->buffer_send_up,l*DIRECTIONS,MPI_DOUBLE,up_rank,TAG_DU,comm->communicator);
	}
	if(c_225 || c_315){//recieve from down_rank
		int l = comm->height - 2 + !c_135 + !c_45;
		MPI_Recv(comm->buffer_recv_down,l*DIRECTIONS,MPI_DOUBLE,down_rank,TAG_DU,comm->communicator,&status);
		for(int x = 0; x < l; x++){
			copyCell(
				comm->buffer_recv_down + x,
				lbm_mesh_get_cell(mesh,x + c_225,comm->height - 2)
			);
		}
	}
	//====================================================
	// CORNERS
	//====================================================
	//****************************************************
	// Communication ( '/' )
	//****************************************************
	if(c_45){
		int rank_45 = up_rank + 1;
		double* cell = lbm_mesh_get_cell(mesh,comm->width - 2,1);
		MPI_Send(cell,DIRECTIONS,MPI_DOUBLE,rank_45,TAG_45,comm->communicator);
	}
	if(c_225){
		int rank_225 = down_rank - 1;
		double* cell = lbm_mesh_get_cell(mesh,0,comm->height - 1);
		MPI_Recv(cell,DIRECTIONS,MPI_DOUBLE,rank_225,TAG_45,comm->communicator,&status);
	}
	//****************************************************
	// Communication ( ./. )
	//****************************************************
	if(c_225){
		int rank_225 = down_rank - 1;
		double* cell = lbm_mesh_get_cell(mesh,1,comm->height - 2);
		MPI_Send(cell,DIRECTIONS,MPI_DOUBLE,rank_225,TAG_225,comm->communicator);
	}
	if(c_45){
		int rank_45 = up_rank + 1;
		double* cell = lbm_mesh_get_cell(mesh,comm->width - 1,0);
		MPI_Recv(cell,DIRECTIONS,MPI_DOUBLE,rank_45,TAG_225,comm->communicator,&status);
	}
	//****************************************************
	// Communication ( '\' )
	//****************************************************
	if(c_135){
		int rank_135 = up_rank - 1;
		double* cell = lbm_mesh_get_cell(mesh,1,1);
		MPI_Send(cell,DIRECTIONS,MPI_DOUBLE,rank_135,TAG_135,comm->communicator);
	}
	if(c_315){
		int rank_315 = down_rank + 1;
		double* cell = lbm_mesh_get_cell(mesh,comm->width - 1,comm->height - 1);
		MPI_Recv(cell,DIRECTIONS,MPI_DOUBLE,rank_315,TAG_135,comm->communicator,&status);
	}
	//****************************************************
	// Communication ( .\. )
	//****************************************************
	if(c_315){
		int rank_315 = down_rank + 1;
		double* cell = lbm_mesh_get_cell(mesh,comm->width - 2,comm->height - 2);
		MPI_Send(cell,DIRECTIONS,MPI_DOUBLE,rank_315,TAG_315,comm->communicator);
	}
	if(c_135){
		int rank_135 = up_rank - 1;
		double* cell = lbm_mesh_get_cell(mesh,0,0);
		MPI_Recv(cell,DIRECTIONS,MPI_DOUBLE,rank_135,TAG_315,comm->communicator,&status);
	}
}
