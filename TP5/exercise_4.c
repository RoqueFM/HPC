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

/****************************************************/
// Helper function to copy cell data
void copy_cell(double* src, double* dst) {
    for (int i = 0; i < DIRECTIONS; i++) {
        dst[i] = src[i];
    }
}

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
    //
    // Calculate the splitting parameters for the current task.
    //
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate the number of tasks along X axis and Y axis.
    // Use MPI_Dims_create for optimal decomposition
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    comm->nb_x = dims[0];
    comm->nb_y = dims[1];

    // Calculate the current task position in the splitting
    comm->rank_x = rank % comm->nb_x;
    comm->rank_y = rank / comm->nb_x;

    // Calculate the local sub-domain size (including ghost cells)
    comm->width = total_width / comm->nb_x + 2;
    comm->height = total_height / comm->nb_y + 2;

    // Calculate the absolute position (in cell number) in the global mesh
    // without accounting for ghost cells
    comm->x = comm->rank_x * (comm->width - 2);
    comm->y = comm->rank_y * (comm->height - 2);

    // Allocate temporary copy buffers
    // We need buffers for top/bottom rows (width * DIRECTIONS per row)
    comm->buffer_recv_down = (double*)malloc(sizeof(double) * DIRECTIONS * comm->width);
    comm->buffer_send_down = (double*)malloc(sizeof(double) * DIRECTIONS * comm->width);
    comm->buffer_recv_up = (double*)malloc(sizeof(double) * DIRECTIONS * comm->width);
    comm->buffer_send_up = (double*)malloc(sizeof(double) * DIRECTIONS * comm->width);

    // Store MPI_COMM_WORLD in the communicator field (or create a cartesian topology)
    comm->communicator = MPI_COMM_WORLD;

    // Optional: Create cartesian communicator
    /*
    int periods[2] = {0, 0}; // Non-periodic mesh as specified
    int coords[2];
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm->communicator);
    MPI_Cart_coords(comm->communicator, rank, 2, coords);
    
    // Verify coordinates (should match our calculation)
    comm->rank_x = coords[0];
    comm->rank_y = coords[1];
    */

    // Debug output
    // lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
    // Free allocated resources
    free(comm->buffer_recv_down);
    free(comm->buffer_send_down);
    free(comm->buffer_recv_up);
    free(comm->buffer_send_up);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
    // Get rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // Calculate neighbor ranks
    // Use MPI_PROC_NULL for border processes
    int left_rank = (comm->rank_x == 0) ? MPI_PROC_NULL : rank - 1;
    int right_rank = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : rank + 1;
    int up_rank = (comm->rank_y == 0) ? MPI_PROC_NULL : rank - comm->nb_x;
    int down_rank = (comm->rank_y == comm->nb_y - 1) ? MPI_PROC_NULL : rank + comm->nb_x;
    
    // Calculate diagonal neighbor ranks
    int up_left_rank = (up_rank != MPI_PROC_NULL && left_rank != MPI_PROC_NULL) ? up_rank - 1 : MPI_PROC_NULL;
    int up_right_rank = (up_rank != MPI_PROC_NULL && right_rank != MPI_PROC_NULL) ? up_rank + 1 : MPI_PROC_NULL;
    int down_left_rank = (down_rank != MPI_PROC_NULL && left_rank != MPI_PROC_NULL) ? down_rank - 1 : MPI_PROC_NULL;
    int down_right_rank = (down_rank != MPI_PROC_NULL && right_rank != MPI_PROC_NULL) ? down_rank + 1 : MPI_PROC_NULL;

    // Define tags for different communications
    const int TAG_LEFT_RIGHT = 1;
    const int TAG_RIGHT_LEFT = 2;
    const int TAG_TOP_BOTTOM = 3;
    const int TAG_BOTTOM_TOP = 4;
    const int TAG_DIAG_UP_LEFT = 5;
    const int TAG_DIAG_UP_RIGHT = 6;
    const int TAG_DIAG_DOWN_LEFT = 7;
    const int TAG_DIAG_DOWN_RIGHT = 8;

    //==================================================
    // STEP 1: LEFT/RIGHT COMMUNICATIONS (contiguous data)
    //==================================================
    
    // Send right column, receive left ghost column
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, comm->width - 2, 1), (comm->height - 2) * DIRECTIONS, MPI_DOUBLE, right_rank, TAG_LEFT_RIGHT,
        lbm_mesh_get_cell(mesh, 0, 1), (comm->height - 2) * DIRECTIONS, MPI_DOUBLE, left_rank, TAG_LEFT_RIGHT,
        comm->communicator, &status
    );
    
    // Send left column, receive right ghost column
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, 1, 1), (comm->height - 2) * DIRECTIONS, MPI_DOUBLE, left_rank, TAG_RIGHT_LEFT,
        lbm_mesh_get_cell(mesh, comm->width - 1, 1), (comm->height - 2) * DIRECTIONS, MPI_DOUBLE, right_rank, TAG_RIGHT_LEFT,
        comm->communicator, &status
    );

    //==================================================
    // STEP 2: TOP/BOTTOM COMMUNICATIONS (non-contiguous data)
    //==================================================
    
    // Copy bottom row (excluding corners) to send buffer
    for (int x = 1; x < comm->width - 1; x++) {
        double *src = lbm_mesh_get_cell(mesh, x, comm->height - 2);
        double *dst = &comm->buffer_send_down[(x-1) * DIRECTIONS];
        copy_cell(src, dst);
    }
    
    // Copy top row (excluding corners) to send buffer
    for (int x = 1; x < comm->width - 1; x++) {
        double *src = lbm_mesh_get_cell(mesh, x, 1);
        double *dst = &comm->buffer_send_up[(x-1) * DIRECTIONS];
        copy_cell(src, dst);
    }
    
    // Send bottom row, receive top ghost row
    MPI_Sendrecv(
        comm->buffer_send_down, (comm->width - 2) * DIRECTIONS, MPI_DOUBLE, down_rank, TAG_TOP_BOTTOM,
        comm->buffer_recv_up, (comm->width - 2) * DIRECTIONS, MPI_DOUBLE, up_rank, TAG_TOP_BOTTOM,
        comm->communicator, &status
    );
    
    // Send top row, receive bottom ghost row
    MPI_Sendrecv(
        comm->buffer_send_up, (comm->width - 2) * DIRECTIONS, MPI_DOUBLE, up_rank, TAG_BOTTOM_TOP,
        comm->buffer_recv_down, (comm->width - 2) * DIRECTIONS, MPI_DOUBLE, down_rank, TAG_BOTTOM_TOP,
        comm->communicator, &status
    );
    
    // Copy received data from buffer to ghost cells (top row)
    for (int x = 1; x < comm->width - 1; x++) {
        double *src = &comm->buffer_recv_up[(x-1) * DIRECTIONS];
        double *dst = lbm_mesh_get_cell(mesh, x, 0);
        copy_cell(src, dst);
    }
    
    // Copy received data from buffer to ghost cells (bottom row)
    for (int x = 1; x < comm->width - 1; x++) {
        double *src = &comm->buffer_recv_down[(x-1) * DIRECTIONS];
        double *dst = lbm_mesh_get_cell(mesh, x, comm->height - 1);
        copy_cell(src, dst);
    }

    //==================================================
    // STEP 3: DIAGONAL COMMUNICATIONS (corners)
    //==================================================
    
    // Top-left corner
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, 1, 1), DIRECTIONS, MPI_DOUBLE, up_left_rank, TAG_DIAG_UP_LEFT,
        lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS, MPI_DOUBLE, up_left_rank, TAG_DIAG_UP_LEFT,
        comm->communicator, &status
    );
    
    // Top-right corner
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, comm->width - 2, 1), DIRECTIONS, MPI_DOUBLE, up_right_rank, TAG_DIAG_UP_RIGHT,
        lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS, MPI_DOUBLE, up_right_rank, TAG_DIAG_UP_RIGHT,
        comm->communicator, &status
    );
    
    // Bottom-left corner
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, 1, comm->height - 2), DIRECTIONS, MPI_DOUBLE, down_left_rank, TAG_DIAG_DOWN_LEFT,
        lbm_mesh_get_cell(mesh, 0, comm->height - 1), DIRECTIONS, MPI_DOUBLE, down_left_rank, TAG_DIAG_DOWN_LEFT,
        comm->communicator, &status
    );
    
    // Bottom-right corner
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, comm->width - 2, comm->height - 2), DIRECTIONS, MPI_DOUBLE, down_right_rank, TAG_DIAG_DOWN_RIGHT,
        lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1), DIRECTIONS, MPI_DOUBLE, down_right_rank, TAG_DIAG_DOWN_RIGHT,
        comm->communicator, &status
    );
}