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
#include <string.h>

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	int rank;
	int size;
	int dims[2] = {0, 0};
	int periods[2] = {0, 0};
	int coords[2];
	int local_width_without_ghost;
	int local_height_without_ghost;
	size_t row_buffer_size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Dims_create(size, 2, dims);
	comm->nb_y = dims[0];comm->nb_x = dims[1];
	if (total_width % comm->nb_x != 0 || total_height % comm->nb_y != 0)
		MPI_Abort(MPI_COMM_WORLD, -1);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm->communicator);
	MPI_Cart_coords(comm->communicator, rank, 2, coords);

comm->rank_y = coords[0];
	comm->rank_x = coords[1];
	local_width_without_ghost = total_width / comm->nb_x;
	local_height_without_ghost = total_height / comm->nb_y;

	comm->width = local_width_without_ghost +2;
	comm->height = local_height_without_ghost +2;
	comm->x = comm->rank_x *local_width_without_ghost;
	comm->y = comm->rank_y *local_height_without_ghost;
	row_buffer_size = comm->width * DIRECTIONS *sizeof(double);
	comm->buffer_send_up = malloc(row_buffer_size);
	comm->buffer_send_down = malloc(row_buffer_size);
	comm->buffer_recv_up = malloc(row_buffer_size);
	comm->buffer_recv_down = malloc(row_buffer_size);

if (comm->buffer_send_up == NULL || comm->buffer_send_down == NULL ||
	    comm->buffer_recv_up == NULL || comm->buffer_recv_down == NULL)
		fatal("Failed to allocate communication buffers for exercise 4.");

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	free(comm->buffer_send_up);
	free(comm->buffer_send_down);
	free(comm->buffer_recv_up);
	free(comm->buffer_recv_down);
	comm->buffer_send_up = NULL;
	comm->buffer_send_down = NULL;
	comm->buffer_recv_up = NULL;
	comm->buffer_recv_down = NULL;
	if (comm->communicator != MPI_COMM_NULL)
		MPI_Comm_free(&comm->communicator);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int coords[2];
	int left, right, up, down;
	int up_left, up_right, down_left, down_right;
	int rank;
	int x;
	const size_t cell_size = DIRECTIONS *sizeof(double);
	const int column_size = comm->height *DIRECTIONS;
	const int row_size = comm->width * DIRECTIONS;

	if (comm->rank_x - 1 < 0 || comm->rank_x - 1 >= comm->nb_x ||
	    comm->rank_y < 0 || comm->rank_y >= comm->nb_y)
		left = MPI_PROC_NULL;
	else
	{	coords[0] = comm->rank_y;
		coords[1] = comm->rank_x - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		left = rank;}

	if (comm->rank_x + 1 < 0 || comm->rank_x + 1 >= comm->nb_x ||
	    comm->rank_y < 0 || comm->rank_y >= comm->nb_y)
		right = MPI_PROC_NULL;
	else
	{	coords[0] = comm->rank_y;
		coords[1] = comm->rank_x + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		right = rank;
	}
	if (comm->rank_x < 0 || comm->rank_x >= comm->nb_x ||
	    comm->rank_y - 1 < 0 || comm->rank_y - 1 >= comm->nb_y)
		up = MPI_PROC_NULL;
	else
	{
		coords[0] = comm->rank_y - 1;
		coords[1] = comm->rank_x;
		MPI_Cart_rank(comm->communicator, coords,&rank);
		up = rank;}

	if (comm->rank_x < 0 || comm->rank_x >= comm->nb_x || comm->rank_y + 1 < 0 || comm->rank_y + 1 >= comm->nb_y){down = MPI_PROC_NULL;}
	else
	{	coords[0] = comm->rank_y + 1;
		coords[1] = comm->rank_x;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		down = rank;
	}
	if (comm->rank_x - 1 < 0 || comm->rank_x - 1 >= comm->nb_x || comm->rank_y - 1 < 0 || comm->rank_y - 1 >= comm->nb_y){up_left = MPI_PROC_NULL;}
	else{
		coords[0] = comm->rank_y - 1;
		coords[1] = comm->rank_x - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		up_left = rank;
	}
	if (comm->rank_x + 1 < 0 || comm->rank_x + 1 >= comm->nb_x || comm->rank_y - 1 < 0 || comm->rank_y - 1 >= comm->nb_y){up_right = MPI_PROC_NULL;}
	else
	{	coords[0] = comm->rank_y - 1;
		coords[1] = comm->rank_x + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		up_right = rank;
	}

	if (comm->rank_x - 1 < 0 || comm->rank_x - 1 >= comm->nb_x || comm->rank_y + 1 < 0 || comm->rank_y + 1 >= comm->nb_y){down_left = MPI_PROC_NULL;}
	else
	{	coords[0] = comm->rank_y + 1;
		coords[1] = comm->rank_x - 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		down_left = rank;
	}
	if (comm->rank_x + 1 < 0 || comm->rank_x + 1 >= comm->nb_x || comm->rank_y + 1 < 0 || comm->rank_y + 1 >= comm->nb_y){down_right = MPI_PROC_NULL;}
	else
	{
		coords[0] = comm->rank_y + 1;
		coords[1] = comm->rank_x + 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		down_right = rank;
	}
	MPI_Sendrecv(lbm_mesh_get_cell(mesh, 1, 0), column_size, MPI_DOUBLE,left,0, lbm_mesh_get_cell(mesh, comm->width - 1, 0),column_size, MPI_DOUBLE,right, 0, comm->communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(lbm_mesh_get_cell(mesh, comm->width - 2, 0),column_size, MPI_DOUBLE, right, 1,lbm_mesh_get_cell(mesh, 0, 0), column_size, MPI_DOUBLE,left,1, comm->communicator, MPI_STATUS_IGNORE);
	for (x = 0; x < comm->width; x++)
		{memcpy(&comm->buffer_send_up[x * DIRECTIONS], lbm_mesh_get_cell(mesh, x, 1), cell_size);}
	for (x = 0; x < comm->width; x++)
		{memcpy(&comm->buffer_send_down[x * DIRECTIONS], lbm_mesh_get_cell(mesh, x, comm->height - 2), cell_size);}
	MPI_Sendrecv(comm->buffer_send_up, row_size, MPI_DOUBLE, up, 2, comm->buffer_recv_down, row_size, MPI_DOUBLE, down, 2,comm->communicator, MPI_STATUS_IGNORE);
	if (down != MPI_PROC_NULL){
		for (x = 0; x < comm->width; x++)
			{memcpy(lbm_mesh_get_cell(mesh, x, comm->height - 1), &comm->buffer_recv_down[x * DIRECTIONS], cell_size);}
	}

	MPI_Sendrecv(comm->buffer_send_down, row_size,MPI_DOUBLE, down, 3, comm->buffer_recv_up, row_size, MPI_DOUBLE,up, 3, comm->communicator, MPI_STATUS_IGNORE);
	if (up != MPI_PROC_NULL){
		for (x = 0; x < comm->width; x++)
			{memcpy(lbm_mesh_get_cell(mesh, x, 0), &comm->buffer_recv_up[x * DIRECTIONS], cell_size);}
	}

	MPI_Sendrecv(lbm_mesh_get_cell(mesh, 1, 1),DIRECTIONS, MPI_DOUBLE, up_left,4,lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1), DIRECTIONS, MPI_DOUBLE,down_right, 4, comm->communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(lbm_mesh_get_cell(mesh, comm->width - 2, 1),DIRECTIONS, MPI_DOUBLE, up_right,5,lbm_mesh_get_cell(mesh, 0, comm->height - 1), DIRECTIONS, MPI_DOUBLE,down_left,5, comm->communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(lbm_mesh_get_cell(mesh, 1, comm->height - 2),DIRECTIONS, MPI_DOUBLE, down_left,6, lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS,MPI_DOUBLE,up_right,6, comm->communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(lbm_mesh_get_cell(mesh, comm->width - 2, comm->height - 2),DIRECTIONS,MPI_DOUBLE, down_right,7,lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS,MPI_DOUBLE,up_left,7, comm->communicator,MPI_STATUS_IGNORE);
}