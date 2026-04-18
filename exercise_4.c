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
static int lbm_comm_rank_from_coords_ex4(lbm_comm_t * comm, int x, int y)
{
	int coords[2];
	int rank;

	if (x < 0 || x >= comm->nb_x || y < 0 || y >= comm->nb_y)
		return MPI_PROC_NULL;

	coords[0] = y;
	coords[1] = x;
	MPI_Cart_rank(comm->communicator, coords, &rank);
	return rank;
}

/****************************************************/
static void lbm_comm_pack_row_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh, int y, double * buffer)
{
	int x;
	const size_t cell_size = DIRECTIONS * sizeof(double);

	for (x = 0; x < comm->width; x++)
		memcpy(&buffer[x * DIRECTIONS], lbm_mesh_get_cell(mesh, x, y), cell_size);
}

/****************************************************/
static void lbm_comm_unpack_row_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh, int y, double * buffer)
{
	int x;
	const size_t cell_size = DIRECTIONS * sizeof(double);

	for (x = 0; x < comm->width; x++)
		memcpy(lbm_mesh_get_cell(mesh, x, y), &buffer[x * DIRECTIONS], cell_size);
}

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

	comm->nb_y = dims[0];
	comm->nb_x = dims[1];

	if (total_width % comm->nb_x != 0 || total_height % comm->nb_y != 0)
		MPI_Abort(MPI_COMM_WORLD, -1);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm->communicator);
	MPI_Cart_coords(comm->communicator, rank, 2, coords);

	comm->rank_y = coords[0];
	comm->rank_x = coords[1];

	local_width_without_ghost = total_width / comm->nb_x;
	local_height_without_ghost = total_height / comm->nb_y;

	comm->width = local_width_without_ghost + 2;
	comm->height = local_height_without_ghost + 2;

	comm->x = comm->rank_x * local_width_without_ghost;
	comm->y = comm->rank_y * local_height_without_ghost;

	row_buffer_size = comm->width * DIRECTIONS * sizeof(double);
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
	int left_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x - 1, comm->rank_y);
	int right_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x + 1, comm->rank_y);
	int up_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x, comm->rank_y - 1);
	int down_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x, comm->rank_y + 1);
	int up_left_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x - 1, comm->rank_y - 1);
	int up_right_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x + 1, comm->rank_y - 1);
	int down_left_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x - 1, comm->rank_y + 1);
	int down_right_rank = lbm_comm_rank_from_coords_ex4(comm, comm->rank_x + 1, comm->rank_y + 1);
	const int column_size = comm->height * DIRECTIONS;
	const int row_size = comm->width * DIRECTIONS;

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 1, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		0,
		lbm_mesh_get_cell(mesh, comm->width - 1, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		0,
		comm->communicator,
		MPI_STATUS_IGNORE
	);

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, comm->width - 2, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		1,
		lbm_mesh_get_cell(mesh, 0, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		1,
		comm->communicator,
		MPI_STATUS_IGNORE
	);

	lbm_comm_pack_row_ex4(comm, mesh, 1, comm->buffer_send_up);
	lbm_comm_pack_row_ex4(comm, mesh, comm->height - 2, comm->buffer_send_down);

	MPI_Sendrecv(
		comm->buffer_send_up,
		row_size,
		MPI_DOUBLE,
		up_rank,
		2,
		comm->buffer_recv_down,
		row_size,
		MPI_DOUBLE,
		down_rank,
		2,
		comm->communicator,
		MPI_STATUS_IGNORE
	);
	if (down_rank != MPI_PROC_NULL)
		lbm_comm_unpack_row_ex4(comm, mesh, comm->height - 1, comm->buffer_recv_down);

	MPI_Sendrecv(
		comm->buffer_send_down,
		row_size,
		MPI_DOUBLE,
		down_rank,
		3,
		comm->buffer_recv_up,
		row_size,
		MPI_DOUBLE,
		up_rank,
		3,
		comm->communicator,
		MPI_STATUS_IGNORE
	);
	if (up_rank != MPI_PROC_NULL)
		lbm_comm_unpack_row_ex4(comm, mesh, 0, comm->buffer_recv_up);

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 1, 1),
		DIRECTIONS,
		MPI_DOUBLE,
		up_left_rank,
		4,
		lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1),
		DIRECTIONS,
		MPI_DOUBLE,
		down_right_rank,
		4,
		comm->communicator,
		MPI_STATUS_IGNORE
	);

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, comm->width - 2, 1),
		DIRECTIONS,
		MPI_DOUBLE,
		up_right_rank,
		5,
		lbm_mesh_get_cell(mesh, 0, comm->height - 1),
		DIRECTIONS,
		MPI_DOUBLE,
		down_left_rank,
		5,
		comm->communicator,
		MPI_STATUS_IGNORE
	);

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 1, comm->height - 2),
		DIRECTIONS,
		MPI_DOUBLE,
		down_left_rank,
		6,
		lbm_mesh_get_cell(mesh, comm->width - 1, 0),
		DIRECTIONS,
		MPI_DOUBLE,
		up_right_rank,
		6,
		comm->communicator,
		MPI_STATUS_IGNORE
	);

	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, comm->width - 2, comm->height - 2),
		DIRECTIONS,
		MPI_DOUBLE,
		down_right_rank,
		7,
		lbm_mesh_get_cell(mesh, 0, 0),
		DIRECTIONS,
		MPI_DOUBLE,
		up_left_rank,
		7,
		comm->communicator,
		MPI_STATUS_IGNORE
	);
}
