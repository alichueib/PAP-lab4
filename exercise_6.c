/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication with non-blocking
//       messages.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - MPI type for non contiguous cells
// NEW:
//     - Non-blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
static int lbm_comm_rank_from_coords_ex6(lbm_comm_t * comm, int x, int y)
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
void lbm_comm_init_ex6(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5
	lbm_comm_init_ex5(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_release_ex6(lbm_comm_t * comm)
{
	//we use the same implementation than ext 5
	lbm_comm_release_ex5(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex6(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int left_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x - 1, comm->rank_y);
	int right_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x + 1, comm->rank_y);
	int up_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x, comm->rank_y - 1);
	int down_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x, comm->rank_y + 1);
	int up_left_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x - 1, comm->rank_y - 1);
	int up_right_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x + 1, comm->rank_y - 1);
	int down_left_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x - 1, comm->rank_y + 1);
	int down_right_rank = lbm_comm_rank_from_coords_ex6(comm, comm->rank_x + 1, comm->rank_y + 1);
	int req_count = 0;
	const int column_size = comm->height * DIRECTIONS;

	MPI_Irecv(
		lbm_mesh_get_cell(mesh, comm->width - 1, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		0,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		1,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 1, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		0,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, comm->width - 2, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		1,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Waitall(req_count, comm->requests, MPI_STATUSES_IGNORE);

	req_count = 0;
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, comm->height - 1),
		1,
		comm->type,
		down_rank,
		2,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, 0),
		1,
		comm->type,
		up_rank,
		3,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 0, 1),
		1,
		comm->type,
		up_rank,
		2,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 0, comm->height - 2),
		1,
		comm->type,
		down_rank,
		3,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Waitall(req_count, comm->requests, MPI_STATUSES_IGNORE);

	req_count = 0;
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1),
		DIRECTIONS,
		MPI_DOUBLE,
		down_right_rank,
		4,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, comm->height - 1),
		DIRECTIONS,
		MPI_DOUBLE,
		down_left_rank,
		5,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, comm->width - 1, 0),
		DIRECTIONS,
		MPI_DOUBLE,
		up_right_rank,
		6,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, 0),
		DIRECTIONS,
		MPI_DOUBLE,
		up_left_rank,
		7,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 1, 1),
		DIRECTIONS,
		MPI_DOUBLE,
		up_left_rank,
		4,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, comm->width - 2, 1),
		DIRECTIONS,
		MPI_DOUBLE,
		up_right_rank,
		5,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 1, comm->height - 2),
		DIRECTIONS,
		MPI_DOUBLE,
		down_left_rank,
		6,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Isend(
		lbm_mesh_get_cell(mesh, comm->width - 2, comm->height - 2),
		DIRECTIONS,
		MPI_DOUBLE,
		down_right_rank,
		7,
		comm->communicator,
		&comm->requests[req_count++]
	);
	MPI_Waitall(req_count, comm->requests, MPI_STATUSES_IGNORE);
}
