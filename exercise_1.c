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
	int rank;
	int size;
	int local_width_without_ghost;

	// Get rank information from the global communicator.
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// 1D decomposition: all ranks are arranged along X only.
	comm->nb_x = size;
	comm->nb_y = 1;

	comm->rank_x = rank;
	comm->rank_y = 0;

	// The global width must be divisible by the number of X sub-domains.
	if (total_width % size != 0)
		MPI_Abort(MPI_COMM_WORLD, -1);

	local_width_without_ghost = total_width / size;

	// Add one ghost column on each side and one ghost row on top/bottom.
	comm->width = local_width_without_ghost + 2;
	comm->height = total_height + 2;

	// Absolute position in the global mesh, excluding ghosts.
	comm->x = rank * local_width_without_ghost;
	comm->y = 0;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int left_rank;
	int right_rank;
	const int column_size = comm->height * DIRECTIONS;

	left_rank = (comm->rank_x == 0) ? MPI_PROC_NULL : comm->rank_x - 1;
	right_rank = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : comm->rank_x + 1;

	// Send our last interior column (x = width - 2) to the right neighbor and
	// receive the left neighbor interior column into our left ghost (x = 0).
	if (left_rank != MPI_PROC_NULL)
		MPI_Recv(
			lbm_mesh_get_cell(mesh, 0, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);
	if (right_rank != MPI_PROC_NULL)
		MPI_Send(
			lbm_mesh_get_cell(mesh, comm->width - 2, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			0,
			MPI_COMM_WORLD
		);

	// Send our first interior column (x = 1) to the left neighbor and
	// receive the right neighbor interior column into our right ghost (x = width - 1).
	if (right_rank != MPI_PROC_NULL)
		MPI_Recv(
			lbm_mesh_get_cell(mesh, comm->width - 1, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			1,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);
	if (left_rank != MPI_PROC_NULL)
		MPI_Send(
			lbm_mesh_get_cell(mesh, 1, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			1,
			MPI_COMM_WORLD
		);
}
