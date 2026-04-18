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
	int rank;
	int left_rank;
	int right_rank;
	const int column_size = comm->height * DIRECTIONS;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	left_rank = (comm->rank_x == 0) ? MPI_PROC_NULL : comm->rank_x - 1;
	right_rank = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : comm->rank_x + 1;

	if ((rank % 2) == 0)
	{
		// Even ranks first receive the left-going column from the right neighbor
		// into their right ghost column.
		MPI_Recv(
			lbm_mesh_get_cell(mesh, comm->width - 1, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);

		// Then they send their right interior column to the right neighbor,
		// which will be received into its left ghost column.
		MPI_Send(
			lbm_mesh_get_cell(mesh, comm->width - 2, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			1,
			MPI_COMM_WORLD
		);

		// Then they send their left interior column to the left neighbor,
		// which will be received into its right ghost column.
		MPI_Send(
			lbm_mesh_get_cell(mesh, 1, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			0,
			MPI_COMM_WORLD
		);

		// Finally they receive the right-going column from the left neighbor
		// into their left ghost column.
		MPI_Recv(
			lbm_mesh_get_cell(mesh, 0, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			1,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);
	}
	else
	{
		// Odd ranks first send their left interior column to the left neighbor,
		// which will be received into its right ghost column.
		MPI_Send(
			lbm_mesh_get_cell(mesh, 1, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			0,
			MPI_COMM_WORLD
		);

		// Then they receive the right-going column from the left neighbor
		// into their left ghost column.
		MPI_Recv(
			lbm_mesh_get_cell(mesh, 0, 0),
			column_size,
			MPI_DOUBLE,
			left_rank,
			1,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);

		// Then they receive the left-going column from the right neighbor
		// into their right ghost column.
		MPI_Recv(
			lbm_mesh_get_cell(mesh, comm->width - 1, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			0,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE
		);

		// Finally they send their right interior column to the right neighbor,
		// which will be received into its left ghost column.
		MPI_Send(
			lbm_mesh_get_cell(mesh, comm->width - 2, 0),
			column_size,
			MPI_DOUBLE,
			right_rank,
			1,
			MPI_COMM_WORLD
		);
	}
}
