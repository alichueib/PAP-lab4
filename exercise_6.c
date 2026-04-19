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
void lbm_comm_init_ex6(lbm_comm_t * comm, int total_width, int total_height){
	lbm_comm_init_ex5(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_release_ex6(lbm_comm_t * comm){
	lbm_comm_release_ex5(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex6(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int coords[2];
	int rank;
	int left, right, up, down;
	int up_left, up_right, down_left, down_right;
	int req_count=0;
	const int column_size=comm->height*DIRECTIONS;
	if(comm->rank_x-1<0|| comm->rank_x-1>=comm->nb_x|| comm->rank_y<0 || comm->rank_y>=comm->nb_y){
		left=MPI_PROC_NULL;}
	else{
		coords[0]=comm->rank_y;
		coords[1]=comm->rank_x-1;
		MPI_Cart_rank(comm->communicator,coords,&rank);
		left=rank;
	}
	if(comm->rank_x+1<0|| comm->rank_x+1>=comm->nb_x || comm->rank_y<0|| comm->rank_y>=comm->nb_y){
		right=MPI_PROC_NULL;}
	else{
		coords[0]=comm->rank_y;
		coords[1]=comm->rank_x+1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		right=rank;
	}

	if(comm->rank_x<0|| comm->rank_x>=comm->nb_x || comm->rank_y-1< 0 || comm->rank_y-1>= comm->nb_y){
		up=MPI_PROC_NULL;}
	else{
		coords[0]=comm->rank_y-1;
		coords[1]=comm->rank_x;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		up = rank;
	}

	if(comm->rank_x<0||comm->rank_x>=comm->nb_x|| comm->rank_y+1<0 || comm->rank_y+1 >= comm->nb_y){
		down=MPI_PROC_NULL;} 
	else{
		coords[0] = comm->rank_y+ 1;
		coords[1]= comm->rank_x;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		down=rank;
	}

	if(comm->rank_x-1<0||comm->rank_x-1>=comm->nb_x||comm->rank_y-1<0||comm->rank_y-1>=comm->nb_y){
		up_left=MPI_PROC_NULL;}
	else{
		coords[0] = comm->rank_y- 1;
		coords[1] = comm-> rank_x- 1;
		MPI_Cart_rank(comm->communicator, coords, &rank);
		up_left= rank;
	}

	if(comm->rank_x+1<0||comm->rank_x + 1 >= comm->nb_x || comm->rank_y-1<0 || comm->rank_y-1>= comm->nb_y){
		up_right= MPI_PROC_NULL;}
	else{
		coords[0]= comm->rank_y-1;
		coords[1] =comm->rank_x+1;
		MPI_Cart_rank(comm-> communicator, coords,&rank);
		up_right = rank;
	}

	if(comm->rank_x -1 <0 || comm ->rank_x- 1 >=comm -> nb_x || comm-> rank_y+1< 0|| comm-> rank_y+ 1 >=comm-> nb_y){
		down_left = MPI_PROC_NULL;}
	else{
		coords[0]=comm->rank_y+1;
		coords[1]=comm->rank_x-1;
		MPI_Cart_rank(comm->communicator,coords,&rank);
		down_left=rank;
	}

	if(comm->rank_x+1<0 || comm-> rank_x+1 >= comm-> nb_x || comm->rank_y +1 < 0 || comm->rank_y+1 >= comm-> nb_y){
		down_right = MPI_PROC_NULL;}
	else{
		coords[0]= comm-> rank_y+1;
		coords[1]= comm-> rank_x+1;
		MPI_Cart_rank(comm->communicator, coords,&rank);
		down_right = rank;
	}

	MPI_Irecv( lbm_mesh_get_cell(mesh,comm->width-1, 0), column_size,MPI_DOUBLE, right,0, comm->communicator,&comm->requests[req_count++] );
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0,0),column_size, MPI_DOUBLE,left, 1,comm->communicator, &comm->requests[req_count++] );
	MPI_Isend( lbm_mesh_get_cell(mesh,1, 0), column_size, MPI_DOUBLE,left,0,comm->communicator,&comm->requests[req_count++] );
	MPI_Isend(lbm_mesh_get_cell(mesh, comm->width-2,0),column_size,MPI_DOUBLE, right,1, comm->communicator, &comm->requests[req_count++] );
	MPI_Waitall(req_count,comm->requests, MPI_STATUSES_IGNORE);

	req_count = 0;
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0,comm->height-1),1,comm->type, down,2,comm->communicator,&comm->requests[req_count++] );
	MPI_Irecv( lbm_mesh_get_cell(mesh,0, 0), 1, comm->type,up,3, comm->communicator,&comm->requests[req_count++] );
	MPI_Isend(lbm_mesh_get_cell(mesh, 0,1),1, comm->type, up,2,comm->communicator, &comm->requests[req_count++] );
	MPI_Isend( lbm_mesh_get_cell(mesh,0,comm->height-2), 1,comm->type,down, 3,comm->communicator,&comm->requests[req_count++] );
	MPI_Waitall( req_count, comm->requests,MPI_STATUSES_IGNORE );

	req_count=0;
	MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width-1,comm->height-1),DIRECTIONS, MPI_DOUBLE,down_right,4,comm->communicator,&comm->requests[req_count++] );
	MPI_Irecv( lbm_mesh_get_cell(mesh,0,comm->height-1), DIRECTIONS,MPI_DOUBLE, down_left,5, comm->communicator,&comm->requests[req_count++] );
	MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width-1,0), DIRECTIONS, MPI_DOUBLE, up_right,6,comm->communicator,&comm->requests[req_count++] );
	MPI_Irecv( lbm_mesh_get_cell(mesh,0,0),DIRECTIONS,MPI_DOUBLE,up_left,7,comm->communicator, &comm->requests[req_count++] );
	MPI_Isend(lbm_mesh_get_cell(mesh, 1,1), DIRECTIONS, MPI_DOUBLE, up_left,4,comm->communicator,&comm->requests[req_count++] );
	MPI_Isend( lbm_mesh_get_cell(mesh,comm->width-2,1),DIRECTIONS,MPI_DOUBLE, up_right,5, comm->communicator,&comm->requests[req_count++] );
	MPI_Isend(lbm_mesh_get_cell(mesh, 1,comm->height-2), DIRECTIONS,MPI_DOUBLE,down_left,6,comm->communicator,&comm->requests[req_count++] );
	MPI_Isend( lbm_mesh_get_cell(mesh,comm->width-2,comm->height-2),DIRECTIONS, MPI_DOUBLE, down_right,7,comm->communicator,&comm->requests[req_count++] );
	MPI_Waitall(req_count, comm->requests, MPI_STATUSES_IGNORE);
	}