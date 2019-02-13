/* =============================================
 * Kostas Evangelou
 * mpi_heat2D_opt.c
 * Parallel Systems
 * MPI Optimization of heat2D problem
 * =============================================*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define  CONVERGENCE 0                 /* Enable Convergence */ 
#define  NXPROB     20                 /* x dimension of problem grid */
#define  NYPROB     20                 /* y dimension of problem grid */
#define  STEPS     100                 /* number of time steps */
#define  DTAG        0                 /* message tag */
#define  UTAG        1                 /* message tag */
#define  LTAG        2                 /* message tag */
#define  RTAG        3                 /* message tag */
#define  NONE       -1                 /* indicates no neighbor */
#define  DONE        4                 /* message tag */
#define  MASTER      0                 /* taskid of first process */
#define  UP          0
#define  DOWN        1
#define  LEFT        2
#define  RIGHT       3

struct Parms {
    float cx;
    float cy;
} parms = {0.1, 0.1};

int main(int argc, char *argv[]) {
    void    inidat(), prtdat(), update();
    float   ***u;                       /* array for grid */
    float   *temp[2];                   /* temp 1d array for fast memory allocation */
    int     taskid,                     /* this task's unique id */
            numtasks,                   /* number of worker processes */
            rc=0,                       /* misc */
            i, j, ix, iy, iz,           /* loop variables */
            dims[2],                    /* cart argument*/
            reorder = 0,                /* Ranking may be reordered (true) or not (false) (logical) */
            periods[2],                 /* Logical  array of size ndims specifying whether the grid is periodic */
            neighbors[4],               /* Neighbors indicator */
            size,                       /* Calculate num of grid */
            sub_y,sub_x,
            sub_table_dim,
            step;
    MPI_Status status[8];
    double start_time = 0.0,
            end_time = 0.0,
            task_time = 0.0,
            reduced_time=0.0;

    MPI_Datatype MPI_row,
                 MPI_column;

/* First, find out my taskid and how many tasks are running */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if ((NYPROB * NXPROB) % numtasks != 0) {
        printf("ERROR: Number of tasks must have remainder 0 when divided with number of workers %d x %d / %d \n", NXPROB, NYPROB, numtasks);
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    /* Mpi Cartesian Grid */

    MPI_Comm cart_comm;
    size = sqrt(numtasks);
    dims[0] = size;
    dims[1] = size;
    periods[0]=periods[1]=0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);

    MPI_Barrier(cart_comm);
    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[UP], &neighbors[DOWN]);
    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[LEFT], &neighbors[RIGHT]);

    if (taskid==0){
        printf ("Up:    %d\n",neighbors[UP]);
        printf ("Down:  %d\n",neighbors[DOWN]);
        printf ("Left:  %d\n",neighbors[LEFT]);
        printf ("Right: %d\n",neighbors[RIGHT]);
    }
    sleep(10);
    //MPI_Abort(MPI_COMM_WORLD, rc);
    //exit(1);
    sub_table_dim = sqrt(NXPROB*NYPROB/numtasks);
    sub_x = sub_table_dim + 2;
    sub_y = sub_table_dim + 2;

    /* Declare MPI Types */
    MPI_Type_contiguous(sub_x, MPI_FLOAT, &MPI_row);
    MPI_Type_commit(&MPI_row);
    MPI_Type_vector(sub_y, 1, sub_x, MPI_FLOAT, &MPI_column);
    MPI_Type_commit(&MPI_column);

    MPI_Request req[8];

    /* Allocate Memory for Table */
    u=(float***)malloc(2 * sizeof(float**));
    for(i=0;i<2;i++){
      u[i] = (float**)malloc(sub_x * sizeof(float*));
      temp[i] = (float*)malloc(sub_y*sub_x*sizeof(float));
      for(j=0;j<sub_x;j++){
        u[i][j] = temp[i] + (j*sub_y);
      }
    }

    /* Init table with 0s */
    for (iz=0; iz<2; iz++)
        for (ix=0; ix<sub_x; ix++)
            for (iy=0; iy<sub_y; iy++)
                u[iz][ix][iy] = 0.0;

    /* Set Random values to sub_table */
    inidat(sub_x, sub_y, taskid, u[0]);

    MPI_Barrier(cart_comm);
    start_time = MPI_Wtime();

#if (CONVERGENCE == 1)
    task_convergence = reduced_convergence = 0;
#endif

    for(step=0; step<=STEPS;step++){
        /* Send and Receive asynchronous the shared values of neighbor */
        if (neighbors[UP] >= 0){
            MPI_Isend(u + iz*sub_x*sub_y + sub_y + 1, sub_table_dim, MPI_FLOAT, neighbors[UP], DTAG, cart_comm, &req[0]);
            MPI_Irecv(u + iz*sub_x*sub_y + 1, sub_table_dim, MPI_FLOAT, neighbors[UP], UTAG, cart_comm, &req[1]);
        }

        if (neighbors[DOWN] >= 0){
            MPI_Isend(u + iz*sub_x*sub_y + sub_table_dim*sub_y + 1, sub_table_dim , MPI_FLOAT, neighbors[DOWN], UTAG, cart_comm, &req[2]);
            MPI_Irecv(u + iz*sub_x*sub_y + (sub_table_dim+1)*sub_y + 1, sub_table_dim , MPI_FLOAT, neighbors[DOWN], DTAG, cart_comm, &req[3]);
        }

        if (neighbors[LEFT] >= 0){
            MPI_Isend(u + iz*sub_x*sub_y + sub_y + 1, 1, MPI_column, neighbors[LEFT], RTAG, cart_comm,&req[4]);
            MPI_Irecv(u + iz*sub_x*sub_y + sub_y, 1, MPI_column, neighbors[LEFT], LTAG, cart_comm, &req[5]);
        }

        if (neighbors[RIGHT] >= 0 ){
            MPI_Isend(u + iz*sub_x*sub_y + sub_y + sub_table_dim, 1, MPI_column, neighbors[RIGHT], LTAG, cart_comm,&req[6]);
            MPI_Irecv(u + iz*sub_x*sub_y + sub_y + sub_table_dim + 1, 1, MPI_column , neighbors[RIGHT], RTAG, cart_comm,&req[7]);
        }

        /* Update inside table while the process wait for neighbor values */
        update(sub_table_dim, u + iz*sub_x*sub_y, u + (1-iz)*sub_x*sub_y);

        /* Wait for neighbor values */
        if(neighbors[UP] >= 0){
            MPI_Wait(&req[0],&status[0]);
            MPI_Wait(&req[1],&status[1]);
        }
        if(neighbors[DOWN] >= 0){
            MPI_Wait(&req[2],&status[2]);
            MPI_Wait(&req[3],&status[3]);
        }
        if(neighbors[LEFT] >= 0){
            MPI_Wait(&req[4],&status[4]);
            MPI_Wait(&req[5],&status[5]);
        }
        if(neighbors[RIGHT] >= 0){
            MPI_Wait(&req[6],&status[6]);
            MPI_Wait(&req[7],&status[7]);
        }

        /* Update outside table with neighboor values */
        update(sub_table_dim, u + iz*sub_x*sub_y, u + (1-iz)*sub_x*sub_y);

        /* Next loop with have to deal with the other table */
        iz = 1 - iz;

#if (CONVERGENCE == 1)
        if (step % CONVERGENCE_n == 0) {
            task_convergence = check_convergence(1, task_X - 2, 1, task_Y - 2, task_Y, *u[iz], *u[1 - iz]);
            MPI_Barrier(cart_comm);
            MPI_Allreduce(&task_convergence, &reduced_convergence, 1, MPI_INT, MPI_LAND, cart_comm);
            if (reduced_convergence == 1) {
                break;
            }
        }
#endif
    }
    /* Stop timer, calculate duration, reduce */
    end_time = MPI_Wtime();
    task_time = end_time - start_time;
    MPI_Barrier(cart_comm);
    MPI_Reduce(&task_time, &reduced_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    if (taskid == 0) {
        printf("Convergence: %d\n", CONVERGENCE);
        printf("u size: [%d][%d]\n", NXPROB, NYPROB);
        printf("tasks: %d\n", numtasks);
        printf("Time elapsed: %f seconds\n", reduced_time);
    }
    /****************/

    /* Cleanup everything */
    free(temp[0]);
    free(temp[1]);
    free(u[0]);
    free(u[1]);
    free(u);
    MPI_Type_free(&MPI_row);
    MPI_Type_free(&MPI_column);
    MPI_Finalize();
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update ( int start_x, int end_x,int start_y, int end_y,int ny, float *u1, float *u2 )
{
    int ix, iy;
    for ( ix = start_x; ix <= end_x; ix++ )
        for ( iy = start_y; iy <= end_y; iy++ )
            * ( u2+ix*ny+iy ) = * ( u1+ix*ny+iy )  +
                                parms.cx * ( * ( u1+ ( ix+1 ) *ny+iy ) +
                                             * ( u1+ ( ix-1 ) *ny+iy ) -
                                             2.0 * * ( u1+ix*ny+iy ) ) +
                                parms.cy * ( * ( u1+ix*ny+iy+1 ) +
                                             * ( u1+ix*ny+iy-1 ) -
                                             2.0 * * ( u1+ix*ny+iy ) );
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, int id, float **u) {
int ix, iy;
int r = rand()%10;
for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     u[ix][iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1))*(r+id);
}