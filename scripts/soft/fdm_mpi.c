/*****************************************************************
 * Author: Pablo Goitia <pablo.goitia@alumnos.unican.es>
 * Project: Finite-Difference Approximation to the Heat Equation
 * Version: MPI
 * Date: Mar-2025
 *
 * Usage: ./fdm [deltaH]
 * Compile: mpicc -o fdm_mpi fdm_mpi.c
 *****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

unsigned int mdf_heat(double *__restrict__ u0,
					  double *__restrict__ u1,
					  const unsigned int npX,
					  const unsigned int npY,
					  const unsigned int npZ,
					  const double deltaH,
					  const double deltaT,
					  const double inErr,
					  const double boundaries,
					  MPI_Comm compute_comm,
					  const int myrank,
					  const int size,
					  const unsigned int first_point,
					  const unsigned int last_point)
{
	register double alpha = deltaT / (deltaH * deltaH);
	int continued = 1;
	register unsigned int steps = 0;

	// Loop indexes
	unsigned int i, j, k, idx;
	unsigned int points_per_slice = last_point - first_point + 1;

	// Boundaries
	register double left;
	register double right;
	register double up;
	register double down;
	register double top;
	register double bottom;

	double *temp;
	double err, maxErr;

	// MPI comms handling variables
	MPI_Status status;
	MPI_Request request_left, request_right;

	while (1)
	{
		steps++;

		err = 0.0f;
		maxErr = 0.0f;

		for (i = first_point; i <= last_point; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = 0; k < npZ; k++)
				{
					left = boundaries;
					right = boundaries;
					up = boundaries;
					down = boundaries;
					top = boundaries;
					bottom = boundaries;

					idx = (i - first_point + 1) * npY * npZ + j * npZ + k;

					if ((i > 0) && (i < (npX - 1)))
					{
						left = u0[idx - npY * npZ];
						right = u0[idx + npY * npZ];
					}
					else if (i == 0)
						right = u0[idx + npY * npZ];
					else
						left = u0[idx - npY * npZ];

					if ((j > 0) && (j < (npY - 1)))
					{
						up = u0[idx - npZ];
						down = u0[idx + npZ];
					}
					else if (j == 0)
						down = u0[idx + npZ];
					else
						up = u0[idx - npZ];

					if ((k > 0) && (k < (npZ - 1)))
					{
						top = u0[idx - 1];
						bottom = u0[idx + 1];
					}
					else if (k == 0)
						bottom = u0[idx + 1];
					else
						top = u0[idx - 1];

					u1[idx] = alpha * (top + bottom + up + down + left + right - (6.0f * u0[idx])) + u0[idx];

					err = fabs(u1[idx] - boundaries);
					if (err > inErr)
						maxErr = err;
					else
						continued = 0;
				}
			}
		}

		// Use MPI_LAND to check if any process has set 'continued' to 0
		MPI_Allreduce(&continued, &continued, 1, MPI_INT, MPI_LAND, compute_comm);

		if (continued == 0)
			break;

		// Exchange first and last point values with neighbour MPI processes
		if (myrank > 0)
		{
			// Send the content of the first point of the slice to the left neighbour
			MPI_Isend(&u1[npZ * npY], npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &request_left);

			// Receive the content of the first point of the slice from the left neighbour
			MPI_Irecv(&u1[0], npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &request_left);
		}

		if (myrank < (size - 1))
		{
			// Send the content of the last point of the slice to the right neighbour
			MPI_Isend(&u1[npZ * npY * points_per_slice], npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &request_right);

			// Receive the content of the last point of the slice from the right neighbour
			MPI_Irecv(&u1[npZ * npY * (points_per_slice + 1)], npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &request_right);
		}

		// Check if the messages have been received
		if (myrank > 0)
			MPI_Wait(&request_left, &status);

		if (myrank < (size - 1))
			MPI_Wait(&request_right, &status);

		// Swap the pointers (the next instant of time is now the current time)
		temp = u0;
		u0 = u1;
		u1 = temp;
	}

	return steps;
}

int main(int ac, char **av)
{
	double start_main, end_main, execution_main;
	start_main = MPI_Wtime();

	// Data structures
	double *u0;
	double *u1;

	// Loop indexes
	unsigned int i;

	// Simulation parameters
	double deltaT = 0.01;
	double deltaH = 0.005f;
	double sizeX = 1.0f;
	double sizeY = 1.0f;
	double sizeZ = 1.0f;
	unsigned int npX, npY, npZ;

	// MPI variables
	int my_global_rank;

	// Subcommunicator for the MPI processes that will compute the points
	MPI_Comm compute_comm;
	int myrank, size;
	int is_computing_proc = 0;

	// Variables for the MPI processes
	unsigned int points_per_slice, remaining_points;
	unsigned int first_point, last_point;

	// Variables for results
	unsigned int steps, max_steps;
	double start_time, end_time, execution_time, max_time;

	// Start MPI communication
	MPI_Init(&ac, &av);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// Input processing
	if (ac > 1)
	{
		deltaH = atof(av[1]);
		if (deltaH <= 0.0f || deltaH >= 1.0f)
		{
			if (my_global_rank == 0)
				fprintf(stderr, "Error: deltaH must be a value in the range (0.0, 1.0)\n");
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

	// Calculate the number of points in each axis
	npX = (unsigned int)(sizeX / deltaH);
	npY = (unsigned int)(sizeY / deltaH);
	npZ = (unsigned int)(sizeZ / deltaH);

	// Discriminate the MPI processes that will compute the points and the remaining ones
	is_computing_proc = (my_global_rank < npX) ? 1 : 0;
	MPI_Comm_split(MPI_COMM_WORLD, is_computing_proc, my_global_rank, &compute_comm);
	MPI_Comm_rank(compute_comm, &myrank);
	MPI_Comm_size(compute_comm, &size);

	if (is_computing_proc)
	{
		/** The three-dimensional space will be divided into a number of slices across the X axis
		 * matching the number of MPI processes, and which are made up of a set of points. */
		points_per_slice = npX / size; // Provisional value

		/** Some slices might have more points than others if npX cannot be exactly divided by
		 * the number of MPI processes, so we will implement the simplest approach for a load
		 * balancing algorithm:
		 *  	Distribute the remaining points among the first MPI processes.
		 *  	Assuming equal workload for each point and equal hardware resources, this would
		 *  	be the fairest algorithm.
		 *
		 *  Obviously, those conditions are not frequently met, so this algorithm is not
		 *  the best one, but it will be enough for this first implementation.
		 */
		remaining_points = npX % size;

		// Calculate the definitive number of points for each MPI process
		if (myrank < remaining_points)
			points_per_slice++;

		/** Calculate the first and last point in the slice for each MPI process,
		 * considering that some processes will have one more point. */
		first_point = 0;
		for (i = 0; i < myrank; i++)
			first_point += (i < npX % size) ? (points_per_slice + 1) : points_per_slice;
		last_point = first_point + points_per_slice - 1;

		// Processes print the number of points in each axis
		printf("I am process %d of %d. points_per_slice=%d of %d. first_point=%d, last_point=%d\n",
			   myrank + 1, size, points_per_slice, npX, first_point, last_point);

		/** Allocating memory for the tri-dimensional space, including:
		 * - The points in each slice, that is, the value of points_per_slice.
		 * - Two more points in the X axis, because the algorithm will need to access the
		 *   left and right neighbour points.
		 */
		u0 = malloc((points_per_slice + 2) * npY * npZ * sizeof(double *));
		u1 = malloc((points_per_slice + 2) * npY * npZ * sizeof(double *));

		// Initial condition: zero in all points
		memset(u0, 0x01, (points_per_slice + 2) * npY * npZ * sizeof(double *));
		memset(u1, 0x02, (points_per_slice + 2) * npY * npZ * sizeof(double *));

		start_time = MPI_Wtime();

		// Each MPI process will compute its assigned points with the finite difference method
		steps = mdf_heat(u0, u1, npX, npY, npZ, deltaH, deltaT, 1e-15, 100.0f, compute_comm, myrank, size, first_point, last_point);

		// Collect the number of steps from all MPI processes and get the maximum value
		MPI_Reduce(&steps, &max_steps, 1, MPI_UNSIGNED, MPI_MAX, 0, compute_comm);

		// Save the execution time of the parallel region
		end_time = MPI_Wtime();
		execution_time = end_time - start_time;

		// Get max execution time across all processes
		MPI_Reduce(&execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, compute_comm);

		// Print the results
		if (myrank == 0)
			printf("Done! in %u steps\nExecution time: %f seconds\n", max_steps, max_time);

		// Free memory
		free(u0);
		free(u1);
	}

	// Free the subcommunicators
	MPI_Comm_free(&compute_comm);

	MPI_Finalize();

	end_main = MPI_Wtime();
	execution_main = end_main - start_main;
	printf("Main execution time: %f seconds\n", execution_main);

	return EXIT_SUCCESS;
}
