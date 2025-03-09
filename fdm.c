/*****************************************************************
/* Author: Pablo Goitia <pablo.goitia@alumnos.unican.es>
/* Project: Finite-Difference Approximation to the Heat Equation
/* Date: Feb-2025
/*
/* Usage: ./fdm [deltaH]
/*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

unsigned int mdf_heat(double ***__restrict__ u0,
					  double ***__restrict__ u1,
					  const unsigned int npX,
					  const unsigned int npY,
					  const unsigned int npZ,
					  const double deltaH,
					  const double deltaT,
					  const double inErr,
					  const double boundaries,
					  const unsigned int first_point,
					  const unsigned int last_point,
					  const int myrank,
					  const int size)
{

	register double alpha = deltaT / (deltaH * deltaH);
	int continued = 1;
	register unsigned int steps = 0;

	// Loop indexes
	unsigned int i, j, k, i_internal;
	unsigned int points_per_slice = last_point - first_point + 1;

	// Boundaries
	register double left;
	register double right;
	register double up;
	register double down;
	register double top;
	register double bottom;

	// MPI Status variable
	MPI_Status status;

	// MPI Pack and Unpack variables
	int position = 0;

	// Buffer for MPI communication
	double *buffer = (double *)malloc(npY * npZ * sizeof(double));
	
	while (continued)
	{
		steps++;

		// Calculate the provisional values for the points in the slice, excluding 'left' and 'right' points
		for (i = 1; i <= points_per_slice; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = 0; k < npZ; k++)
				{
					up = boundaries;
					down = boundaries;
					top = boundaries;
					bottom = boundaries;

					if ((j > 0) && (j < (npY - 1)))
					{
						up = u0[i][j - 1][k];
						down = u0[i][j + 1][k];
					}
					else if (j == 0)
						down = u0[i][j + 1][k];
					else
						up = u0[i][j - 1][k];

					if ((k > 0) && (k < (npZ - 1)))
					{
						top = u0[i][j][k - 1];
						bottom = u0[i][j][k + 1];
					}
					else if (k == 0)
						bottom = u0[i][j][k + 1];
					else
						top = u0[i][j][k - 1];

					// Calculates a provisional value for u1[i][j][k]
					u1[i][j][k] = alpha * (top + bottom + up + down - (6.0f * u0[i][j][k])) + u0[i][j][k];
				}
			}
		}

		// Exchange first and last point values with neighbour MPI processes
		position = 0;

		if (myrank > 0)
		{
			// Send the content of the first point of the slice to the left neighbour
			for (i = 0; i < npY; i++)
				MPI_Pack(&u0[1][i][0], npZ, MPI_DOUBLE, buffer, npY * npZ * sizeof(double), &position, MPI_COMM_WORLD);

			MPI_Send(buffer, npY * npZ, MPI_DOUBLE, myrank - 1, steps, MPI_COMM_WORLD);
		}

		position = 0;
		
		if (myrank < (size - 1))
		{
			// Send the content of the last point of the slice to the right neighbour
			for (i = 0; i < npY; i++)
				MPI_Pack(&u0[points_per_slice][i][0], npZ, MPI_DOUBLE, buffer, npY * npZ * sizeof(double), &position, MPI_COMM_WORLD);

			MPI_Send(buffer, npY * npZ, MPI_DOUBLE, myrank + 1, steps, MPI_COMM_WORLD);
		}

		position = 0;

		if (myrank > 0)
		{
			// Receive the content of the first point of the slice from the left neighbour
			MPI_Recv(buffer, npY * npZ, MPI_DOUBLE, myrank - 1, steps, MPI_COMM_WORLD, &status);
			for (i = 0; i < npY; i++)
				MPI_Unpack(buffer, npY * npZ * sizeof(double), &position, &u0[0][i][0], npZ, MPI_DOUBLE, MPI_COMM_WORLD);
		}

		position = 0;
		
		if (myrank < (size - 1))
		{
			// Receive the content of the last point of the slice from the right neighbour
			MPI_Recv(buffer, npY * npZ, MPI_DOUBLE, myrank + 1, steps, MPI_COMM_WORLD, &status);
			for (i = 0; i < npY; i++)
				MPI_Unpack(buffer, npY * npZ * sizeof(double), &position, &u0[points_per_slice + 1][i][0], npZ, MPI_DOUBLE, MPI_COMM_WORLD);
		}

		i_internal = 1;	// Internal index for the slice, which also avoids the 'left' point

		// Calculate the definitive values for the points in the slice, including 'left' and 'right' points
		for (i = first_point; i <= last_point; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = 0; k < npZ; k++)
				{
					left = boundaries;
					right = boundaries;

					// Checks if the point is not in the borders of the TOTAL X axis
					if ((i > 0) && (i < (npX - 1)))
					{
						left = u0[i_internal - 1][j][k];
						right = u0[i_internal + 1][j][k];
					}
					else if (i == 0)
						right = u0[i_internal + 1][j][k];
					else
						left = u0[i_internal - 1][j][k];

					u1[i_internal][j][k] += alpha * (left + right);
				}
			}
			i_internal++;
		}

		double ***ptr = u0;
		u0 = u1;
		u1 = ptr;

		double err = 0.0f;
		double maxErr = 0.0f;
		for (i = 1; i <= points_per_slice; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = 0; k < npZ; k++)
				{
					err = fabs(u0[i][j][k] - boundaries);
					if (err > inErr)
                      	maxErr = err;
					else
						continued = 0;
				}
			}
		}
		
		// if (myrank == 0)
		// 	printf ("err = %.4g > inErr = %.4g\n", err, inErr);

		// Use MPI_LAND to check if any process has set 'continued' to 0
		MPI_Allreduce(&continued, &continued, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
		/** Better idea: to signal the other processes when a 0 has been found in one of them. */
	}

	return steps;
}

int main(int ac, char **av)
{
	double ***u0;
	double ***u1;

	double deltaT = 0.01;
	double deltaH = 0.005f;
	double sizeX = 1.0f;
	double sizeY = 1.0f;
	double sizeZ = 1.0f;

	// MPI variables
	int myrank, size;
	MPI_Status status;

	// Start MPI communication
	MPI_Init(&ac, &av);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Input processing (common for all MPI processes)
	if (ac > 1)
	{
		deltaH = atof(av[1]);
		if (deltaH <= 0.0f || deltaH >= 1.0f)
		{
			if (myrank == 0)
				fprintf(stderr, "Error: deltaH must be (0.0, 1.0)\n");
			exit(EXIT_FAILURE);
		}
	}

	// Calculate the number of points in each axis
	unsigned int npX = (unsigned int)(sizeX / deltaH);
	unsigned int npY = (unsigned int)(sizeY / deltaH);
	unsigned int npZ = (unsigned int)(sizeZ / deltaH);

	/** The three-dimensional space will be divided into a number of slices across the X axis
	 * matching the number of MPI processes, and which are made up of a set of points. */
	unsigned int points_per_slice = npX / size;

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
	unsigned int remaining_points = npX % size;

	/** Calculate the number of points for each MPI process. Processes with 
	 * myrank < remaining_points will have one more point. */
	if (myrank < remaining_points)
		points_per_slice++;

	/** Calculate the first and last point in the slice for each MPI process, 
	 * considering that some processes will have one more point. */
	unsigned int first_point = 0;
	for (int i = 0; i < myrank; i++)
		first_point += (i < remaining_points) ? (points_per_slice + 1) : points_per_slice;
	unsigned int last_point = first_point + points_per_slice - 1;

	// TODO: Idea:
	// To avoid the subsequent execution for processes with 'points_per_slice' = 0, i.e. when num_procs > npX.

	// Allocating memory for the tri-dimensional space
	// Memory allocation for X axis
	u0 = (double ***)malloc((points_per_slice + 2) * sizeof(double **));
	u1 = (double ***)malloc((points_per_slice + 2) * sizeof(double **));

	// Processes print the number of points in each axis
	printf("I am process %d of %d. p(%u, %u, %u). points_per_slice=%d. first_point=%d, last_point=%d\n", myrank, size, npX, npY, npZ, points_per_slice, first_point, last_point);

	// Memory allocation for Y axis
	for (unsigned int i = 0; i < (points_per_slice + 2); i++)
	{
		u0[i] = (double **)malloc(npY * sizeof(double *));
		u1[i] = (double **)malloc(npY * sizeof(double *));
	}

	/** Memory allocation for Z axis.
	 * In this case, we will allocate memory for the points in each slice, that is,
	 * the value of points_per_slice.
	 * We will also allocate memory for the two contiguous points in the X axis, because
	 * the algorithm will need to access the values of that points in the neighbour slices.
	 */
	for (unsigned int i = 0; i < (points_per_slice + 2); i++)
	{
		for (unsigned int j = 0; j < npY; j++)
		{
			double *aux0 = (double *)malloc(npZ * sizeof(double));
			double *aux1 = (double *)malloc(npZ * sizeof(double));
			// initial condition - zero in all points
			memset(aux0, 0x01, npZ * sizeof(double));
			memset(aux1, 0x02, npZ * sizeof(double));
			u0[i][j] = aux0;
			u1[i][j] = aux1;
		}
	}

	// Each MPI process will compute its own points with the finite difference method
	unsigned int steps, max_steps;
	steps = mdf_heat(u0, u1, npX, npY, npZ, deltaH, deltaT, 1e-15, 100.0f, first_point, last_point, myrank, size);

	// Collect the number of steps from all MPI processes and get the maximum value
	MPI_Reduce(&steps, &max_steps, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);

	if (myrank == 0)
		fprintf(stdout, "Done! in %u steps\n", max_steps);

	// Free memory
	for (unsigned int i = 0; i < (points_per_slice + 2); i++)
	{
		for (unsigned int j = 0; j < npY; j++)
		{
			free(u0[i][j]);
			free(u1[i][j]);
		}
	}

	for (unsigned int i = 0; i < (points_per_slice + 2); i++)
	{
		free(u0[i]);
		free(u1[i]);
	}

	free(u0);
	free(u1);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
