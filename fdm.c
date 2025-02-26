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
					  const unsigned int last_point)
{

	register double alpha = deltaT / (deltaH * deltaH);
	register int continued = 1;
	register unsigned int steps = 0;

	// Loop indexes
	unsigned int i, j, k;
	unsigned int points_per_slice = last_point - first_point + 1;

	while (continued)
	{
		steps++;
		for (i = 0; i < npZ; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = first_point; k < last_point; k++)
				{
					register double left = boundaries;
					register double right = boundaries;
					register double up = boundaries;
					register double down = boundaries;
					register double top = boundaries;
					register double bottom = boundaries;
					
					// k point is inside the X axis, and not in the borders
					if ((k > 0) && (k < (npX - 1)))
					{
						// k needs to read information from the left slide last point
						// [...]
						// k needs to read information from the right slide first point
						// [...]

						// Pre-existing code
						left = u0[i][j][k - 1];
						right = u0[i][j][k + 1];
					}
					// k is at the left border of the X axis (k = 0)
					else if (k == 0)
						// left = boundaries;
						right = u0[i][j][k + 1];
					// k is at the right border of the X axis (k = npX - 1)
					else
						left = u0[i][j][k - 1];
						// right = boundaries;

					if ((j > 0) && (j < (npY - 1)))
					{
						up = u0[i][j - 1][k];
						down = u0[i][j + 1][k];
					}
					else if (j == 0)
						down = u0[i][j + 1][k];
					else
						up = u0[i][j - 1][k];

					if ((i > 0) && (i < (npZ - 1)))
					{
						top = u0[i - 1][j][k];
						bottom = u0[i + 1][j][k];
					}
					else if (i == 0)
						bottom = u0[i + 1][j][k];
					else
						top = u0[i - 1][j][k];

					// u1[i][j][k] = alpha * (top + bottom + up + down + left + right - (6.0f * u0[i][j][k])) + u0[i][j][k];
					/** Calculates a provisional value for u1[i][j][k] until left and right values are received 
					 * from the neighbour MPI processes. */
					u1[i][j][k] = alpha * (top + bottom + up + down - (6.0f * u0[i][j][k])) + u0[i][j][k];
					// printf("u1[%d][%d][%d] = %f\n", i, j, k, u1[i][j][k]);
				}
			}
		}

		double ***ptr = u0;
		u0 = u1;
		u1 = ptr;

		double err = 0.0f;
		double maxErr = 0.0f;
		for (i = 0; i < npZ; i++)
		{
			for (j = 0; j < npY; j++)
			{
				for (k = 0; k < npX; k++)
				{
					err = fabs(u0[i][j][k] - boundaries);
					if (err > inErr)
                    {
                      if (err != maxErr)
                         printf ("err = %.15g > inErr = %.15g\n", err, inErr);
                      maxErr = err;
                    }
					else
						continued = 0;
				}
			}
		}
		// printf("err = %f\n", err);
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

	// Allocating memory for the tri-dimensional space
	// Memory allocation for Z axis
	u0 = (double ***)malloc(npZ * sizeof(double **));
	u1 = (double ***)malloc(npZ * sizeof(double **));

	// Memory allocation for Y axis
	for (unsigned int i = 0; i < npZ; i++)
	{
		u0[i] = (double **)malloc(npY * sizeof(double *));
		u1[i] = (double **)malloc(npY * sizeof(double *));
	}

	/** Memory allocation for X axis.
	 * In this case, we will allocate memory for the points in each slice, that is,
	 * the value of points_per_slice.
	 * We will also allocate memory for the two contiguous points in the X axis, because
	 * the algorithm will need to access the values of that points in the neighbour slices.
	 */
	for (unsigned int i = 0; i < npZ; i++)
	{
		for (unsigned int j = 0; j < npY; j++)
		{
			double *aux0 = (double *)malloc(points_per_slice * sizeof(double));
			double *aux1 = (double *)malloc(points_per_slice * sizeof(double));
			// initial condition - zero in all points
			memset(aux0, 0x01, (points_per_slice + 2) * sizeof(double));
			memset(aux1, 0x02, (points_per_slice + 2) * sizeof(double));
			u0[i][j] = aux0;
			u1[i][j] = aux1;
		}
	}

	// Each MPI process will compute its own points with the finite difference method
	unsigned int steps = mdf_heat(u0, u1, npX, npY, npZ, deltaH, deltaT, 1e-15, 100.0f, first_point, last_point);
	unsigned int total_steps;

	// Collect the number of steps from all MPI processes and sum them
	MPI_Reduce(&steps, &total_steps, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	fprintf(stdout, "Done! in %u steps\n", total_steps);

	// Free memory
	for (unsigned int i = 0; i < npZ; i++)
	{
		for (unsigned int j = 0; j < npY; j++)
		{
			free(u0[i][j]);
			free(u1[i][j]);
		}
	}

	for (unsigned int i = 0; i < npZ; i++)
	{
		free(u0[i]);
		free(u1[i]);
	}

	free(u0);
	free(u1);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
