#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

double f(double x);

int main (int argc, char *argv[]) {

    int numprocs; /* Number of processors */
    int procnum; /* Processor number */

    /* Initialize MPI */
    MPI_Init(&argc, &argv);

    /* Find this processor number */
    MPI_Comm_rank(MPI_COMM_WORLD, &procnum);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* Master section */
    if ((procnum == 0) && (numprocs > 1)) {

        double a, b, c, h, fa, fb, fc, S, precision, partialSum, integral = 0;
        int i, workUnits, idleProcs, workers = numprocs - 1;
        double workPool[10000][4];
        MPI_Status status;

        /* idle procs array */
        int freeProcs[numprocs - 1];
        int start = 0, end = numprocs - 2;
        for (i = 0; i < numprocs - 1; i++) {
            freeProcs[i] = i + 1;
        }

        if ((argv[1] != NULL) && (argv[2] != NULL) && (argv[3] != NULL)) {
            a = atoi(argv[1]);
            b = atoi(argv[2]);
            precision = atof(argv[3]);
        } else {
            a = -249;
            b = 250;
            precision = 0.00000000000001;
        }

        c = (double) (a + b)/2;
        h = b - a;
        fa = f(a);
        fb = f(b);
        fc = f(c);
        S = ((double) (h/6))*(fa + 4*fc + fb);

        workPool[0][0] = a;
        workPool[0][1] = b;
        workPool[0][2] = precision;
        workPool[0][3] = S;

        workUnits = 1;
        idleProcs = numprocs - 1;

        double starttime, endtime;
        starttime = MPI_Wtime();

        while ((workUnits) || (idleProcs < workers)) {

            while ((workUnits) && (idleProcs)) {
                MPI_Send(workPool[--workUnits], 4, MPI_DOUBLE, freeProcs[start], 0, MPI_COMM_WORLD);
                start = (start + 1) % workers;
                idleProcs--;
            }

            MPI_Recv(workPool[workUnits], 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == 1) {
                MPI_Recv(workPool[++workUnits], 4, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                workUnits++;
            }
            idleProcs++;
            end = (end + 1) % workers;
            freeProcs[end] = status.MPI_SOURCE;

        }

        for (i = 0; i < workers; i++) {
            MPI_Send(workPool[workUnits], 4, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD);
        }

        for (i = 0; i < numprocs - 1; i++) {
            MPI_Recv(&partialSum, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            integral += partialSum;
        }

        endtime = MPI_Wtime();
        printf("interval [%f, %f], precision %E\n", a, b, precision);
        printf("procs: %d, ", numprocs);
        printf("integral: %.16f, ", integral);
        printf("running time %f seconds\n", endtime - starttime);

    /* Slave section */
    } else {

        /*
            workUnit[0] - a;
            workUnit[1] - b;
            workUnit[2] - epsilon;
            workUnit[3] - S;
        */

        double partialSum = 0;
        double workUnit[4];
        double workUnit1[4];
        double workUnit2[4];
        double c, h, d, e, fa, fb, fc, fd, fe, Sleft, Sright, S2;
        MPI_Status status;

        MPI_Recv(workUnit, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        while (status.MPI_TAG == 0) {

            c = (double) (workUnit[0] + workUnit[1])/2;
            h = workUnit[1] - workUnit[0];
            d = (double) (workUnit[0] + c)/2;
            e = (double) (c + workUnit[1])/2;
            fa = f(workUnit[0]);
            fb = f(workUnit[1]);
            fc = f(c);
            fd = f(d);
            fe = f(e);
            Sleft = ((double) (h/12))*(fa + 4*fd + fc);
            Sright = ((double) (h/12))*(fc + 4*fe + fb);
            S2 = Sleft + Sright;

            if (fabs(S2 - workUnit[3]) <= 15*workUnit[2]) {
                partialSum += S2;
                MPI_Send(workUnit, 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            } else if ((double) workUnit[2]/2 != 0) {
                workUnit1[0] = workUnit[0]; workUnit1[1] = c; workUnit1[2] = (double) workUnit[2]/2; workUnit1[3] = Sleft;
                workUnit2[0] = c; workUnit2[1] = workUnit[1]; workUnit2[2] = (double) workUnit[2]/2; workUnit2[3] = Sright;
                MPI_Send(workUnit1, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Send(workUnit2, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            }
            MPI_Recv(workUnit, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        }
        MPI_Send(&partialSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}

double f(double x) {

    return sin(x);

}
