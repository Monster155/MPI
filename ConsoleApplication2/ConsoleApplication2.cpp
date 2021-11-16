#include <algorithm>
#include <cstdio>
#include "mpi.h"

namespace Task1
{
    int main(int argc, char** argv)
    {
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        printf("Hello from %d thread", rank);
        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

namespace Task2
{
    int main(int argc, char** argv)
    {
        int arr[] = {15, 234, 765, 4, 564, 823, 45, 24, 312, 21, 34, 165, 34, 6};
        constexpr int arr_size = 14;
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        const int start = arr_size / size * rank + rank;
        int end = arr_size / size * (rank + 1) + rank;
        if (end >= arr_size) end = arr_size - 2;

        printf("start %d end %d\n", start, end);

        int max = INT_MIN;
        for (int i = start; i <= end; i++)
            if (max < arr[i])
                max = arr[i];

        printf("from %d max %d\n", rank, max);

        if (rank)
            MPI_Send(&max, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        else
        {
            for (int i = 0; i < size - 1; ++i)
            {
                int max2;
                MPI_Status st;
                MPI_Recv(&max2, 1,MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                if (max < max2)
                    max = max2;
            }
            printf("max: %d", max);
        }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

namespace Task3
{
    int main(int argc, char** argv)
    {
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        double dotsMax = 1e6;
        const double rand_max = RAND_MAX;
        int dotsInside = 0;

        srand(rank);
        for (int i = 0; i < dotsMax; ++i)
        {
            double x = rand() / rand_max;
            double y = rand() / rand_max;
            // printf("x = %f y = %f", x, y);

            if (1 > x * x + y * y)
                dotsInside++;
        }

        printf("in %d thread %d dots from %.f is inside\n", rank, dotsInside, dotsMax);
        double pi = dotsInside / dotsMax * 4;

        printf("Pi in %d thread = %f\n", rank, pi);

        if (rank)
            MPI_Send(&pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        else
        {
            double piSum = pi;
            for (int i = 0; i < size - 1; ++i)
            {
                double pi2;
                MPI_Status st;
                MPI_Recv(&pi2, 1,MPI_DOUBLE, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                piSum += pi2;
            }

            printf("pi: %f", piSum / size);
        }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

namespace Task4
{
    int main(int argc, char** argv)
    {
        int arr[] = {15, 234, 765, 4, 564, 823, 45, 24, 312, 21, 34, 165, 34, 6};
        int sC[] = {4, 4, 3, 3};
        int sD[] = {0, 4, 8, 11};
        constexpr int arr_size = 14;
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        const int rsize = sC[rank];

        int rbuf[arr_size];
        MPI_Scatterv(arr, sC, sD, MPI_INT, &rbuf, arr_size, MPI_INT, 0, MPI_COMM_WORLD);

        // printf("%d : ", rank);
        // for (int i = 0; i < rsize; ++i)
        //     printf("%d ", rbuf[i]);

        int sum = 0;
        for (int i = 0; i < rsize; i++)
            if (rbuf[i] > 0)
                sum += rbuf[i];
        printf("sum in %d thread = %d\n", rank, sum);

        if (rank)
            MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        else
        {
            double average = sum;
            for (int i = 0; i < size - 1; ++i)
            {
                int sum2;
                MPI_Status st;
                MPI_Recv(&sum2, 1,MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                average += sum2;
            }
            printf("average: %f", average / arr_size);
        }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

namespace Task5
{
    int main(int argc, char** argv)
    {
        int arr[]{15, 234, 765, 4, 564, 823, 45, 24, 312, 21, 34, 165, 34, 6};
        int brr[]{45, 765, 823, 21, 34, 312, 165, 34, 24, 4, 564, 15, 234, 9};
        const int n = 14;
        int size, rank;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int* displs = (int*)malloc(size);
        int* scounts = (int*)malloc(size);
        int* rAbuf = (int*)malloc(n);
        int* rBbuf = (int*)malloc(n);

        for (int i = 0; i < size; i++)
        {
            if (i + 1 == size)
                scounts[i] = n - (i * (n / size));
            else
                scounts[i] = n / size;

            displs[i] = i * (n / size);
        }

        MPI_Scatterv(arr, scounts, displs, MPI_INT, rAbuf, n, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(brr, scounts, displs, MPI_INT, rBbuf, n, MPI_INT, 0, MPI_COMM_WORLD);

        int multiplySum = 0;
        for (int i = 0; i < scounts[rank]; i++)
            multiplySum += rAbuf[i] * rBbuf[i];

        if (rank)
            MPI_Send(&multiplySum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        else
        {
            for (int i = 1; i < size; ++i)
            {
                int ms;
                MPI_Status st;
                MPI_Recv(&ms, 1,MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                multiplySum += ms;
            }
            printf("scalar: %d", multiplySum);
        }

        MPI_Finalize();
        return MPI_SUCCESS;
    }
}

namespace Task6
{
    int main(int argc, char** argv)
    {
        const int n = 10, m = 10;
        int size, rank;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Status status;

        int* matrix = (int*)malloc(n * m);
        int* rbuf = (int*)malloc(n * m);
        int* displs = (int*)malloc(size);
        int* scounts = (int*)malloc(size);

        if (!rank)
        {
            printf("Matix:\n");
            for (int i = 0; i < n * m; i++)
            {
                matrix[i] = rand();
                printf("%d ", matrix[i]);
                if ((i + 1) % m == 0)
                    printf("\n");
            }
            printf("\n");
        }

        for (int i = 0; i < size; i++)
        {
            if (i + 1 == size)
                scounts[i] = m * (n - (i * (n / size)));
            else
                scounts[i] = n * m / size;

            displs[i] = m * n * i / size;
        }

        // printf("%d %d", rank, scounts[rank]);

        MPI_Scatterv(matrix, scounts, displs, MPI_INT, rbuf, n * m, MPI_INT, 0, MPI_COMM_WORLD);

        int max = rbuf[0], min = rbuf[0];
        for (int i = 1; i < scounts[rank]; i++)
        {
            if (rbuf[i] > max)
                max = rbuf[i];
            if (rbuf[i] < min)
                min = rbuf[i];
        }

        if (!rank)
        {
            int ms;
            MPI_Recv(&ms, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (max < ms)
                max = ms;
            MPI_Recv(&ms, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            if (min > ms)
                min = ms;
            printf("max: %d, min: %d", max, min);
        }
        else
        {
            MPI_Send(&max, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&min, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }

        MPI_Finalize();
        return MPI_SUCCESS;
    }
}

namespace Task7
{
    int main(int argc, char** argv)
    {
        const int n = 10, m = 10;
        int size, rank;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Status status;

        int* matrix = (int*)malloc(n * m);
        int* vector = (int*)malloc(n);
        int *m_displs = (int*)malloc(size), *v_displs = (int*)malloc(size);
        int *m_scounts = (int*)malloc(size), *v_scounts = (int*)malloc(size);

        if (!rank)
        {
            printf("Matix:\n");
            for (int i = 0; i < n * m; i++)
            {
                matrix[i] = rand();
                printf("%d ", matrix[i]);
                if ((i + 1) % m == 0)
                    printf("\n");
            }
            printf("\n");

            printf("Vector:\n");
            for (int i = 0; i < n; i++)
            {
                vector[i] = rand();
                printf("%d ", vector[i]);
            }
            printf("\n");
        }

        for (int i = 0; i < size; i++)
        {
            if (i + 1 == size)
            {
                m_scounts[i] = m * (n - (i * (n / size)));
                v_scounts[i] = n - (i * (n / size));
            }
            else
            {
                m_scounts[i] = n * m / size;
                v_scounts[i] = n / size;
            }

            m_displs[i] = m * n * i / size;
            v_displs[i] = i * (n / size);
        }

        // int rMbuf[n * m], rVbuf[n];
        // MPI_Scatterv(matrix, m_scounts, m_displs, MPI_INT, &rMbuf, n * m, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Scatterv(vector, v_scounts, v_displs, MPI_INT, &rVbuf, n, MPI_INT, 0, MPI_COMM_WORLD);
        //
        // int result = 0;
        // for (int i = m_displs[rank]; i < m_displs[rank] + m_scounts[i]; i++)
        //     result += rMbuf[i];
        //
        // if (rank)
        //     MPI_Send(&result, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
        // else
        // {
        //     for (int i = 0; i < size - 1; ++i)
        //     {
        //         int ms[arr_size];
        //         MPI_Status st;
        //         MPI_Recv(&ms, rsize,MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &st);
        //         for (int j = 0; j < rsize; ++j)
        //             crr[j] += ms[j];
        //     }
        //     for (int i = 0; i < rsize; ++i)
        //         printf("%d ", crr[i]);
        // }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
} // TODO fix

namespace Task8
{
    int main(int argc, char** argv)
    {
        int n = 10;
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int* array = (int*)malloc(n);
        if (!rank)
        {
            printf("Start array:\n");
            for (int i = 0; i < n; ++i)
            {
                array[i] = rand();
                printf("%d ", array[i]);
            }
            printf("\n");
        }

        int* cropped_arr;
        int crop_size = rank + 1 == size ? n - (rank * (n / size)) : n / size;
        cropped_arr = (int*)malloc(crop_size);

        if (!rank)
        {
            for (int i = 1; i < size; ++i)
            {
                int cropped_size = i + 1 == size ? n - (i * (n / size)) : n / size;
                const int start = n / size * i;

                for (int j = 0; j < cropped_size; ++j)
                    cropped_arr[j] = array[start + j];

                MPI_Send(cropped_arr, cropped_size, MPI_INT, i, 0,MPI_COMM_WORLD);
            }

            for (int i = 0; i < crop_size; ++i)
                cropped_arr[i] = array[i];

            int* result = (int*)malloc(n);
            int k = 0;
            for (k = 0; k < crop_size; k++)
                result[k] = cropped_arr[k];
            for (int i = 1; i < size; ++i)
            {
                int cropped_size = i + 1 == size ? n - (i * (n / size)) : n / size;
                MPI_Status st;
                MPI_Recv(cropped_arr, cropped_size, MPI_INT, i, MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                for (int j = 0; j < cropped_size; j++, k++)
                    result[k] = cropped_arr[j];
            }

            printf("Result:\n");
            for (int i = 0; i < n; ++i)
                printf("%d ", result[i]);
        }
        else
        {
            MPI_Status st;
            MPI_Recv(cropped_arr, crop_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &st);
            printf("Thread %d: ", rank);
            for (int i = 0; i < crop_size; ++i)
                printf("%d ", cropped_arr[i]);

            MPI_Send(cropped_arr, crop_size, MPI_INT, 0, 0,MPI_COMM_WORLD);
        }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

namespace Task9
{
    int main(int argc, char** argv)
    {
        int arr[] = {15, 234, 765, 4, 564, 823, 45, 24, 312, 21, 34, 165};
        constexpr int arr_size = 12, cropped_size = 4;
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int cropped_arr[cropped_size];
        if (!rank)
        {
            for (int i = 1; i < size; ++i)
            {
                for (int j = 0; j < cropped_size; ++j)
                    cropped_arr[j] = arr[j + i * cropped_size];

                MPI_Send(cropped_arr, cropped_size, MPI_INT, i, 0,MPI_COMM_WORLD);
            }

            MPI_Status st;
            int result[arr_size];
            for (int j = 0; j < cropped_size; ++j)
                result[j] = arr[arr_size - j - 1];

            for (int i = size - 1; i > 0; i--)
            {
                MPI_Recv(cropped_arr, cropped_size, MPI_INT, i, MPI_ANY_TAG,MPI_COMM_WORLD, &st);
                for (int j = 0; j < cropped_size; ++j)
                    result[i * cropped_size + j] = cropped_arr[j];
            }

            for (int i = 0; i < arr_size; ++i)
                printf("%d ", result[i]);
        }
        else
        {
            MPI_Status st;
            MPI_Recv(cropped_arr, cropped_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &st);
            printf("Thread %d: ", rank);
            for (int i = 0; i < cropped_size; ++i)
                printf("%d ", cropped_arr[i]);

            int reverse[cropped_size];
            for (int i = 0; i < cropped_size; ++i)
                reverse[i] = cropped_arr[cropped_size - 1 - i];

            printf("Reverse: ");
            for (int i = 0; i < cropped_size; ++i)
                printf("%d ", reverse[i]);

            MPI_Send(reverse, cropped_size, MPI_INT, 0, 0,MPI_COMM_WORLD);
        }

        MPI_Finalize();

        return MPI_SUCCESS;
    }
}

int main(int argc, char* argv[])
{
    Task8::main(argc, argv);
}
