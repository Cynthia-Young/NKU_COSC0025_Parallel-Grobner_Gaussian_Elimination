#include <iostream>
#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include<iomanip>
#include <omp.h>
#include <fstream>
#include <arm_neon.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
using namespace std;


//------------------------------------------ �߳̿��Ʊ��� ------------------------------------------
typedef struct
{
    int t_id; 
} threadParam_t;  

const int THREAD_NUM = 4; 

//barrier
pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;


// --------------------------------------:---- ------------------------------------------
const int N=500;
const int L = 100;
const int LOOP = 1;

float ddata[N][N];
float matrix[N][N];

void init_data(int x);
void init_matrix();
void calculate_serial();
void calculate_simd();
void calculate_openmp();
void calculate_pthread();
void calculate_openmp_SIMD();
void calculate_openmp_autoSIMD();
void calculate_openmp_static();
void calculate_openmp_dynamic();
void calculate_openmp_guided();
void openmp_row();
void openmp_col();
void print_matrix();

// ------------------------------------------ ��ʼ��data ------------------------------------------

void init_data()
{
	for (int i = 0; i < N; i++)
        for (int j = i; j < N; j++)
            ddata[i][j] = rand() * 1.0 / RAND_MAX * L;
    for (int i = 0; i < N - 1; i++)
        for (int j = i + 1; j < N; j++)
            for (int k = 0; k < N; k++)
                ddata[j][k] += ddata[i][k];
}

void init_matrix()
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix[i][j] = ddata[i][j];
}

int main()
{
    struct timeval start;
    struct timeval end;
    float time = 0;
    cout<<THREAD_NUM<<endl;
    init_data();
    cout<<"N:"<<N<<endl;
    // ====================================== serial_row ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_serial();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "serial_row:" << time / LOOP << "ms" << endl;
	    //print_matrix();

        // ====================================== SIMD ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_simd();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "simd:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== pthread ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_pthread();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "thread:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== openmp_row ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_openmp();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "openmp_row:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== openmp_SIMD(��д) ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_openmp_SIMD();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "openmp_SIMD_arti:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== openmp_autoSIMD ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_openmp_autoSIMD();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "openmp_SIMD_auto:" << time / LOOP << "ms" << endl;
	    //print_matrix();

        // ====================================== openmp_static ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            calculate_openmp_static();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        cout << "openmp_static:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== openmp_dynamic ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            calculate_openmp_dynamic();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
       cout << "openmp_dynamic:" << time / LOOP << "ms" << endl;
	    //print_matrix();
        // ====================================== openmp_guided ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            gettimeofday(&start, NULL);
            //calculate_openmp_guided();
            gettimeofday(&end, NULL);
			time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
        }
        //cout << "openmp_guided:" << time / LOOP << "ms" << endl;
	    //print_matrix();
   //system("pause");
    return 0;
}

void calculate_serial()
{ 
    for (int k = 0; k < N; k++)
    {
        for (int j = k + 1; j < N; j++)
        {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;
        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}


// neon
void calculate_simd()
{

    for (int k = 0; k < N; k++)
    {
        float32x4_t Akk = vmovq_n_f32(matrix[k][k]);
        int j;
        for (j = k + 1; j + 3 < N; j += 4)
        {
            float32x4_t Akj = vld1q_f32(matrix[k] + j);
            Akj = vdivq_f32(Akj, Akk);
            vst1q_f32(matrix[k] + j, Akj);
        }
        for (; j < N; j++)
        {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;
        //除法
        for (int i = k + 1; i < N; i++)
        {
            float32x4_t Aik = vmovq_n_f32(matrix[i][k]);
            for (j = k + 1; j + 3 < N; j += 4)
            {
                float32x4_t Akj = vld1q_f32(matrix[k] + j);
                float32x4_t Aij = vld1q_f32(matrix[i] + j);
                float32x4_t AikMulAkj = vmulq_f32(Aik, Akj);
                Aij = vsubq_f32(Aij, AikMulAkj);
                vst1q_f32(matrix[i] + j, Aij);
            }
            for (; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}


void *threadFunc3_barrier(void *param)
{
	threadParam_t *p = (threadParam_t*)param;
	int t_id = p->t_id;

	for (int k = 0; k <N; ++k)
	{
		
		if (t_id == 0)
		{
			for (int j=k+1; j<N; j++)
			{
				matrix[k][j ] = matrix[k][j] / matrix[k][k];
			}
			matrix[k][k] = 1.0;
		}


		pthread_barrier_wait(&barrier_Division);
		for(int i=k+1+t_id; i <N; i += THREAD_NUM)
		{
			
			for (int j = k + 1; j <N; ++j)
			{
				matrix[i ][ j ] = matrix[i][j] - matrix[i][k] * matrix[k][j];
			}
			matrix[i][k]=0.0;
		}
		pthread_barrier_wait(&barrier_Elimination);
	}
	pthread_exit(NULL);
    return NULL;
}

void calculate_pthread() {
	//barrier
	pthread_barrier_init(&barrier_Division, NULL, THREAD_NUM);
	pthread_barrier_init(&barrier_Elimination, NULL, THREAD_NUM);

	pthread_t handles[THREAD_NUM];
	threadParam_t param[THREAD_NUM];
	for(int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], NULL, threadFunc3_barrier, (void *)(&param[t_id]));
	}
	for(int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		pthread_join(handles[t_id],NULL);
	}


	pthread_barrier_destroy(&barrier_Division);
	pthread_barrier_destroy(&barrier_Elimination);
}


// ====================================== openMP======================================

void calculate_openmp()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for
        for (i = k + 1; i < N; i++)
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}
//openMP+SIMD

void calculate_openmp_SIMD()
{
    int i, j, k;
    float tmp;
#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k, tmp) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {
            float32x4_t Akk = vmovq_n_f32(matrix[k][k]);
        int j;
        for (j = k + 1; j + 3 < N; j += 4)
        {
            float32x4_t Akj = vld1q_f32(matrix[k] + j);
            Akj = vdivq_f32(Akj, Akk);
            vst1q_f32(matrix[k] + j, Akj);
        }
        for (; j < N; j++)
        {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;

        }
#pragma omp for
           for(int i=k+1;i<N;i++)
           {
                float32x4_t Aik = vmovq_n_f32(matrix[i][k]);
            for (j = k + 1; j + 3 < N; j += 4)
            {
                float32x4_t Akj = vld1q_f32(matrix[k] + j);
                float32x4_t Aij = vld1q_f32(matrix[i] + j);
                float32x4_t AikMulAkj = vmulq_f32(Aik, Akj);
                Aij = vsubq_f32(Aij, AikMulAkj);
                vst1q_f32(matrix[i] + j, Aij);
            }
            for (; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
            }
    }
}
//openMPSIMD

void calculate_openmp_autoSIMD()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {
#pragma omp simd
            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for
        for (i = k + 1; i < N; i++)
        {
#pragma omp simd
            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}
//openmp static
void calculate_openmp_static()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for schedule(static,N/(3*THREAD_NUM))
        for (i = k + 1; i < N; i++)
        {
            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }

            matrix[i][k] = 0;
        }
    }
}
//openmp dynamic
void calculate_openmp_dynamic()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for schedule(dynamic,100)
        for (i = k + 1; i < N; i++)
        {
            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}
//openmp guided
void calculate_openmp_guided()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for schedule(guided)
        for (i = k + 1; i < N; i++)
        {
            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}

void openmp_row()
{
    int i, j, k;

#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1.0;
        }
#pragma omp for
        for (i = k + 1; i < N; i++)
        {

            for (j = k + 1; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}


void print_matrix()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}
