#include <iostream>
#include <pthread.h>
#include <semaphore.h>
#include<windows.h>
#include <stdlib.h>
#include<iomanip>
#include <omp.h>
#include <fstream>
#include <nmmintrin.h> //SSSE4.2
using namespace std;


//------------------------------------------ 线程控制变量 ------------------------------------------
typedef struct
{
    int t_id; //线程编号
} threadParam_t;  //传参打包

//barrier 定义
pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;


// ------------------------------------------ 全局计算变量 ------------------------------------------
int N ; 
int Step[6]={500,1000,1500,2000,2500,3000};
float** ddata;
float** matrix;
const int L = 100; 
const int LOOP = 1; 
const int THREAD_NUM = 16; //线程数
ofstream res_stream;
void init_data(int x);
void init_matrix();
void calculate_serial();
void calculate_serial_column();
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

// ------------------------------------------ 初始化data ------------------------------------------
// 保证每次数据都是一致的
void init_data(int x)
{
    N=x;
    ddata=new float*[N];
    for(int i=0;i<N;i++)
        {
            ddata[i]=new float[N];
        }
    for (int i = 0; i < N; i++)
        for (int j = i; j < N; j++)
            ddata[i][j]=0;    
    for (int i = 0; i < N; i++)
        for (int j = i; j < N; j++)
            ddata[i][j] = rand() * 1.0 / RAND_MAX * L;
    for (int i = 0; i < N - 1; i++)
        for (int j = i + 1; j < N; j++)
            for (int k = 0; k < N; k++)
                ddata[j][k] += ddata[i][k];
}

// 用data初始化matrix
void init_matrix()
{
    matrix=new float*[N];
    for(int i=0;i<N;i++)
        {
            matrix[i]=new float[N];
        }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix[i][j] = ddata[i][j];
}

int main()
{
    long long head, tail,freq;
    float seconds;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    float time = 0;
    //res_stream.open("result.txt", ios::out);
    for(int ii=0;ii<6;ii++)
    {
        init_data(Step[ii]);
        cout<<"N:"<<N<<endl;
        //res_stream<<"N:"<<N<<endl;
    // ====================================== serial_row ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_serial();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "serial_row:" << time / LOOP << "ms" << endl;
        //res_stream << "serial_row:" << time / LOOP << "ms" << endl;
	    //print_matrix();     
        // ====================================== SIMD ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_simd();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "simd:" << time / LOOP << "ms" << endl;
        //res_stream << "simd:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== thread ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_pthread();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "thread:" << time / LOOP << "ms" << endl;
        //res_stream << "thread:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== openmp_row ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            calculate_openmp();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        cout << "openmp_row:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_row:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== openmp_SIMD(手写) ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_openmp_SIMD();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "openmp_SIMD_arti:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_SIMD_arti:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== openmp_autoSIMD ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_openmp_autoSIMD();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "openmp_SIMD_auto:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_SIMD_auto:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
  
        // ====================================== openmp_static ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_openmp_static();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "openmp_static:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_static:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== openmp_dynamic ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_openmp_dynamic();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "openmp_dynamic:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_dynamic:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
        // ====================================== openmp_guided ======================================
        time = 0;
        for (int i = 0; i < LOOP; i++)
        {
            init_matrix();
		    //print_matrix();
            QueryPerformanceCounter((LARGE_INTEGER *)&head);
            //calculate_openmp_guided();
            QueryPerformanceCounter((LARGE_INTEGER *)&tail);
            time+=(tail - head) * 1000.0 / freq;
        }
        //cout << "openmp_guided:" << time / LOOP << "ms" << endl;
        //res_stream << "openmp_guided:" << time / LOOP << "ms" << endl;
	    //print_matrix(); 
    }
   // res_stream.close();
   system("pause");
    return 0;
}

// 串行算法
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


// neon并行算法
void calculate_simd()
{
    for (int k = 0; k < N; k++)
    {
            for (int j = k + 1; j < N; j++)
        {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;
           for(int i=k+1;i<N;i++)
           {
                __m128 vaik= _mm_set1_ps(matrix[i][k]);
                int j=k+1;
                for(j;j+4<=N;j=j+4)
                {
                    __m128 vakj=_mm_loadu_ps(matrix[k]+j);
                    __m128 vaij=_mm_loadu_ps(matrix[i]+j);
                    __m128 vx=_mm_mul_ps(vakj,vaik);
                    vaij=_mm_sub_ps(vaij,vx);
                    _mm_storeu_ps(matrix[i]+j,vaij);
                }
                for(j;j<N;j++)
                {
                    matrix[i][j]=matrix[i][j]-matrix[k][j]*matrix[i][k];
                }
                matrix[i][k]=0;
            }
    }
}

//线程函数定义
void *threadFunc3_barrier(void *param) 
{
	threadParam_t *p = (threadParam_t*)param;
	int t_id = p->t_id;
	
	for (int k = 0; k <N; ++k)
	{
		// t_id 为 0 的线程做除法操作，其它工作线程先等待
		// 这里只采用了一个工作线程负责除法操作，同学们可以尝试采用多个工作线程完成除法操作
		if (t_id == 0)
		{
			for (int j=k+1; j<N; j++)
			{
				matrix[k][j ] = matrix[k][j] / matrix[k][k];
			}
			matrix[k][k] = 1.0;
		}
		
		//第一个同步点
		pthread_barrier_wait(&barrier_Division);
		
		//循环划分任务（同学们可以尝试多种任务划分方式）
		for(int i=k+1+t_id; i <N; i += THREAD_NUM)
		{
			//消去
			for (int j = k + 1; j <N; ++j)
			{
				matrix[i ][ j ] = matrix[i][j] - matrix[i][k] * matrix[k][j];
			}
			matrix[i][k]=0.0;
		}
		
		// 第二个同步点
		pthread_barrier_wait(&barrier_Elimination);
	}
	pthread_exit(NULL);
    return NULL;
}


void calculate_pthread() {
	//初始化 barrier
	pthread_barrier_init(&barrier_Division, NULL, THREAD_NUM);
	pthread_barrier_init(&barrier_Elimination, NULL, THREAD_NUM);
	//创建线程
	pthread_t handles[THREAD_NUM];// 创建对应的 Handle
	threadParam_t param[THREAD_NUM];// 创建对应的线程数据结构
	for(int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], NULL, threadFunc3_barrier, (void *)(&param[t_id]));
	}
	for(int t_id = 0; t_id < THREAD_NUM; t_id++)
	{
		pthread_join(handles[t_id],NULL);
	}
	
	//销毁所有的 barrier
	pthread_barrier_destroy(&barrier_Division);
	pthread_barrier_destroy(&barrier_Elimination);
}


// ====================================== openMP======================================


// 指导手册示例
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
//手写openMP+SIMD
void calculate_openmp_SIMD()
{
    int i, j, k;
    float tmp;
#pragma omp parallel num_threads(THREAD_NUM) default(none) private(i, j, k, tmp) shared(matrix, N)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {
            __m128 vt=_mm_set1_ps(matrix[k][k]);
            int j=k+1;
            for(j=k+1;j+4<=N;j=j+4)
            {
                __m128 va = _mm_loadu_ps(matrix[k]+j);
                va=_mm_div_ps(va,vt);
                _mm_storeu_ps(matrix[k]+j,va);
            }
           for(;j<N;j++)
           {
                matrix[k][j]==matrix[k][j]/matrix[k][k];
           }
           matrix[k][k]=1.0;

        }
#pragma omp for
           for(int i=k+1;i<N;i++)
           {
                __m128 vaik= _mm_set1_ps(matrix[i][k]);
                int j=k+1;
                for(j;j+4<=N;j=j+4)
                {
                    __m128 vakj=_mm_loadu_ps(matrix[k]+j);
                    __m128 vaij=_mm_loadu_ps(matrix[i]+j);
                    __m128 vx=_mm_mul_ps(vakj,vaik);
                    vaij=_mm_sub_ps(vaij,vx);
                    _mm_storeu_ps(matrix[i]+j,vaij);
                }
                for(j;j<N;j++)
                {
                    matrix[i][j]=matrix[i][j]-matrix[k][j]*matrix[i][k];
                }
                matrix[i][k]=0;
            }
    }
}
//openMP自动SIMD
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
//openmp static分配
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
//openmp dynamic分配
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
//openmp guided分配
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
//openmp水平划分
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