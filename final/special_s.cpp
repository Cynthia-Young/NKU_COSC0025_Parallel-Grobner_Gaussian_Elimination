#include <iostream>
#include <fstream>
#include <sstream>
#include<bitset>
#include<windows.h>
//#include<sys/time.h>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>
#include <omp.h>
#include <immintrin.h>

using namespace std;

#pragma comment(lib, "D:\\Visual Studio 2022\\VC\\lib\\x64\\pthreadVC2.lib")
typedef struct
{
    int t_id; //线程编号
    int i;
} threadParam_t;  //传参打包

const int THREAD_NUM = 16; //线程数

//barrier 定义
pthread_barrier_t barrier_a;
pthread_barrier_t barrier_b;
pthread_barrier_t barrier_c;



const int Columnnum = 3799;    //矩阵列数 130 254 562 1011 2362 3799 8399 23045
//const int Columnnum = 8399;
const int Rnum = 2759;    //非零消元子 22 106 170 539 1226 2759 6375 18748
//const int Rnum = 6375;
const int Enum = 1953;    //被消元行 8 53 53 263 453 1953 4535 14325
//const int Enum =4535;
const int ArrayColumn = 119;           //5 8 18 32 74 119 263 721
//8 = 254 / 32 ceiling  每一个位置用1bit表示，int型有32 bit，所以需要8个整形
//const int ArrayColumn = 263;
const int leftbit = 9;              //30 2 14 13 6 9 17 26
//32 - (1 + 253 % 32) = 2
//const int leftbit = 17;                //32 - (1 + 43576 % 32) = 7
//unsigned
int R[Columnnum][ArrayColumn];  //消元子矩阵 最多不超过columnnum行
//unsigned 
int E[Enum][ArrayColumn];       //被消元行矩阵 一共enum行

int First[Enum];

bitset<32> MyBit(0);

int Find_First(int index) {
    int j = 0;
    int cnt = 0;
    while (E[index][j] == 0) {
        j++;
        if (j == ArrayColumn) break;
    }
    if (j == ArrayColumn) return -1; //如果全为0 返回-1
    unsigned int tmp = E[index][j]; //
    while (tmp != 0) {
        tmp = tmp >> 1;
        cnt++;
    }//此时得到j是首元在第几个int里，cnt代表是这个int里面的第几个bit
    return Columnnum - 1 - ((j + 1) * 32 - cnt - leftbit);
}

void Init_R() {
    unsigned int a;
    ifstream infile("C:\\Users\\vivia\\Desktop\\消元子6.txt");

    char fin[5000] = { 0 };
    int index;
    while (infile.getline(fin, sizeof(fin))) { //按行取出
        std::stringstream line(fin);
        bool flag = 0;
        while (line >> a) { //一个个取数
            if (flag == 0) { //第一个数
                index = a; //取到的第一个数即是消元子的索引 如首元为819，那么它在消元子中就排在第819.
                flag = 1;
            }
            int k = a % 32; //分组后是第几个
            int j = a / 32; //它在第几组
            int temp = 1 << k; //一共32位
            R[index][ArrayColumn - 1 - j] += temp; //把对应的那1bit设为1.
        }
    }
}

void Init_E() {
    unsigned int a;
    ifstream infile("C:\\Users\\vivia\\Desktop\\被消元行6.txt");

    char fin[5000] = { 0 };
    int index = 0;
    while (infile.getline(fin, sizeof(fin))) {
        std::stringstream line(fin);
        int flag = 0;
        while (line >> a) {
            if (flag == 0) {
                First[index] = a; //被消元行按读入顺序排序，用First数组记录每个被消元行的首元
                flag = 1;
            }
            int k = a % 32;
            int j = a / 32;
            int temp = 1 << k;
            E[index][ArrayColumn - 1 - j] += temp;
        }
        index++;
    }
}

bool Is_NULL(int index) {
    for (int j = 0; j < ArrayColumn; j++) {
        if (R[index][j] != 0) return 0;
    }
    return 1;
}

void Set_R(int eindex, int rindex) {
    for (int j = 0; j < ArrayColumn; j++) {
        R[rindex][j] = E[eindex][j];
    }
}

void XOR(int eindex, int rindex) {
    for (int j = 0; j < ArrayColumn; j++) {
        E[eindex][j] = E[eindex][j] ^ R[rindex][j];
    }
}
//-----------------------------------------串行算法----------------------------------------------------
void Serial() {
    for (int i = 0; i < Enum; i++) { //对被消元行做循环
        while (First[i] != -1) {  //如果不是全0的被消元行 
            if (!Is_NULL(First[i])) { //如果此被消元行的首元 有对应的消元子
                XOR(i, First[i]); //做减法（异或）
                First[i] = Find_First(i); //重新设置这个被消元行的首元
            }
            else {
                Set_R(i, First[i]); //将这个被消元行升为消元子
                break;
            }
        }
    }
}


//---------------------------------------------pthread---------------------------------------------------------
//线程函数定义
void *thread_xor(void *param)
{
    threadParam_t *p = (threadParam_t*)param;
    int t_id = p->t_id;
    int i = p->i;
    for (int j = t_id; j < ArrayColumn; j += THREAD_NUM) {

        E[i][j] = E[i][j] ^ R[First[i]][j];
    }
    pthread_barrier_wait(&barrier_a);
    if (t_id == 0)
    {
        First[i] = Find_First(i); //重新设置这个被消元行的首元
    }
    pthread_barrier_wait(&barrier_b);

    pthread_exit(NULL);
    return NULL;
}

void* thread_set(void* param)
{
    threadParam_t* p = (threadParam_t*)param;
    int t_id = p->t_id;
    int i = p->i;

    for (int j = t_id; j < ArrayColumn; j += THREAD_NUM) {
        R[First[i]][j] = E[i][j];
    }
    pthread_barrier_wait(&barrier_c);

    pthread_exit(NULL);
    return NULL;
}

void Parallel_pthread(){
    //初始化 barrier
    pthread_barrier_init(&barrier_a, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_b, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_c, NULL, THREAD_NUM);


    //创建线程
    pthread_t handles[THREAD_NUM];// 创建对应的 Handle
    threadParam_t param[THREAD_NUM];// 创建对应的线程数据结构

    for (int i = 0; i < Enum; i++) { //对被消元行做循环
        while (First[i] != -1) {  //如果不是全0的被消元行
            if (!Is_NULL(First[i])) { //如果此被消元行的首元 有对应的消元子
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    param[t_id].t_id = t_id;
                    param[t_id].i = i;
                    pthread_create(&handles[t_id], NULL, thread_xor, (void*)(&param[t_id]));
                }
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_join(handles[t_id], NULL);
                }
            }
            else {
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_create(&handles[t_id], NULL, thread_set, (void*)(&param[t_id]));
                }
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_join(handles[t_id], NULL);
                }
                break;
            }
        }
    }


    

    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_a);
    pthread_barrier_destroy(&barrier_b);
    pthread_barrier_destroy(&barrier_c);
}

//---------------------------------------------pthread+simd---------------------------------------------------------
void* thread_xor_simd(void* param)
{
    threadParam_t* p = (threadParam_t*)param;
    int t_id = p->t_id;
    int i = p->i;
    for (int j = t_id; j < ArrayColumn; j += THREAD_NUM) {

        E[i][j] = E[i][j] ^ R[First[i]][j];
    }
    pthread_barrier_wait(&barrier_a);
    if (t_id == 0)
    {
        First[i] = Find_First(i); //重新设置这个被消元行的首元
    }
    pthread_barrier_wait(&barrier_b);

    pthread_exit(NULL);
    return NULL;
}

void* thread_set_simd(void* param)
{
    threadParam_t* p = (threadParam_t*)param;
    int t_id = p->t_id;
    int i = p->i;

    for (int j = t_id; j < ArrayColumn; j += THREAD_NUM) {
        R[First[i]][j] = E[i][j];
    }
    pthread_barrier_wait(&barrier_c);

    pthread_exit(NULL);
    return NULL;
}

void Parallel_pthread_simd() {
    //初始化 barrier
    pthread_barrier_init(&barrier_a, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_b, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_c, NULL, THREAD_NUM);


    //创建线程
    pthread_t handles[THREAD_NUM];// 创建对应的 Handle
    threadParam_t param[THREAD_NUM];// 创建对应的线程数据结构

    for (int i = 0; i < Enum; i++) { //对被消元行做循环
        while (First[i] != -1) {  //如果不是全0的被消元行
            if (!Is_NULL(First[i])) { //如果此被消元行的首元 有对应的消元子
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    param[t_id].t_id = t_id;
                    param[t_id].i = i;
                    pthread_create(&handles[t_id], NULL, thread_xor_simd, (void*)(&param[t_id]));
                }
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_join(handles[t_id], NULL);
                }
            }
            else {
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_create(&handles[t_id], NULL, thread_set_simd, (void*)(&param[t_id]));
                }
                for (int t_id = 0; t_id < THREAD_NUM; t_id++)
                {
                    pthread_join(handles[t_id], NULL);
                }
                break;
            }
        }
    }




    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_a);
    pthread_barrier_destroy(&barrier_b);
    pthread_barrier_destroy(&barrier_c);
}


void Parallel_simd()
{
    int i, j, k = 0;
    for (int i = 0; i < Enum; i++) { //对被消元行做循环
        while (First[i] != -1) {  //如果不是全0的被消元行 
            if (!Is_NULL(First[i])) { //如果此被消元行的首元 有对应的消元子
                for (j = 0; j + 3 < ArrayColumn; j += 4)
                {
                    __m128 Rfj = _mm_loadu_ps((float*)R[First[i]] + j);
                    __m128 Eij = _mm_loadu_ps((float*)E[i] + j);
                    Eij = _mm_xor_ps(Rfj, Eij);
                    _mm_storeu_ps((float*)E[i] + j, Eij);
                }
                for (j; j < ArrayColumn; j++)
                {
                    E[i][j] = E[i][j] ^ R[First[i]][j];
                }
                First[i] = Find_First(i);
            }
            else {
                for (j = 0; j + 3 < ArrayColumn; j += 4)
                {
                    __m128 Eij = _mm_loadu_ps((float*)E[i] + j);
                    _mm_storeu_ps((float*)R[First[i]] + j, Eij);
                }
                for (j; j < ArrayColumn; j++)
                {
                    R[First[i]][j] = E[i][j];
                }

                break;
            }
        }
    }
}

void Parallel_openmp()
{
    int i, j, k = 0;
    for (i = 0; i < Enum; i++)
    { //对被消元行做循环
        while (First[i] != -1)
        {//如果不是全0的被消元行
            if (!Is_NULL(First[i]))
            {   //如果此被消元行的首元 有对应的消元子
                //XOR(i,First[i]); //做减法（异或）
#pragma omp parallel for num_threads(THREAD_NUM) private(j) shared(E,R)
                for (j = 0; j < ArrayColumn; j++)
                {
                    E[i][j] = E[i][j] ^ R[First[i]][j];
                }
                First[i] = Find_First(i); //重新设置这个被消元行的首元
            }
            else {
                //Set_R(i,First[i]); //将这个被消元行升为消元子
#pragma omp parallel for num_threads(THREAD_NUM) private(j) shared(E,R)
                for (j = 0; j < ArrayColumn; j++)
                {
                    R[First[i]][j] = E[i][j];
                }
                break;
            }
        }
    }
}

void Parallel_openmp_simd()
{
    int i, j, k = 0;
    for (i = 0; i < Enum; i++)
    {
        while (First[i] != -1)
        {
            if (!Is_NULL(First[i]))
            {
#pragma omp parallel for num_threads(THREAD_NUM) private(j) shared(E,R)
                for (j = 0; j + 3 < ArrayColumn; j += 4)
                {
                    __m128 Rfj = _mm_loadu_ps((float *)R[First[i]] + j);
                    __m128 Eij = _mm_loadu_ps((float*)E[i] + j);
                    Eij = _mm_xor_ps(Rfj, Eij);
                    _mm_storeu_ps((float*)E[i] + j, Eij);
                }
                for (j; j < ArrayColumn; j++)
                {
                    E[i][j] = E[i][j] ^ R[First[i]][j];
                }
                First[i] = Find_First(i);
            }
            else {
#pragma omp parallel for num_threads(THREAD_NUM) private(j) shared(E,R)
                for (j = 0; j + 3 < ArrayColumn; j += 4)
                {
                    __m128 Eij = _mm_loadu_ps((float*)E[i] + j);
                    _mm_storeu_ps((float*)R[First[i]] + j,  Eij);
                }
                for (j; j < ArrayColumn; j++)
                {
                    R[First[i]][j] = E[i][j];
                }

                break;
            }
        }
    }
}


void Print() {//Print the answer
    for (int i = 0; i < Enum; i++) {
        if (First[i] == -1) {
            cout << endl;
            continue;
        }
        for (int j = 0; j < ArrayColumn; j++) {
            if (E[i][j] == 0) continue;
            MyBit = E[i][j];//MENTION: bitset manipulates from the reverse direction
            for (int k = 31; k >= 0; k--) {
                if (MyBit.test(k)) {
                    cout << 32 * (ArrayColumn - j - 1) + k << ' ';
                }
            }
        }
        cout << endl;
    }
}

int main()
{

    long long head, tail, freq;
    float seconds;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    float time = 0;


    Init_R();
    Init_E();
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for(int i=0;i<5;i++)
    //Serial();
    //Parallel_simd();
    //Parallel_openmp_simd();
    //Parallel_openmp();
    Parallel_pthread();

    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    time += (tail - head) * 1000.0 / freq;
    cout << "Enum: " << Enum << ", Time: " << time/5 << "ms" << endl;
    //Print();

    return 0;
}
