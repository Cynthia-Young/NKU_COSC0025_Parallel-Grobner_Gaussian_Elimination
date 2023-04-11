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
    int t_id; //�̱߳��
    int i;
} threadParam_t;  //���δ��

const int THREAD_NUM = 16; //�߳���

//barrier ����
pthread_barrier_t barrier_a;
pthread_barrier_t barrier_b;
pthread_barrier_t barrier_c;



const int Columnnum = 3799;    //�������� 130 254 562 1011 2362 3799 8399 23045
//const int Columnnum = 8399;
const int Rnum = 2759;    //������Ԫ�� 22 106 170 539 1226 2759 6375 18748
//const int Rnum = 6375;
const int Enum = 1953;    //����Ԫ�� 8 53 53 263 453 1953 4535 14325
//const int Enum =4535;
const int ArrayColumn = 119;           //5 8 18 32 74 119 263 721
//8 = 254 / 32 ceiling  ÿһ��λ����1bit��ʾ��int����32 bit��������Ҫ8������
//const int ArrayColumn = 263;
const int leftbit = 9;              //30 2 14 13 6 9 17 26
//32 - (1 + 253 % 32) = 2
//const int leftbit = 17;                //32 - (1 + 43576 % 32) = 7
//unsigned
int R[Columnnum][ArrayColumn];  //��Ԫ�Ӿ��� ��಻����columnnum��
//unsigned 
int E[Enum][ArrayColumn];       //����Ԫ�о��� һ��enum��

int First[Enum];

bitset<32> MyBit(0);

int Find_First(int index) {
    int j = 0;
    int cnt = 0;
    while (E[index][j] == 0) {
        j++;
        if (j == ArrayColumn) break;
    }
    if (j == ArrayColumn) return -1; //���ȫΪ0 ����-1
    unsigned int tmp = E[index][j]; //
    while (tmp != 0) {
        tmp = tmp >> 1;
        cnt++;
    }//��ʱ�õ�j����Ԫ�ڵڼ���int�cnt���������int����ĵڼ���bit
    return Columnnum - 1 - ((j + 1) * 32 - cnt - leftbit);
}

void Init_R() {
    unsigned int a;
    ifstream infile("C:\\Users\\vivia\\Desktop\\��Ԫ��6.txt");

    char fin[5000] = { 0 };
    int index;
    while (infile.getline(fin, sizeof(fin))) { //����ȡ��
        std::stringstream line(fin);
        bool flag = 0;
        while (line >> a) { //һ����ȡ��
            if (flag == 0) { //��һ����
                index = a; //ȡ���ĵ�һ����������Ԫ�ӵ����� ����ԪΪ819����ô������Ԫ���о����ڵ�819.
                flag = 1;
            }
            int k = a % 32; //������ǵڼ���
            int j = a / 32; //���ڵڼ���
            int temp = 1 << k; //һ��32λ
            R[index][ArrayColumn - 1 - j] += temp; //�Ѷ�Ӧ����1bit��Ϊ1.
        }
    }
}

void Init_E() {
    unsigned int a;
    ifstream infile("C:\\Users\\vivia\\Desktop\\����Ԫ��6.txt");

    char fin[5000] = { 0 };
    int index = 0;
    while (infile.getline(fin, sizeof(fin))) {
        std::stringstream line(fin);
        int flag = 0;
        while (line >> a) {
            if (flag == 0) {
                First[index] = a; //����Ԫ�а�����˳��������First�����¼ÿ������Ԫ�е���Ԫ
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
//-----------------------------------------�����㷨----------------------------------------------------
void Serial() {
    for (int i = 0; i < Enum; i++) { //�Ա���Ԫ����ѭ��
        while (First[i] != -1) {  //�������ȫ0�ı���Ԫ�� 
            if (!Is_NULL(First[i])) { //����˱���Ԫ�е���Ԫ �ж�Ӧ����Ԫ��
                XOR(i, First[i]); //�����������
                First[i] = Find_First(i); //���������������Ԫ�е���Ԫ
            }
            else {
                Set_R(i, First[i]); //���������Ԫ����Ϊ��Ԫ��
                break;
            }
        }
    }
}


//---------------------------------------------pthread---------------------------------------------------------
//�̺߳�������
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
        First[i] = Find_First(i); //���������������Ԫ�е���Ԫ
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
    //��ʼ�� barrier
    pthread_barrier_init(&barrier_a, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_b, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_c, NULL, THREAD_NUM);


    //�����߳�
    pthread_t handles[THREAD_NUM];// ������Ӧ�� Handle
    threadParam_t param[THREAD_NUM];// ������Ӧ���߳����ݽṹ

    for (int i = 0; i < Enum; i++) { //�Ա���Ԫ����ѭ��
        while (First[i] != -1) {  //�������ȫ0�ı���Ԫ��
            if (!Is_NULL(First[i])) { //����˱���Ԫ�е���Ԫ �ж�Ӧ����Ԫ��
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


    

    //�������е� barrier
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
        First[i] = Find_First(i); //���������������Ԫ�е���Ԫ
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
    //��ʼ�� barrier
    pthread_barrier_init(&barrier_a, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_b, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_c, NULL, THREAD_NUM);


    //�����߳�
    pthread_t handles[THREAD_NUM];// ������Ӧ�� Handle
    threadParam_t param[THREAD_NUM];// ������Ӧ���߳����ݽṹ

    for (int i = 0; i < Enum; i++) { //�Ա���Ԫ����ѭ��
        while (First[i] != -1) {  //�������ȫ0�ı���Ԫ��
            if (!Is_NULL(First[i])) { //����˱���Ԫ�е���Ԫ �ж�Ӧ����Ԫ��
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




    //�������е� barrier
    pthread_barrier_destroy(&barrier_a);
    pthread_barrier_destroy(&barrier_b);
    pthread_barrier_destroy(&barrier_c);
}


void Parallel_simd()
{
    int i, j, k = 0;
    for (int i = 0; i < Enum; i++) { //�Ա���Ԫ����ѭ��
        while (First[i] != -1) {  //�������ȫ0�ı���Ԫ�� 
            if (!Is_NULL(First[i])) { //����˱���Ԫ�е���Ԫ �ж�Ӧ����Ԫ��
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
    { //�Ա���Ԫ����ѭ��
        while (First[i] != -1)
        {//�������ȫ0�ı���Ԫ��
            if (!Is_NULL(First[i]))
            {   //����˱���Ԫ�е���Ԫ �ж�Ӧ����Ԫ��
                //XOR(i,First[i]); //�����������
#pragma omp parallel for num_threads(THREAD_NUM) private(j) shared(E,R)
                for (j = 0; j < ArrayColumn; j++)
                {
                    E[i][j] = E[i][j] ^ R[First[i]][j];
                }
                First[i] = Find_First(i); //���������������Ԫ�е���Ԫ
            }
            else {
                //Set_R(i,First[i]); //���������Ԫ����Ϊ��Ԫ��
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
