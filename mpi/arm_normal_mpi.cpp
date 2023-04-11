#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <arm_neon.h>



//#define _PRINT
#define _TEST

using namespace std;

// ============================================== ������� ==============================================
int N;
const int L = 100;
int LOOP = 1;
float **origin_data;
float **matrix = nullptr;

int NUM_THREADS = 8;

ofstream res_stream;

void init_data();
void init_matrix();
void calculate_serial();
double calculate_MPI_block();
double calculate_MPI_cycle();
double calculate_MPI_pipeline();
double calculate_MPI_SIMD();
double calculate_MPI_OMP();
double calculate_MPI_OMP_SIMD();
void print_matrix();
void test(int);
void print_result(double);


// ��ʼ������
void init_data() {
    origin_data = new float* [N], matrix = new float* [N];
    auto* tmp = new float[N * N];
    for (int i = 0; i < N; i++) {
        origin_data[i] = new float[N], matrix[i] = tmp + i * N;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0, origin_data[i][j] = 0;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            origin_data[i][j] = rand() * 1.0 / RAND_MAX * L;
        }
    }
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int k = 0; k < N; k++) {
                origin_data[j][k] += origin_data[i][k];
            }
        }
    }
}

// ��data��ʼ��matrix����֤ÿ�ν��м����������һ�µ�
void init_matrix() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = origin_data[i][j];
        }
    }
}


int main() {
    MPI_Init(nullptr, nullptr);
#ifdef _TEST
    res_stream.open("result.csv", ios::out);
    LOOP = 5;
    for (int i = 500; i <= 3000; i += 500)
        test(i);
    res_stream.close();
#endif
#ifdef _PRINT
    test(10);
#endif
    MPI_Finalize();
    return 0;
}


// ���Ժ���
void test(int n) {
    N = n;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        cout << "=================================== " << N << " ===================================" << endl;
        res_stream << N;
    }
    struct timeval start {};
    struct timeval end {};
    double time = 0;
    init_data();
    // ====================================== serial ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        init_matrix();
        gettimeofday(&start, nullptr);
        calculate_serial();
        gettimeofday(&end, nullptr);
        time += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000;
    }
    if (rank == 0) {
        cout << "serial:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_block ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_block();
    }
    if (rank == 0) {
        cout << "MPI_block:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_cycle ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_cycle();
    }
    if (rank == 0) {
        cout << "MPI_cycle:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_cycle ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_pipeline();
    }
    if (rank == 0) {
        cout << "MPI_pipeline:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_SIMD ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_SIMD();
    }
    if (rank == 0) {
        cout << "MPI_SIMD:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_OMP ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_OMP();
    }
    if (rank == 0) {
        cout << "MPI_OMP:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    // ====================================== MPI_OMP_SIMD ======================================
    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_OMP_SIMD();
    }
    if (rank == 0) {
        cout << "MPI_OMP_SIMD:" << time / LOOP << "ms" << endl;
        print_result(time);
    }
    if (rank == 0) {
        res_stream << endl;
    }
}



// �����㷨
void calculate_serial() {
    for (int k = 0; k < N; k++) {
        for (int j = k + 1; j < N; j++) {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;
        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}

// MPI_block �����㷨
double calculate_MPI_block() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = ceil(N * 1.0 / size);
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            int start = i * task_num;
            int end = (i + 1) * task_num;
            if (i == size - 1){
                end = N;
            }
            MPI_Send(&matrix[start][0], (end - start) * N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        if (rank != size - 1) {
            MPI_Recv(&matrix[rank * task_num][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&matrix[rank * task_num][0], (N - rank * task_num) * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
    }

    // ����Ԫ����
    int start = rank * task_num;
    int end = (rank + 1) * task_num < N ? (rank + 1) * task_num : N;
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k >= start && k < end) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            for (int p = 0; p < size; p++) {
                if (p != rank) {
                    MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
                }
            }
        }
            // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // ������Ԫ����
        for (int i = max(k + 1, start); i < end; i++) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI_cycle �����㷨
double calculate_MPI_cycle() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto *buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k % size == rank) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            for (int p = 0; p < size; p++) {
                if (p != rank) {
                    MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                }
            }
        }
            // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
        for (int i = begin; i > k; i -= size) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI pipeline �����㷨
double calculate_MPI_pipeline() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto *buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    int pre_proc = (rank + (size - 1)) % size;
    int next_proc = (rank + 1) % size;
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k % size == rank) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            MPI_Send(&matrix[k][0], N, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
        }
            // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, pre_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (next_proc != k % size) {
                MPI_Send(&matrix[k][0], N, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
            }
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
        for (int i = begin; i > k; i -= size) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI SIMD �����㷨
double calculate_MPI_SIMD() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto *buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k % size == rank) {
            // float Akk = matrix[k][k];
            float32x4_t Akk = vmovq_n_f32(matrix[k][k]);
            int j;
            for (j = k + 1; j + 3 < N; j += 4)
            {
                float32x4_t Akj = vld1q_f32(matrix[k] + j);
                Akj = vdivq_f32(Akj, Akk);
                vst1q_f32(matrix[k] + j, Akj);
            }
            // ���д����β
            for (; j < N; j++)
            {
                matrix[k][j] = matrix[k][j] / matrix[k][k];
            }
            matrix[k][k] = 1;
            for (int p = 0; p < size; p++) {
                if (p != rank) {
                    MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                }
            }
        }
            // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
        for (int i = begin; i > k; i -= size) {
            // float Aik = matrix[i][k];
            float32x4_t Aik = vmovq_n_f32(matrix[i][k]);
            int j;
            for (j = k + 1; j + 3 < N; j += 4)
            {
                float32x4_t Akj = vld1q_f32(matrix[k] + j);
                float32x4_t Aij = vld1q_f32(matrix[i] + j);
                float32x4_t AikMulAkj = vmulq_f32(Aik, Akj);
                Aij = vsubq_f32(Aij, AikMulAkj);
                vst1q_f32(matrix[i] + j, Aij);
            }
            // ���д����β
            for (; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI OMP �����㷨
double calculate_MPI_OMP() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto *buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k) shared(matrix, N, size, rank)
    for (k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
#pragma omp single
        {
            if (k % size == rank) {
                for (j = k + 1; j < N; j++) {
                    matrix[k][j] /= matrix[k][k];
                }
                matrix[k][k] = 1;
                for (int p = 0; p < size; p++) {
                    if (p != rank) {
                        MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                    }
                }
            }
                // ������̽��ճ����еĽ��
            else {
                MPI_Recv(&matrix[k][0], N, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
#pragma omp for schedule(simd : guided)
        for (i = begin; i > k; i -= size) {
            for (j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI OMP SIMD �����㷨
double calculate_MPI_OMP_SIMD(){
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto *buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
        // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    int i,j,k;
#pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k) shared(matrix, N, size, rank)
    for (k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
#pragma omp single
        {
            if (k % size == rank) {
                // float Akk = matrix[k][k];
                float32x4_t Akk = vmovq_n_f32(matrix[k][k]);
                for (j = k + 1; j + 3 < N; j += 4)
                {
                    float32x4_t Akj = vld1q_f32(matrix[k] + j);
                    Akj = vdivq_f32(Akj, Akk);
                    vst1q_f32(matrix[k] + j, Akj);
                }
                // ���д����β
                for (; j < N; j++)
                {
                    matrix[k][j] = matrix[k][j] / matrix[k][k];
                }
                matrix[k][k] = 1;
                for (int p = 0; p < size; p++) {
                    if (p != rank) {
                        MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                    }
                }
            }
                // ������̽��ճ����еĽ��
            else {
                MPI_Recv(&matrix[k][0], N, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
#pragma omp for schedule(simd : guided)
        for (i = begin; i > k; i -= size) {
            // float Aik = matrix[i][k];
            float32x4_t Aik = vmovq_n_f32(matrix[i][k]);
            for (j = k + 1; j + 3 < N; j += 4)
            {
                float32x4_t Akj = vld1q_f32(matrix[k] + j);
                float32x4_t Aij = vld1q_f32(matrix[i] + j);
                float32x4_t AikMulAkj = vmulq_f32(Aik, Akj);
                Aij = vsubq_f32(Aij, AikMulAkj);
                vst1q_f32(matrix[i] + j, Aij);
            }
            // ���д����β
            for (; j < N; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// ��ӡ����
void print_matrix() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}



// �����ӡ
void print_result(double time) {
#ifdef _TEST
    res_stream << "," << time / LOOP;
#endif
#ifdef _PRINT
    print_matrix();
#endif
}
