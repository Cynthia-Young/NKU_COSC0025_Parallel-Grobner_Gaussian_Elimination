//普通高斯 SSE 并行 x86
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2
#include <iostream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<windows.h>
#include <stdlib.h>
#include<iomanip>
using namespace std;


int main()
{
	int N,i,j,k;
	int step=10;
	for(N=10;N<=5000;N+=step)
	{
		float **m=NULL;;
		m=new float*[N];
		for(i=0;i<N;i++)
			m[i]=new float[N];
		for(i=0;i<N;i++)
		{
			for(j=0;j<i;j++)
				m[i][j]=0;
			m[i][i]=1.0;
			for(j=i+1;j<N;j++)
				m[i][j]=rand();
		}
		for(k=0;k<N;k++)
			for(i=k+1;i<N;i++)
				for(j=0;j<N;j++)
					m[i][j]+=m[k][j];
		//测试用例生成

       	float seconds;
        float T;
      	long long head, tail,freq;
		QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        for(int k=0;k<N;k++)
        {
            __m128 vt=_mm_set1_ps(m[k][k]);
            int j=k+1;
            for(j=k+1;j+4<=N;j=j+4)
            {
                __m128 va = _mm_loadu_ps(m[k]+j);
                va=_mm_div_ps(va,vt);
                _mm_storeu_ps(m[k]+j,va);
            }
           for(;j<N;j++)
           {
                m[k][j]==m[k][j]/m[k][k];
           }
           m[k][k]=1.0;
           for(int i=k+1;i<N;i++)
           {
                __m128 vaik= _mm_set1_ps(m[i][k]);
                j=k+1;
                for(j;j+4<=N;j=j+4)
                {
                    __m128 vakj=_mm_loadu_ps(m[k]+j);
                    __m128 vaij=_mm_loadu_ps(m[i]+j);
                    __m128 vx=_mm_mul_ps(vakj,vaik);
                    vaij=_mm_sub_ps(vaij,vx);
                    _mm_storeu_ps(m[i]+j,vaij);
                }
                for(j;j<N;j++)
                {
                    m[i][j]=m[i][j]-m[k][j]*m[i][k];
                }
                m[i][k]=0;
            }
        }
        QueryPerformanceCounter((LARGE_INTEGER *)&tail);
		T=(tail - head) * 1000.0 / freq;
        cout<<N<<":"<<T<<endl;

        if(N==100) step=100;
		if(N==1000) step=1000;
		//测试10，20,...,100,200,...,1000,2000,...,5000
    }
    return 0;
}
