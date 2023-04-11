#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <arm_neon.h>
#include<ctime>
//VLD1.8 {D0,D1}, [R4:128]!

using namespace std;

int main()
{
	int N,i,j,k;
	int step=10;
	float time_use=0;
	struct timeval start;
	struct timeval end;

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
			srand (static_cast <unsigned> (time(0)));
			for(j=i+1;j<N;j++)
			{
				m[i][j]=(float)rand()/100;
			}
		}
		for(k=0;k<N;k++)
			for(i=k+1;i<N;i++)
				for(j=0;j<N;j++)
					m[i][j]+=m[k][j];
		//测试用例生成

		//cout<<&m[0][0]<<endl; 字节寻址 低位为0代表对齐

		/*if(N<50)
		{
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
				cout<<m[i][j]<<" ";
			cout<<endl;
		}
		cout<<endl;
		}*/
		
		int a=N%4;

		gettimeofday(&start,NULL);
		for(k=0;k<N;k++)
		{
			float32x4_t vt=vmovq_n_f32(m[k][k]);
			for(j=k+1;j+4<=N;j+=4)
			{
				float32x4_t va=vld1q_f32(m[k]+j);
				va=vdivq_f32(va,vt);
				vst1q_f32(m[k]+j,va);
			}
			for(j=N-a;j<N;j++)
			{
				m[k][j]=m[k][j]/m[k][k];
			}
			m[k][k]=1.0;
			//除法
			for(i=k+1;i<N;i++)
			{
				float32x4_t vaik=vmovq_n_f32(m[i][k]);
				for(j=k+1;j+4<=N;j+=4)
				{
					float32x4_t vakj=vld1q_f32(m[k]+j);
					float32x4_t vaij=vld1q_f32(m[i]+j);
					float32x4_t vx=vmulq_f32(vakj,vaik);
					vaij=vsubq_f32(vaij,vx);
					vst1q_f32(m[i]+j,vaij);
				}
				for(j=N-((N-k-1)%4);j<N;j++)
				{
					m[i][j]=m[i][j]-m[k][j]*m[i][k];
				}
				m[i][k]=0;
			}
			//消去
		}
		//实现普通高斯消去并行NEON

		gettimeofday(&end,NULL);
		double timeuse=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)*0.001;
		cout<<N<<" "<<timeuse<<endl;
		//计时:单位ms
		/*if(N<50)
		{
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
				cout<<m[i][j]<<" ";
			cout<<endl;
		}
		cout<<endl;
		}*/
		if(N==100) step=100;
		if(N==1000) step=1000;
		//测试10，20,...,100,200,...,1000,2000,...,5000
	}

}
