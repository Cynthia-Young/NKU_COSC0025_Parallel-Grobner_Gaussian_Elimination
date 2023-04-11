#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

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
			for(j=i+1;j<N;j++)
				m[i][j]=rand();
		}
		for(k=0;k<N;k++)
			for(i=k+1;i<N;i++)
				for(j=0;j<N;j++)
					m[i][j]+=m[k][j];
		//测试用例生成

		gettimeofday(&start,NULL);
		for(k=0;k<N;k++)
		{
			for(j=k+1;j<N;j++)
				m[k][j]=m[k][j]/m[k][k];
			m[k][k]=1.0;
			for(i=k+1;i<N;i++)
			{
				for(j=k+1;j<N;j++)
					m[i][j]=m[i][j]-m[i][k]*m[k][j];
				m[i][k]=0;
			}
		}
		//实现普通高斯消去

		gettimeofday(&end,NULL);
		double timeuse=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)*0.001;
		cout<<N<<" "<<timeuse<<endl;
		//计时：单位ms

		if(N==100) step=100;
		if(N==1000) step=1000;
		//测试10，20,...,100,200,...,1000,2000,...,5000
	}
}