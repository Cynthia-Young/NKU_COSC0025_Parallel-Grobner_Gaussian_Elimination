#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

using namespace std;

int main()
{
    int n,i,j;
	int step=10;
	float time_use=0;
	struct timeval start;
	struct timeval end;
	for(n=10;n<=10000;n+=step)
	{
		int **b=NULL;
		b=new int *[n];
		for(i=0;i<n;i++)
		{
			b[i]=new int[n];
		}
		int *a=new int[n];
		int *sum=new int[n];
		//构造矩阵和向量

		for(i=0;i<n;i++)
			a[i]=i;
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
				b[i][j]=i+j;
		}
		//初始化

		int count=100;
		gettimeofday(&start,NULL);
		while(count!=0){
			count--;
			for(i = 0; i < n; i++)
                sum[i] = 0.0;
            for(j = 0; j < n; j++)
                for(i = 0; i < n; i++)
                sum[i] += b[j][i] * a[j ];
		}
		gettimeofday(&end,NULL);
		double timeuse=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)*0.001;
		cout<<n<<" "<<timeuse<<"ms"<<endl;
		if(n==100)step=100;
		if(n==1000)step=500;
		if(n==5000)step=1000;
	}
	return 0;
}