#include <iostream>
#include <windows.h>
//#include <stdlib .h>

using namespace std;

int main()
{
    int n,i,j;
	int step=10;
	long long head, tail , freq ;
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

		int counti=100;
		QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
		QueryPerformanceCounter((LARGE_INTEGER *)&head);
		while(counti!=0){
			counti--;
			for(i=0;i<n;i++)
			{
				sum[i]=0.0;
				for(j=0;j<n;j++)
					sum[i]+=b[j][i]*a[j];
			}
		}
		QueryPerformanceCounter((LARGE_INTEGER *)&tail );
		cout<<n<<" "<<( tail-head)*1000.0 / freq<<"ms"<<endl;
		if(n==100)step=100;
		if(n==1000)step=500;
		if(n==5000)step=1000;


	}
	return 0;
}
