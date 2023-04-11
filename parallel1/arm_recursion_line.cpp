#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <cmath>
using namespace std;
int main()
{
    int n,i,j;
        int step=2;
        float time_use=0;
        struct timeval start;
        struct timeval end;
        for(n=2;n<=pow(2.0,25.0);n*=step)
        {
                int *a=new int[n];
                int *b=new int[n];
                for(i=0;i<n;i++)
                        a[i]=i;

                int count=100;
                gettimeofday(&start,NULL);
                while(count!=0){
                        count--;
                        for(j=n;j>1;j/=2) // log(n)¸ö²½Öè
                                for(i=0;i<j/2;i++)
                                        if (j==n)
                                             b[i]=a[i*2]+a[i*2+1];
                                        else
                                             b[i]=b[i*2]+b[i*2+1];
                }
                gettimeofday(&end,NULL);
                double timeuse=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)*0.001;
                cout<<n<<" "<<b[0]<<" "<<timeuse<<"ms"<<endl;
        }
        return 0;
}
