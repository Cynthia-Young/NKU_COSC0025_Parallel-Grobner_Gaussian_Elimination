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
                int sum,sum1,sum2;
                //¹¹ÔìÏòÁ¿

                for(i=0;i<n;i++)
                        a[i]=i;

                int count=100;
                gettimeofday(&start,NULL);
                while(count!=0){
                        count--;
                        sum=0;
                        sum1=0;
                        sum2=0;
                        for(i=0;i<n;i+=2){
                                sum1+=a[i];
                                sum2+=a[i+1];
                        }
                        sum =sum1+sum2;
                }
                gettimeofday(&end,NULL);
                double timeuse=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)*0.001;
                cout<<n<<" "<<sum<<" "<<timeuse<<"ms"<<endl;
        }
        return 0;
}
