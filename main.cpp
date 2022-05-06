#include <iostream>
#include<time.h>
#include <pthread.h>
#include<arm_neon.h>
#include<semaphore.h>
using namespace std;
int n=50;
int thread_count=7;
int thread_num=thread_count+1;
float **A;
float **B;
float **C;
float **D;
float **E;
typedef struct
{
    int t_id;//线程id
    int k;
}thread_art;
void m_reset(int n)
{

    A=new float*[n];
    for(int i=0;i<n;i++)
        A[i]=new float[n];
  for(int i=0;i<n;i++)
  {
      for(int j=0;j<i;j++)
        A[i][j]=0;
      A[i][i]=1.0;
      for(int j=i+1;j<n;j++)
      {
         A[i][j]=rand();
      }

  }
  for(int k=0;k<n;k++)
  {

      for(int i=k+1;i<n;i++)
      {
          for(int j=0;j<n;j++)
            A[i][j]+=A[k][j];
      }
  }

}
void deepcopy()
{
    B=new float*[n];
    for(int i=0;i<n;i++)
        B[i]=new float[n];
    C=new float*[n];
    for(int i=0;i<n;i++)
        C[i]=new float[n];
    D=new float*[n];
    for(int i=0;i<n;i++)
        D[i]=new float[n];
    E=new float*[n];
    for(int i=0;i<n;i++)
        E[i]=new float[n];
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            B[i][j]=A[i][j];
            C[i][j]=A[i][j];
            D[i][j]=A[i][j];
            E[i][j]=A[i][j];
        }
    }
}
void chuanxing()
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    for (int k = 0; k < n; k++)
    {
		for (int j = k + 1; j < n; j++)
		{
			C[k][j] /= C[k][k];
		}
		C[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
			{
				C[i][j] -= C[i][k] * C[k][j];
			}
			C[i][k] = 0;
		}

	}


	timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);


}
void *thread_funcd(void *parm)
{
    thread_art*p=(thread_art*)parm;
    int k=p->k;
    int id=p->t_id;
    for(int i=k+id+1;i<n;i+=thread_count)
    {
        float32x4_t vaik=vdupq_n_f32(E[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&E[k][j]);
                float32x4_t vaij=vld1q_f32(&E[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&E[i][j], vaij);

            }
            while(j<n)
            {
                E[i][j]-=E[k][j]*E[i][k];
                j++;
            }
            E[i][k]=0;

    }
    pthread_exit(nullptr);
}
void dynamic()
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    int k;
    for (k = 0; k<n; k++)
    {
		for (int j = k + 1; j < n; j++)
		{
			E[k][j] /= E[k][k];
		}
		E[k][k] = 1.0;
        thread_art*param=new thread_art[thread_count];
        pthread_t*handles=new pthread_t[thread_count];
        for(int i=0;i<thread_count;i++)
        {
            param[i].t_id=i;
            param[i].k=k;

        }
        for(int i=0;i<thread_count;i++)
        {
            pthread_create(&handles[i],nullptr,thread_funcd,(void*)&param[i]);
        }
        for(int i=0;i<thread_count;i++)
        {
            pthread_join(handles[i],nullptr);
        }

	}

	timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
}
sem_t sem_main;
sem_t *sem_workerstart=new sem_t[thread_count];
sem_t *sem_workerend=new sem_t[thread_count];
//sem_t sem_workerstart;
//sem_t sem_workerend;
void *thread_func(void *parm)
{
    thread_art*p=(thread_art*)parm;
    int id=p->t_id;
    for(int k=0;k<n;++k)
    {
        sem_wait(&sem_workerstart[id]);
        for(int i=k+1+id;i<n;i+=thread_count)
        {
           float32x4_t vaik=vdupq_n_f32(A[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&A[k][j]);
                float32x4_t vaij=vld1q_f32(&A[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&A[i][j], vaij);

            }
            while(j<n)
            {
                A[i][j]-=A[k][j]*A[i][k];
                j++;
            }
            A[i][k]=0;

        }
        sem_post(&sem_main);
        sem_wait(&sem_workerend[id]);
    }
    pthread_exit(nullptr);
}
void signal()
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    int flag=sem_init(&sem_main,0,0);
    if(flag!=0)
    {
        cout<<"error";
        return;
    }
    for(int i=0;i<thread_count;i++)
    {
        int flag1=sem_init(&sem_workerstart[i],0,0);
        int flag2=sem_init(&sem_workerend[i],0,0);
        if(flag1==-1||flag2==-1)
        {
            cout<<"error";
            return;
        }
    }
    //sem_init(&sem_workerend,0,0);
    pthread_t *handles=new pthread_t[thread_count];
    thread_art *param=new thread_art[thread_count];
    for(int i=0;i<thread_count;i++)
    {
        param[i].t_id=i;
        pthread_create(&handles[i],nullptr,thread_func,(void*)&param[i]);
    }
    for(int k=0;k<n;++k)
    {
        for(int j=k+1;j<n;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;
        for(int t_id=0;t_id<thread_count;t_id++)
        {
            //cout<<"开始第"<<k<<"轮第"<<t_id<<"消去"<<endl;
            sem_post(&sem_workerstart[t_id]);
            //sem_wait(&sem_main);
        }
        for(int i=0;i<thread_count;i++)
            sem_wait(&sem_main);
        for(int i=0;i<thread_count;i++)
        {
            //cout<<"开始下一轮"<<endl;
            sem_post(&sem_workerend[i]);
        }
    }
    for(int i=0;i<thread_count;i++)
        pthread_join(handles[i],nullptr);
    //sem_destroy(&sem_workerend);
    for(int i=0;i<thread_count;i++)
        sem_destroy(&sem_workerstart[i]);
    sem_destroy(&sem_main);
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
}
sem_t sem_leader;
sem_t *sem_Division=new sem_t[thread_num-1];
sem_t *sem_Elimination=new sem_t[thread_num-1];
void *thread_func2(void *parm)
{
    thread_art*p=(thread_art*)parm;
    int id=p->t_id;
    for(int k=0;k<n;++k)
    {
        if(id==0)
        {
            for(int j=k+1;j<n;j++)
                B[k][j]/=B[k][k];
            B[k][k]=1.0;
        }
        else
            sem_wait(&sem_Division[id-1]);
        if(id==0)
        {
            for(int i=0;i<thread_num-1;i++)
                sem_post(&sem_Division[i]);
        }

            for(int i=k+1+id;i<n;i+=thread_num)
            {
                float32x4_t vaik=vdupq_n_f32(B[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&B[k][j]);
                float32x4_t vaij=vld1q_f32(&B[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&B[i][j], vaij);

            }
            while(j<n)
            {
                B[i][j]-=B[k][j]*B[i][k];
                j++;
            }
            B[i][k]=0;

            }

        if(id==0)
        {
            for(int i=0;i<thread_num-1;i++)
                sem_wait(&sem_leader);
            for(int i=0;i<thread_num-1;i++)
                sem_post(&sem_Elimination[i]);
        }
        else
        {
            sem_post(&sem_leader);
            sem_wait(&sem_Elimination[id-1]);

        }
    }
    pthread_exit(nullptr);
}

void signal_xunhuan()
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    sem_init(&sem_leader,0,0);
    for(int i=0;i<thread_num-1;i++)
    {
        sem_init(&sem_Division[i],0,0);
        sem_init(&sem_Elimination[i],0,0);
    }
    pthread_t *handles=new pthread_t[thread_num];
    thread_art *param=new thread_art[thread_num];
    for(int i=0;i<thread_num;i++)
    {
        param[i].t_id=i;
        pthread_create(&handles[i],nullptr,thread_func2,(void*)&param[i]);
    }
    for(int i=0;i<thread_num;i++)
        {
            pthread_join(handles[i],nullptr);

        }

    for(int i=0;i<thread_num-1;i++)
    {
        sem_destroy(&sem_workerend[i]);
        sem_destroy(&sem_workerstart[i]);
    }
    sem_destroy(&sem_main);
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
}


pthread_barrier_t Divsion;
pthread_barrier_t Elimination;
void *thread_func3(void *parm)
{
    thread_art*p=(thread_art*)parm;
    int id=p->t_id;
    for(int k=0;k<n;++k)
    {
        if(id==0)
        {
            for(int j=k+1;j<n;j++)
                D[k][j]/=D[k][k];
            D[k][k]=1.0;
        }
        pthread_barrier_wait(&Divsion);
        for(int i=k+1+id;i<n;i+=thread_num)
        {
            float32x4_t vaik=vdupq_n_f32(D[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&D[k][j]);
                float32x4_t vaij=vld1q_f32(&D[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&D[i][j], vaij);

            }
            while(j<n)
            {
                D[i][j]-=D[k][j]*D[i][k];
                j++;
            }
            D[i][k]=0;

        }
        pthread_barrier_wait(&Elimination);

    }
    pthread_exit(nullptr);
}
void pthre_barrier()
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    pthread_barrier_init(&Divsion, nullptr, thread_num);
    pthread_barrier_init(&Elimination, nullptr, thread_num);

    pthread_t *handles=new pthread_t[thread_num];
    thread_art *param=new thread_art[thread_num];
    for(int i=0;i<thread_num;i++)
    {
        param[i].t_id=i;
        pthread_create(&handles[i],nullptr,thread_func3,(void*)&param[i]);
    }
    for(int i=0;i<thread_num;i++)
        {
            pthread_join(handles[i],nullptr);

        }

    pthread_barrier_destroy(&Divsion);
    pthread_barrier_destroy(&Elimination);
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
}
int main()
{
    while(n<2000)
    {
        cout<<n<<"         ";
        m_reset(n);
        deepcopy();
        chuanxing();
        dynamic();
        signal();
        signal_xunhuan();
        pthre_barrier();
        cout<<endl;
        n+=100;
    }
    return 0;
}
