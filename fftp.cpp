#include "fftp.h"

#include <QDebug>
#include <QTime>

struct RecombineStep {
    bigint N=0,M=0;
    bigint *indices1=0;
    bigint *indices2=0;
    double *factors_re=0;
    double *factors_im=0;

    void create(bigint N0,bigint M0);
    void apply(double *input,double *output);
    void destroy();
};

class FftpPrivate {
public:
    Fftp *q;

    bigint m_N=0;
    bigint *m_rearrange_indices=0;
    QList<RecombineStep> m_recombine_steps;

    bigint rearrange_index(bigint N,bigint i);
};

Fftp::Fftp()
{
    d=new FftpPrivate;
    d->q=this;
}

Fftp::~Fftp()
{
    if (d->m_rearrange_indices)
        free(d->m_rearrange_indices);
    delete d;
}

bool is_power_of_two(bigint N) {
    if (N==0) return false;
    if (N==1) return true;
    if (N==2) return true;
    if ((N%2)==1) return false;
    return is_power_of_two(N/2);
}

void Fftp::initialize(bigint N)
{
    if (!is_power_of_two(N)) {
        printf("Cannot run Fftp of length that is not a power of two. Aborting.\n");
        exit(-1);
    }
    d->m_N=N;
    d->m_rearrange_indices=(bigint *)malloc(sizeof(bigint)*N);
    for (bigint i=0; i<N; i++) {
        //inefficient?
        d->m_rearrange_indices[i]=d->rearrange_index(N,i);
    }
    d->m_recombine_steps.clear();
    for (bigint M=2; M<=N; M*=2) {
        RecombineStep S;
        d->m_recombine_steps << S;
        d->m_recombine_steps[d->m_recombine_steps.count()-1].create(N,M);
    }
}

void Fftp::run(double *input, double *output)
{
    bigint N=d->m_N;
    QTime timer_tot; timer_tot.start();

    double *working=(double *)malloc(sizeof(double)*N*2);

    //rearrange
    {
        QTime timer; timer.start();
        //cache-inefficient?
        for (bigint i=0; i<N; i++) {
            bigint i1=d->m_rearrange_indices[i];
            working[i1*2]=input[i*2]; //real
            working[i1*2+1]=input[i*2+1]; //imag
        }
        qDebug() << QString("Rearrange: %1 msec (total elapsed: %4 msec)").arg(timer.elapsed()).arg(timer_tot.elapsed());
    }
    double *in=working;
    double *out=output;
    for (int i=0; i<d->m_recombine_steps.count(); i++) {
        QTime timer; timer.start();
        d->m_recombine_steps[i].apply(in,out);
        qDebug() << QString("Recombine step %1 of %2: %3 msec (total elapsed: %4 msec)").arg(i+1).arg(d->m_recombine_steps.count()).arg(timer.elapsed()).arg(timer_tot.elapsed());
        //swap the input/output pointers
        double *hold=in;
        in=out; out=hold;
    }
    //move data to output if necessary
    if (in!=output) {
        for (bigint i=0; i<N; i++) {
            output[i*2]=in[i*2]; //real
            output[i*2+1]=in[i*2+1]; //imag
        }
    }

    free(working);
}

void RecombineStep::create(bigint N0, bigint M0)
{
    N=N0;
    M=M0;
    double *cos_thetas=(double *)malloc(sizeof(double)*(M/2));
    double *sin_thetas=(double *)malloc(sizeof(double)*(M/2));
    //the following could be sped up using multiplications instead of cos/sin
    for (bigint k=0; k<M/2; k++) {
        double theta=-2*M_PI*k*1.0/M;
        cos_thetas[k]=cos(theta);
        sin_thetas[k]=sin(theta);
    }

    indices1=(bigint*)malloc(sizeof(bigint)*N);
    indices2=(bigint*)malloc(sizeof(bigint)*N);
    factors_re=(double*)malloc(sizeof(double)*N);
    factors_im=(double*)malloc(sizeof(double)*N);
    for (bigint i=0; i<N; i+=M) {
        for (bigint k=0; k<M/2; k++) {
            double cos_theta=cos_thetas[k];
            double sin_theta=sin_thetas[k];

            indices1[i+k]=i+k;
            indices2[i+k]=i+k+M/2;
            factors_re[i+k]=cos_theta;
            factors_im[i+k]=sin_thetas[k];

            indices1[i+k+M/2]=i+k;
            indices2[i+k+M/2]=i+k+M/2;
            factors_re[i+k+M/2]=-cos_theta;
            factors_im[i+k+M/2]=-sin_theta;
        }
    }
    free(cos_thetas);
    free(sin_thetas);
}

void RecombineStep::apply(double *input, double *output)
{
    for (bigint i=0; i<N; i++) {
        bigint i1=indices1[i];
        bigint i2=indices2[i];
        double re=factors_re[i];
        double im=factors_im[i];
        //bigint i1=0,i2=0;
        //double re=1,im=0;
        output[2*i]=input[2*i1]+input[2*i2]*re-input[2*i2+1]*im; //real
        output[2*i+1]=input[2*i1+1]+input[2*i2]*im+input[2*i2+1]*re; //imag
        //output[2*i]=input[2*i1]+input[2*i2]*re; //real
        //output[2*i+1]=input[2*i1+1]+input[2*i2]*im; //imag
    }
}

void RecombineStep::destroy()
{
    free(this->indices1);
    free(this->indices2);
    free(this->factors_re);
    free(this->factors_im);
}

bigint FftpPrivate::rearrange_index(bigint N, bigint i)
{
    //inefficient?
    bigint ret=0;
    bigint M=N/2;
    while (i>0) {
        if (i%2==1)
            ret+=M;
        i/=2;
        M/=2;
    }
    return ret;
}
