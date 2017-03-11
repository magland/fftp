#include "fftp.h"
#include <stdio.h>
#include <stdlib.h>
#include <QDebug>
#include <QTime>

int main(int argc,char *argv[]) {
    bigint N=1048576/2;
    //N=16;
    double *input=(double *)malloc(sizeof(double)*N*2);
    double *output=(double *)malloc(sizeof(double)*N*2);

    for (bigint i=0; i<N; i++) {
        input[2*i]=2*(i%2)-1; //real
        input[2*i+1]=0; //imag
    }
    input[2*4]=14.2;

    Fftp X;

    {
        QTime timer; timer.start();
        X.initialize(N);
        double elapsed_sec=timer.elapsed()*1.0/1000;
        double rate=N/elapsed_sec;
        qDebug() << QString("Elapsed for init: %1 sec, %2 pts/sec").arg(elapsed_sec).arg(rate);
    }

    {
        QTime timer; timer.start();
        X.run(input,output);
        double elapsed_sec=timer.elapsed()*1.0/1000;
        double rate=N/elapsed_sec;
        qDebug() << QString("Elapsed for fft: %1 sec, %2 pts/sec").arg(elapsed_sec).arg(rate);
    }


    if (N<=32) {
        X.run(output,input);
        for (bigint i=0; i<N; i++) {
            input[i*2]/=N;
            input[i*2+1]/=N;
        }

        for (bigint i=0; i<N; i++) {
            printf("%lld: %g + %g i ---- %g + %g i\n",i,input[2*i],input[2*i+1],output[2*i],output[2*i+1]);
        }
    }
    free(input);
    free(output);

    return 0;
}
