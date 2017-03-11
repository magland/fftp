#ifndef FFTP_H
#define FFTP_H

typedef long long bigint;

class FftpPrivate;
class Fftp
{
public:
    friend class FftpPrivate;
    Fftp();
    virtual ~Fftp();

    void initialize(bigint N);
    void run(double *input,double *output); //length 2xN, complex N-vectors

private:
    FftpPrivate *d;
};

#endif // FFTP_H
