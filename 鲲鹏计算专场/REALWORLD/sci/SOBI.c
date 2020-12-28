#include <stdio.h>
#include <mkl.h>
#include <string.h>
#include <mkl_scalapack.h>
#include <math.h>
#include <time.h>

#define INPUT_FILE "DATA.bin"
#define OUTPUT_FILE "data.dat"

const int M = 81626;
const int N = 128;

void load(const char *path, float *data)
{
    FILE *input = fopen(path, "rb");
    float *_data = (float *)mkl_malloc(sizeof(float) * M * N, 4 * sizeof(float));
    fread(_data, sizeof(float), N * M, input);
    mkl_somatcopy('R', 'T', M, N, 1, _data, N, data, M);
    mkl_free(_data);
    fclose(input);
}

void save(const char *path, double *data)
{
    mkl_dimatcopy('R', 'T', N, N, 1, data, N, N);
    FILE *fp = fopen(path, "wb+");
    for (int j = 1; j < N*N; j++) {
        fwrite(data + j, sizeof(double), 1, fp);
    }
    fclose(fp);
}

void check(const int *tau, int tau_size)
{
    //not implemented
    return;
}

int main(int argc, char *argv[])
{
    time_t start = time(NULL);
    float *eeg_data = (float *)mkl_malloc(sizeof(float) * M * N, 4 * sizeof(float));
    load(INPUT_FILE, eeg_data);
    //after load eeg_data : 2d array of N x M

    const int tau[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                        12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50,
                        55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
                        125, 150, 175, 200, 225, 250, 275, 300, 325, 350 };
    const int fs = 1000;
    const double jthresh = 1e-5;
    const int lblk = 1000;
    const int max_tau = 350;
    int Ntau = sizeof(tau) / sizeof(int);

//section 1
    check(tau, Ntau);

    const int t_end = M - lblk - max_tau;

    int row_to_copy = N;
    int col_to_copy = lblk + max_tau;
    const int xeeg_ncol = lblk + max_tau;
    float *xeeg = (float *)mkl_malloc(sizeof(float) * N * xeeg_ncol, sizeof(float));
    //rtau_all : 关联矩阵 Ntau x N x N
    float *_rtau_all = (float *)mkl_calloc(Ntau * N * N, sizeof(float), sizeof(float));

    int tt = 0;
    for (int t = 0; t <= t_end; t += lblk, tt += lblk)
    {
        //copy submatrix   means : [XEEG]=EEG_DATA(:,t0:t0+lblk+taumax-1);
        LAPACKE_slacpy(LAPACK_ROW_MAJOR, 'A', row_to_copy, col_to_copy, eeg_data + t, M, xeeg, col_to_copy);

        for (int i = 0; i < N; i++)
        {
            float mean = 0;
            for (int j = 0; j < xeeg_ncol; j++)
            {
                mean += xeeg[i * xeeg_ncol + j];
            }
            mean /= xeeg_ncol;
            for (int j = 0; j < xeeg_ncol; j++)
                xeeg[i * xeeg_ncol + j] -= mean;
        }

        float *rr = (float *)mkl_calloc(N * N, sizeof(float), 64);
        float *rrt = (float *)mkl_malloc(sizeof(float) * N * N, sizeof(float));
        float *sub_xeeg1 = (float *)mkl_malloc(sizeof(float) * N * lblk, 64);
        float *sub_xeeg2 = (float *)mkl_malloc(sizeof(float) * N * lblk, 64);
        __assume_aligned(rr, 64);
        __assume_aligned(sub_xeeg1, 64);
        __assume_aligned(sub_xeeg2, 64);
        LAPACKE_slacpy(LAPACK_ROW_MAJOR, 'A', row_to_copy, lblk, xeeg, xeeg_ncol, sub_xeeg1, lblk);
        for (int i = 0; i < Ntau; i++)
        {
            int taui = tau[i];
            LAPACKE_slacpy(LAPACK_ROW_MAJOR, 'A', row_to_copy, lblk, xeeg + taui, xeeg_ncol, sub_xeeg2, lblk);
            //means  RR  = (X(:,1:lblk)*X(:,1+taui:lblk+taui)');
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, lblk, 1, sub_xeeg1, lblk, sub_xeeg2, lblk, 0, rr, N);
            mkl_somatcopy('R', 'T', N, N, 1, rr, N, rrt, N);
            //means  Rtau_all(:,:,i)=Rtau_all(:,:,i) + 0.5*(RR+RR');
            cblas_saxpy(N * N, 1, rrt, 1, rr, 1);
            cblas_saxpy(N * N, 0.5, rr, 1, _rtau_all + i * N * N, 1);
        }
        mkl_free(rr);
        mkl_free(rrt);
        mkl_free(sub_xeeg1);
        mkl_free(sub_xeeg2);
    }
    mkl_free(xeeg);

    cblas_sscal(N * N * Ntau, 1 / (double)tt, _rtau_all, 1);

//section 2

    //now Ntau = 41
    Ntau = Ntau - 1;

    float *r0t = (float *)mkl_malloc(N * N * sizeof(float), sizeof(float));
    //means  rtau_all(1,:,:) = 0.5(rtau_all(1,:,:)+rtau_all(1,:,:)^T)
    mkl_somatcopy('R', 'T', N, N, 1, _rtau_all, N, r0t, N);
    cblas_saxpy(N * N, 1, r0t, 1, _rtau_all, 1);
    cblas_sscal(N * N, 0.5, _rtau_all, 1);
    mkl_free(r0t);
    double *rtau_all = (double *)mkl_malloc(sizeof(double) * N * N * (Ntau + 1), sizeof(double));
    for (int i = 0; i < N * N * (Ntau + 1); i++)
        rtau_all[i] = _rtau_all[i];
    mkl_free(_rtau_all);

    double *lam0 = (double *)mkl_malloc(sizeof(double) * N, sizeof(double));
    double *lam1 = (double *)mkl_malloc(sizeof(double) * N, sizeof(double));
    double *s0 = (double *)mkl_malloc(sizeof(double) * N * N, sizeof(double));
    double *s1 = (double *)mkl_malloc(sizeof(double) * N * N, sizeof(double));
    lapack_int *isuppz = (lapack_int *)mkl_malloc(sizeof(lapack_int) * N * 2, sizeof(lapack_int));
    lapack_int discard;
    //means  [S0,lam0]=eig(rtau_all(1,:,:))
    LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', N, rtau_all, N,
                   0, 0, 0, 0, -1.0, &discard, lam0, s0, N, isuppz);
//	double *e = (double *)mkl_malloc(sizeof(double) * (N - 1), sizeof(double));
//	double *temp_tau = (double *)mkl_malloc(sizeof(double) * N, sizeof(double));
//	LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', N, rtau_all, N, lam0, e, temp_tau);
//	LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', N, lam0, e, s0, N);
//	mkl_free(e);
//	mkl_free(temp_tau);
    cblas_dcopy(N, lam0, 1, lam1, -1);
    //means   S1=flip(S0, 2);
    for (int i = 0; i < N; i++)
    {
        cblas_dcopy(N, s0 + i * N, 1, s1 + i * N, -1);
    }

    //means  B=(S1*diag(sqrt(1 ./(lam1+1e-20))))'
    for (int i = 0; i < N; i++)
    {
        lam1[i] = lam1[i] + 1e-20;
    }
    //here lam0 used to store the result
    vdInvSqrt(N, lam1, lam0);
    mkl_free(lam1);
    double *diagOfLam = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));
    vdUnpackI(N, lam0, diagOfLam, N + 1);

    double *B = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));
    // I delay B's transposition to line 162
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, s1, N, diagOfLam, N, 0, B, N);
    mkl_free(diagOfLam);
    mkl_free(lam0);
    mkl_free(s0);
    mkl_free(s1);
    mkl_free(eeg_data);

    //A : N x (N * Ntau)
    double *_A = (double *)mkl_malloc(sizeof(double) * N * N * Ntau, sizeof(double));
    double *temp1 = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));
    double *temp2 = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));
    for (int i = 0; i < Ntau; i++)
    {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1, B, N, rtau_all + N * N * (i + 1), N, 0, temp1, N);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, temp1, N, B, N, 0, temp2, N);
        LAPACKE_dlacpy(LAPACK_ROW_MAJOR, 'A', N, N, temp2, N, _A + N * i, N * Ntau);
    }
    mkl_free(temp1);
    mkl_free(temp2);
    double *B_mult = mkl_malloc(sizeof(double) * N * N, sizeof(double));
    // B_mult = trans(B)
    mkl_domatcopy('R', 'T', N, N, 1, B, N, B_mult, N);
    mkl_free(B);
    mkl_free(rtau_all);

    // note: I think in source code 'i' is spelling error, and after I change 'i' to '1', it runs faster and tests right
    const double b[9] = { 1, 0, 0, 0, 1, 1, 0, -1, 1 };
    const double bt[9] = { 1, 0, 0, 0, 1, -1, 0, 1, 1 };

    const int lda = N * Ntau;
    int encore = 1;
    double *V = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));

    double smax = 0;
    double c, s;
    double final_result[4];

    for (int i = 0; i < N; i++)
        V[i * N + i] = 1;

    double *A = (double *)mkl_malloc(sizeof(double) * N * N * Ntau, sizeof(double));
    for (int j = 0; j < N * Ntau; j++)
    {
        int destj = Ntau * (j % N) + j / N;
        for (int i = 0; i < N; i++)
        {
            A[i * lda + destj] = _A[i * lda + j];
        }
    }

    mkl_free(_A);

    int iter_count = 0;
    while (encore)
    {
        encore = 0;
        smax = 0;
        for (int p = 0; p < N - 1; p++)
        {
            for (int q = p + 1; q < N; q++)
            {
                memset(final_result, 0, 4 * sizeof(double));
                for (int i = 0; i < Ntau; i++)
                {
                    double a = A[p * N * Ntau + p * Ntau + i] - A[q * N * Ntau + q * Ntau + i];
                    double b = A[p * N * Ntau + q * Ntau + i] + A[q * N * Ntau + p * Ntau + i];
                    final_result[0 * 2 + 0] += a * a;
                    final_result[0 * 2 + 1] += a * b;
                    final_result[1 * 2 + 0] += b * a;
                    final_result[1 * 2 + 1] += b * b;
                }
                double ton = final_result[0] - final_result[3];
                double toff = final_result[1] + final_result[2];
                double dist = sqrt((ton*ton) + (toff*toff));
                double theta = 0.5 * atan2(toff, ton + dist);
                c = cos(theta);
                s = sin(theta);

                smax = fmax(s, smax);

                if (fabs(s) > jthresh)
                {
                    encore = 1;

                    for (int i = 0; i < N; i++)
                    {
                        double temp = V[p * N + i];
                        V[p * N + i] = temp * c
                                       + V[q * N + i] * s;
                        V[q * N + i] = temp * (-s)
                                       + V[q * N + i] * c;
                    }
                    int pfrom = p * Ntau, pto = (p + 1) * Ntau;
                    int qfrom = q * Ntau, qto = (q + 1) * Ntau;
                    int base;
                    double tempAip;
                    {
                        for (int i = 0; i < lda; i++)
                        {
                            double temp = A[p*lda + i];
                            A[p*lda + i] = c * temp
                                           + s * A[q*lda + i];
                            A[q*lda + i] = (-s) * temp
                                           + c * A[q*lda + i];
                        }
                        for (int i = 0; i < N; i++)
                        {
                            if(i != q && i != p){
                                for (int ip = pfrom, iq = qfrom; ip < pto; ip++, iq++)
                                {
                                    base = i * lda;
                                    tempAip = A[base + ip];
                                    A[base + ip] = c * tempAip + s * A[base + iq];
                                    A[base + iq] = -s * tempAip + c * A[base + iq];
                                }
                            }
                        }
                    }
                    for (int ip = pfrom, iq = qfrom; ip < pto; ip++, iq++)
                    {
                        base = p * lda;
                        tempAip = A[base + ip];
                        A[base + ip] = c * tempAip + s * A[base + iq];
                        A[base + iq] = -s * tempAip + c * A[base + iq];
                    }
                    for (int ip = pfrom, iq = qfrom; ip < pto; ip++, iq++)
                    {
                        base = q * lda;
                        tempAip = A[base + ip];
                        A[base + ip] = c * tempAip + s * A[base + iq];
                        A[base + iq] = -s * tempAip + c * A[base + iq];
                    }

                }
            }// for : q
        }// for : p
        printf("itercount : %d   eps : %lf\n", iter_count, smax);
        ++iter_count;
    }//while encore
    mkl_free(isuppz);
    mkl_free(A);

    double *unscaled_W = (double *)mkl_calloc(N * N, sizeof(double), sizeof(double));
    double *row = (double *)mkl_malloc(sizeof(double) * N, sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, V, N, B_mult, N, 0, unscaled_W, N);
    double sum;
    double scaling_factor;
    for (int i = 0; i < N; i++)
    {
        vdSqr(N, unscaled_W + i * N, row);
        sum = cblas_dasum(N, row, 1);
        scaling_factor = sqrt(sum);
        cblas_dscal(N, 1 / scaling_factor, unscaled_W + i * N, 1);
    }
    // now unscaled_W scaled
    // output unscaled_W here
    save(OUTPUT_FILE, unscaled_W);

    mkl_free(V);
    mkl_free(unscaled_W);
    mkl_free(row);
    time_t end = time(NULL);
    printf("done : %d s\n", end - start);
}
