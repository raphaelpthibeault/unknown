#include <stdio.h>

#define KK 100                                  /* the long lag */
#define LL 37                                   /* the short lag */
#define MM (1L << 30)                           /* the modulus */
#define TT 70                                   /* guaranteed separation between streams */
#define mod_diff(x,y) (((x) - (y)) & (MM - 1))  /* (x-y) mod MM */
#define is_odd(x) ((x)&1)                       /* the units bit of x */
#define mod_sum(x,y) (((x)+(y)) - (int)((x)+(y))) /* (x+y) mod 1.0 */

long ran_x[KK];                                 /* the generator state for random 30-bit integers */
double ran_u[KK];                               /* the generator state for random fractions */

/* 30-bit integers */
void
ran_array(long aa[], int n) {                   /* put n new values in aa */
    register int i,j;

    for (j = 0; j < KK; ++j)
        aa[j] = ran_x[j];
    for (; j < n; ++j)
        aa[j] = mod_diff(aa[j-KK], aa[j-LL]);
    for (i = 0; i < LL; ++i, ++j)
        ran_x[i] = mod_diff(aa[j-KK], aa[j-LL]);
    for (; i < KK; ++i, ++j)
        ran_x[i] = mod_diff(aa[j-KK], ran_x[i-LL]);
}

void
ran_start(long seed) {                          /* use this to setup ran_array */
    register int t,j;
    long x[KK+KK-1];                            /* the preparation buffer */
    register long ss = (seed + 2) & (MM - 2);
    for (j = 0; j < KK; ++j) {
        x[j] = ss;
        ss <<=1;                                /* bootstrap the buffer */
        if (ss >= MM) ss -= MM - 2;             /* cyclic shift 29 bits */
    }
    ++x[1];                                     /* make x[1] (and only [x]) odd */
    for (ss = seed&(MM-1), t = TT-1; t; ) {
        for (j = KK-1; j > 0; --j) {            /* "square" */
            x[j+j] = x[j];
            x[j+j-1] = 0;
        }
        for (j = KK+KK-2; j >= KK; --j) {
            x[j-(KK-LL)] = mod_diff(x[j-(KK-LL)], x[j]);
            x[j-KK] = mod_diff(x[j-KK], x[j]);
        }
        if (is_odd(ss)) {                       /* "multiply by z" */
            for (j = KK; j> 0; --j)
                x[j] = x[j-1];
            x[0] = x[KK];                       /* shift the buffer cyclically */
            x[LL] = mod_diff(x[LL], x[KK]);
        }
        if (ss)
            ss >>= 1;
        else
            t--;
    }
    for (j = 0; j < LL; ++j)
        ran_x[j+KK-LL] = x[j];
    for (; j < KK; ++j)
        ran_x[j-LL] = x[j];
    for (j = 0; j < 10; ++j)
        ran_array(x, KK+KK-1);                  /* warm it up */
}


/* double-precision random fractions in the range [0..1) */
void
ranf_array(double aa[], int n) {                /* aa gets n random fractions */
    register int i,j;

    for (j = 0; j < KK; ++j)
        aa[j] = ran_u[j];
    for (; j < n; ++j)
        aa[j] = mod_sum(aa[j-KK], aa[j-LL]);
    for (i = 0; i < LL; ++i, ++j)
        ran_u[i] = mod_sum(aa[j-KK], aa[j-LL]);
    for (; i < KK; ++i, ++j)
        ran_u[i] = mod_sum(aa[j-KK], ran_u[i-LL]);
}


void
ranf_start(long seed) {
    register int t,s,j;
    double u[KK+KK-1];
    double ulp = (1.0/(1L<<30))/(1L<<22);       /* 2^(-52) */
    double ss = 2.0 * ulp * ((seed & 0x3fffffff)+2);

    for (j = 0; j < KK; ++j) {                  /* bootstrap the buffer */
        u[j] = ss;
        ss += ss;
        if (ss >= 1.0) ss -= 1.0 - 2*ulp;       /* cyclic shift of 51 bits */
    }
    u[1] += ulp;                                /* make u[1] (and only u[1]) odd */
    for (s = seed & 0x3fffffff, t = TT-1; t; ) {
        for (j = KK-1; j > 0; --j) {            /* "square" */
            u[j+j] = u[j];
            u[j+j-1] = 0.0;
        }
        for (j = KK+KK-2; j>= KK; --j) {
            u[j-(KK-LL)] = mod_sum(u[j-(KK-LL)], u[j]);
            u[j-KK] = mod_sum(u[j-KK], u[j]);
        }
        if (is_odd(s)) {                        /* multiply by z */
            for (j = KK; j > 0; --j)
                u[j] = u[j-1];
            u[0] = u[KK];                       /* shift the buffer cyclically */
            u[LL] = mod_sum(u[LL], u[KK]);
        }
        if (s)
            s >>= 1;
        else
            --t;
    }

    for (j = 0; j < LL; ++j)
        ran_u[j+KK-LL] = u[j];
    for (; j < KK; ++j)
        ran_u[j-LL] = u[j];
    for (j = 0; j < 10; ++j)
        ranf_array(u, KK+KK-1);                  /* warm everything up */
}

int main() {
    register int m;
    long a[2009];
    double b[2009];

    // integers
    ran_start(310952L);
    for (m = 0; m < 2009; ++m)
        ran_array(a, 1009);
    printf("%ld\n", ran_x[0]);

    ran_start(310952L);
    for (m = 0; m < 1009; ++m)
        ran_array(a, 2009);
    printf("%ld\n", ran_x[0]);

    // fractions
    ranf_start(310952L);
    for (m = 0; m < 2009; ++m)
        ranf_array(b, 1009);
    printf("%.20f\n", ran_u[0]);

    ranf_start(310952L);
    for (m = 0; m < 1009; ++m)
        ranf_array(b, 2009);
    printf("%.20f\n", ran_u[0]);

    return 0;
}
