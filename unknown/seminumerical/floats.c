#include <stdio.h>
#include <stdint.h>

/*
meant to be the implementation of floating point in MMIXX, Knuth's architecture. 
*/

#define MSB(n)				    (1ULL << ((n) - 1))

typedef uint32_t tetra;
typedef uint64_t octa;

static const octa inf_octa =    0x7FF0000000000000; 
static const octa max_float =   0x7FF0000000000000; // same as inf_octa
static const octa min_float =   0x0010000000000000;
static const octa standard_NaN =0x7FF8000000000000;

#define ROUND_OFF		        1		// code in rA: 01 - round towards 0
#define ROUND_UP	        	2		// code in rA: 10 - round towards +inf
#define ROUND_DOWN		        3		// code in rA: 11 - round towards -inf
#define ROUND_NEAR		        4		// code in rA: 00 - and to even in case of ties

#define FLOAT_INEXACT           (1 << 0)	// X_BIT
#define FLOAT_DIV_BY_0          (1 << 1)	// Z_BIT
#define FLOAT_UNDERFLOW	        (1 << 2)	// U_BIT
#define FLOAT_OVERFLOW          (1 << 3)	// O_BIT
#define FLOAT_INVALID_OP		(1 << 4)	// I_BIT
#define FLOAT_TO_FIX_OVERFLOW	(1 << 5)	// W_BIT
#define INT_OVERFLOW			(1 << 6)	// V_BIT
#define INT_DIV_CHECK			(1 << 7)	// D_BIT

#define ZERO_EXPONENT           (-1000)	// zero is assumed to have this exponent


typedef enum {
    zro, num, inf , nan
} ftype;

typedef struct {
    octa f;                     // the normalized fraction part
    int e;                      // the raw exponent
    int s;                      // the sign
} Float;

octa
f_round(octa x, int sign, int rounding) {
    switch(rounding) {
        case ROUND_DOWN:
            if(sign == '-')
                x += 3;
            break;
        case ROUND_UP:
            if(sign != '-')
                x += 3;
            break;
        case ROUND_OFF:
            break;
        case ROUND_NEAR:
            x += (x & 4) ? 2 : 1;
            break;
	}
	return x;
}


octa
fpack(Float f, int rounding, unsigned *exceptions) {
    octa o;
    if (f.e > 0x7FD) {          // 0b11111111101
        f.e = 0x7FF;            // 0b11111111111
        o = 0;
    } else {
        if (f.e < 0) {          // negative exponent => unnormalized
            if (f.e < 54) {
                /* if float is larger than precision, then the shift
                 * is 0, then +1 because of the sticky bit */
                o = 1;
            } else {
                octa oo;
                o = f.f >> -f.e;
                oo = o << -f.e;
                if (oo != f.f)  // sticky bit
                    o |= 1; 
            }
            f.e = 0;

        } else {
            o = f.f;
        }
    }

    /* Round and return the result */
    // use the last 2 bits of inexact detection
    if (o & 0b11) {
        *exceptions |= FLOAT_INEXACT;
    }
    // round
    o = f_round(o, f.s, rounding);
    // shift fraction back (see unpack) and add exponent
    o >>= 2;
    // add 1 to e by the addition because the fraction has bit 2^54 set
    o += (octa)f.e << 52;

    if (o >= max_float)
        *exceptions |= FLOAT_OVERFLOW + FLOAT_INEXACT;
    else if (o < min_float)
        *exceptions |= FLOAT_UNDERFLOW;

    if (f.s == '-')
        o |= MSB(64);

    return 0;
}

tetra
sfpack(Float f, int rounding, unsigned *exceptions) {

}



int main() {
    printf("Hello, world!\n");
    return 0;
}
