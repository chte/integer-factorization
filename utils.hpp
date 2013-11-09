#include <gmpxx.h>
#include <stdint.h>
#include <cmath>

#define ln2 log(2);

// Get the i'th bit in row.
inline int32_t check_bit(uint32_t i, uint64_t *row) {
        return (row[i / sizeof(uint64_t)] & (1 << (i % sizeof(uint64_t)))) != 0;
}

// Set the i'th bit in row to 1.
inline void set_bit(uint32_t i, uint64_t *row) {
        row[i / sizeof(uint64_t)] |= (1 << (i % sizeof(uint64_t)));
}

// Set the i'th bit in row to 0.
inline void clear_bit(uint32_t i, uint64_t *row) {
        row[i / sizeof(uint64_t)] &= ~(1 << (i % sizeof(uint64_t)));
}

// Toggle the i'th bit in row.
inline void flip_bit(uint32_t i, uint64_t *row) {
        row[i / sizeof(uint64_t)] ^= (1 << (i % sizeof(uint64_t)));
}

/* Get GCD for y and z */
/* property of Carl E Data */
void gcd(mpz_class& x, mpz_class& y, mpz_class& z){
    mpz_gcd(x.get_mpz_t(),y.get_mpz_t(),z.get_mpz_t());
}

/* Milley-Rabin prime test */
bool is_prime(mpz_class& x){
    return mpz_probab_prime_p(x.get_mpz_t(),25);
}

// Modular exponentiation using the right-to-left binary method.
uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t m) {
        uint64_t r = 1;
        while(b > 0) {
                if(b & 1)
                        r = r * a % m;
                b >>= 1;
                a = a * a % m;
        }
        return r;
}

int32_t legendre_symbol(uint32_t a, uint32_t p) {
        unsigned long t = pow_mod(a, (p - 1) / 2, p);
        return t > 1 ? -1 : t;
}

// Solve the congruence x^2 = n (mod p).
void shanks_tonelli(uint32_t n, uint32_t p, uint32_t *R) {
        if(p == 2) {
                R[0] = n;
                R[1] = n;
                return;
        }

        uint64_t s = 0;
        uint64_t q = p - 1;

        while(q % 2 == 0) {
                q /= 2;
                ++s;
        }

        uint64_t z = 2;
        while(legendre_symbol(z, p) != -1)
                ++z;

        uint64_t c = pow_mod(z, q, p);              // c is used for successive powers of n to update
                                                    // both z and t

        uint64_t r = pow_mod(n, (q + 1) / 2, p);    // r is a guess of the square root that gets better
                                                    // with each iteration.

        uint64_t t = pow_mod(n, q, p);              // t is the "fudge factor" - by how much we're off
                                                    // with the guess. The invariant x^2 = ab (mod p)
                                                    // is maintained throughout the loop.

        uint64_t m = s;                             // M is the exponent - decreases with each update

        while(t % p != 1) {
                int32_t i = 1;
                while(pow_mod(t, pow(2, i), p) != 1)
                        ++i;
        
                uint64_t b = pow_mod(c,  pow(2, m - i - 1), p);
                r = r * b % p;
                t = t * b * b % p;
                c = b * b % p;
                m = i;
        }

        R[0] = r;
        R[1] = p - r;
}