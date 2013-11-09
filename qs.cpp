#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stack>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <time.h>
#include "utils.hpp"

using namespace std;

// How large chunks should we sieve at a time?
const uint32_t SIEVE_STEP = 65536;

// Max amount prime numbers generated
const uint32_t PRIME_BOUND = 200000;

// Linear relations, how many extra smooth numbers to be generated
const uint32_t LINEAR_RELATIONS = 10;

// Node containing smooth numbers
struct smooth_node {
        uint32_t X;
        vector<uint32_t> smooth_numbers;
        
        smooth_node(uint32_t X, vector<uint32_t> smooth_numbers) {
                this->X = X;
                this->smooth_numbers = smooth_numbers;
        }
};

vector<mpz_class> primes;
stack<mpz_class> factors;

void handle_perfect_power(mpz_class& factor);
bool find_small_factors(mpz_class& factor); // Trial division
void era_sieve();   //Sieve Of Eratosthenes
mpz_class quadratic_sieve(mpz_class &N);

int main() {
        srand(4711);

        // Generate primes with Sieve of Eratosthenes
        era_sieve();

        mpz_class N;
        // Read numbers from stdin until end of file.
        while(cin >> N) {
                factors.push(N);

                while(!factors.empty()) {
                        mpz_class factor = factors.top();
                        factors.pop();

                        // Print if factor is prime
                        if(is_prime(factor)) {
                                cout << factor << endl;
                                continue;
                        } else {
                                // Check for small factors, before running quadratic sieve.
                                if(find_small_factors(factor))
                                    continue;

                                // Quadratic sieve doesn't handle perferct powers very well, handle those separately.
                                if(mpz_perfect_power_p(factor.get_mpz_t())) {
                                    handle_perfect_power(factor);

                                } else {
                                    mpz_class f = quadratic_sieve(factor);

                                    factors.push(f);
                                    factors.push(factor / f);
                                }
                        }
                }

                cout << endl;
        }

        return 0;
}




// A quadratic sieve implementation for integers up to 100 bits. N must be composite.
mpz_class quadratic_sieve(mpz_class &N) {
        vector<uint32_t> factor_base;
        vector<uint32_t> log_primes;

        mpz_class sqrt_N = sqrt(N);

        /* Setting up a Factor Base and a Sieving Interval */       
        uint32_t B;
        {   // Self initialize smoothsness bound
                // Approximate natural logarithm of N.
                float log_N = mpz_sizeinbase(N.get_mpz_t(), 2) * ln2;

                // The optimal smoothsness bound is exp((0.5 + o(1)) * sqrt(log(n)*log(log(n)))).
                B = (uint32_t)ceil(exp((0.5+0.068) * sqrt(log_N * log(log_N)))+550); // 2^(-1/2)
                // B = (B<300) ? 300 : B;
        }

        // Criteria 1:  The primes in factor base must be primes such that 
        //              the Legendre symbol (n/p) = 1.
        // 
        // Criteria 2:  These primes must be less than a bound B.
        //
        // Ref: http://www.cs.virginia.edu/crab/QFS_Simple.pdf 

        {   // Make factor base by sieving.
                char *sieve = new char[B + 1];
                memset(sieve, 1, B + 1);
                for(unsigned long p = 2; p <= B; ++p) {
                        if(!sieve[p])
                                continue;

                        if(mpz_kronecker_ui(N.get_mpz_t(), p) == 1){
                                factor_base.push_back(p);
                                log_primes.push_back(log(p)/log(2));
                        }

                        for(unsigned long i = p; i <= B; i += p)
                                sieve[i] = 0;
                }
                delete[] sieve;
        }


        /* Sieving */
        float *log_Q = new float[SIEVE_STEP];

        // The sieve boundary.
        uint32_t lower_x = 0;
        uint32_t upper_x = SIEVE_STEP;

        // Calculate sieve index (where to start the sieve) for each factor base number.
        uint32_t **sieve_start_indexes = new uint32_t*[2];
        sieve_start_indexes[0] = new uint32_t[factor_base.size()];
        sieve_start_indexes[1] = new uint32_t[factor_base.size()];
        for(uint32_t p = 0; p < factor_base.size(); ++p) {
                // Uses Tonelli Shanks to solve the x^2 = n (mod p) and get the start and end indexes for sieving
                // of the factor base. Start index (x) and end index (N-x).
                uint32_t ssidx[2];
                mpz_class N_mod_P = N % mpz_class(factor_base[p]);
                shanks_tonelli(N_mod_P.get_ui(), factor_base[p], ssidx);

                //start index
                mpz_class a = (ssidx[0] - sqrt_N) % factor_base[p];
                a += factor_base[p];
                sieve_start_indexes[0][p] = a.get_ui();

                //end index
                mpz_class b = (ssidx[1] - sqrt_N) % factor_base[p];
                b += factor_base[p];
                sieve_start_indexes[1][p] = b.get_ui();
        }

        float cheat_log = 0;
        uint32_t estimate = 1;

        // We have a smooth number when all of its prime factors are small. Specifally, 
        // we say m is B-smooth if all of its prime factors are <= B
        vector<smooth_node> smooths;

        // Once we have a set of primes for our factor base, we begin to take numbers
        // x from our sieving interval and calculate Q(x), and check to see if it factors
        // completely over our factor base. If it factors, it is said to have smoothness. If
        // it does not, we throw it out, and we go on to the next element of our sieving
        // interval.

        // Sieve in chunk size to save memory. Sieves until there is enough smooths numbers.
        // NOTE TO SELF: If quadratic sieving gets stuck try use more linear relations.
        while(smooths.size() < (factor_base.size() + LINEAR_RELATIONS)) {
                // Generate our log_Q vector for the sieve, containing log approximations that fit in machine words.
                for(uint32_t t = 1; t < SIEVE_STEP; ++t) {
                        // Log estimate is expensive so dont do it for every Q(x).
                        if(estimate <= (t + lower_x)) {
                                mpz_class y = (sqrt_N + t + lower_x) * (sqrt_N + t + lower_x) - N;

                                // Estimation of 2 logarithm
                                cheat_log = mpz_sizeinbase(y.get_mpz_t(), 2);

                                // The higher t gets, the less the logarithm of log_Q[t] changes.
                                estimate = estimate * 1.8 + 1;
                        }
                        log_Q[t] = cheat_log;
                }

                // Perform log sieve.
                for(uint32_t p = 0; p < factor_base.size(); ++p) {
                        float lg = log_primes[p];
                        for(uint32_t t = 0; t < 2; ++t) {
                                while(sieve_start_indexes[t][p] < upper_x) {
                                        log_Q[sieve_start_indexes[t][p] - lower_x] -= lg;
                                        sieve_start_indexes[t][p] += factor_base[p];
                                }
                                // prime number 2 has one modular root so,
                                // do an early break.
                                if(factor_base[p] == 2)
                                        break;
                        }
                }


                {   // Factor using trail division all values of belonging to null space.
                    // I.e. all values whose logarithms were reduced to approximately zero.
                        float threshold = log(factor_base.back()) / ln2;
                        for(uint32_t i = 0; i < SIEVE_STEP; ++i) {
                                if(fabs(log_Q[i]) < threshold) {
                                        mpz_class y = (sqrt_N + i + lower_x) * (sqrt_N + i + lower_x) - N;
                                        vector<uint32_t> smooth_vector;
                                        for(uint32_t p = 0; p < factor_base.size(); ++p) {
                                                while(mpz_divisible_ui_p(y.get_mpz_t(), factor_base[p])) {
                                                        mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), factor_base[p]);
                                                        smooth_vector.push_back(p);
                                                }
                                        }
 
                                        if(y == 1) {
                                                // Found a V with B-smooth numbers.
                                                smooths.push_back(smooth_node(i + lower_x, smooth_vector));
                                                // Break free from loop if enough smooths numbers has been found.
                                                if(smooths.size() >= (factor_base.size() + LINEAR_RELATIONS))
                                                        break;
                                        } else {
                                                // Not B-smooth numbers, skip.
                                        }
                                }
                        }
                }

                lower_x += SIEVE_STEP;
                upper_x += SIEVE_STEP;
        }

        /* Initial matrix A for gauss elimination */
        uint64_t **A = new uint64_t*[factor_base.size()];

        // The amount of words needed to accomodate a row in the augmented A.
        int row_words = (smooths.size() + sizeof(uint64_t)) / sizeof(uint64_t);

        for(uint32_t i = 0; i < factor_base.size(); ++i) {
                A[i] = new uint64_t[row_words];
                memset(A[i], 0, row_words * sizeof(uint64_t));
        }

        for(uint32_t s = 0; s < smooths.size(); ++s) {
                // For each factor in the smooths number, add the factor to the corresponding element in the A.
                for(uint32_t p = 0; p < smooths[s].smooth_numbers.size(); ++p)
                        flip_bit(s, A[smooths[s].smooth_numbers[p]]);
        }

        // Perform gauss elimination. 
        // Augmented A will have dimension factor_base.size() x (smooths.size() + 1).
        {
                uint32_t i = 0, j = 0;
                while(i < factor_base.size() && j < (smooths.size() + 1)) {
                        uint32_t max_i = i;

                        // Find pivot element.
                        for(uint32_t k = i + 1; k < factor_base.size(); ++k) {
                                if(check_bit(j, A[k]) == 1) {
                                        max_i = k;
                                        break;
                                }
                        }
                        if(check_bit(j, A[max_i]) == 1) {
                                swap(A[i], A[max_i]);
                                
                                for(uint32_t u = i + 1; u < factor_base.size(); ++u) {
                                        if(check_bit(j, A[u]) == 1) {
                                                for(int32_t w = 0; w < row_words; ++w)
                                                        A[u][w] ^= A[i][w];
                                        }
                                }
                                ++i;
                        }
                        ++j;
                }
        }

        mpz_class a;
        mpz_class b;

        // A initialize array to hold a copy of A.
        // Back-substitution on will be performed on the copy.
        uint64_t **back_A = new uint64_t*[factor_base.size()];
        for(uint32_t i = 0; i < factor_base.size(); ++i)
                back_A[i] = new uint64_t[row_words];

        uint32_t *x = new uint32_t[smooths.size()];
        uint32_t *power = new uint32_t[factor_base.size()];

        do {    // Search for a non-trivial factor
                // Deep-copy the gauss eliminated A.
                for(uint32_t i = 0; i < factor_base.size(); ++i)
                        memcpy(back_A[i], A[i], row_words * sizeof(uint64_t));

                // Empty x-vector.
                memset(x, 0, smooths.size() * sizeof(uint32_t));

                {  // Get x by back-substituting A that is in row echelon form.
                        int32_t i = factor_base.size() - 1;

                        while(i >= 0) {
                                // Count non-zero elements in current row.
                                int32_t count = 0;
                                int32_t current = -1;
                                for(uint32_t c = 0; c < smooths.size(); ++c) {
                                        count += check_bit(c, back_A[i]);
                                        current = check_bit(c, back_A[i]) ? c : current;
                                }

                                // Empty row, advance to next.
                                if(count == 0) {
                                        --i;
                                        continue;
                                }

                                // If row has more than 1 non-zero element the system is undetermined.
                                // current x can be choosen freely. Avoid trivial solution by not setting x to 0. 
                                x[current] = count > 1 ? rand() % 2 : check_bit(smooths.size(), back_A[i]);

                                for(int32_t u = 0; u <= i; ++u) {
                                        if(check_bit(current, back_A[u]) == 1) {
                                                if(x[current] == 1)
                                                        flip_bit(smooths.size(), back_A[u]);
                                                clear_bit(current, back_A[u]);
                                        }
                                }

                                if(count == 1)
                                        --i;
                        }
                }
                a = 1;
                b = 1;

                // The way to combine the factor base to get our square.
                memset(power, 0, sizeof(uint32_t) * factor_base.size());
                for(uint32_t i = 0; i < smooths.size(); ++i) {
                        if(x[i] == 1) {
                                for(uint32_t p = 0; p < smooths[i].smooth_numbers.size(); ++p)
                                        ++power[smooths[i].smooth_numbers[p]];
                                b *= (smooths[i].X + sqrt_N);
                        }
                }

                for(uint32_t p = 0; p < factor_base.size(); ++p) {
                        for(uint32_t i = 0; i < (power[p] / 2); ++i)
                                a *= factor_base[p];               
                }

        // If   a =  b (mod N) 
        // or   a = -b (mod N) 
        // then solution is trivial, run the loop again to find a new a and b.
        } while(a % N == b % N || a % N == (-b) % N + N);

        //We have the desired identity (a+b)(a-b) = 0 (mod N)
        //Compute GCD of n with the difference of a and b.
        b -= a;

        mpz_class factor;
        gcd(factor, b, N);

        for(uint32_t i = 0; i < factor_base.size(); ++i) {
                delete[] A[i];
                delete[] back_A[i];
        }
        delete[] A;
        delete[] back_A;
        delete[] sieve_start_indexes[0];
        delete[] sieve_start_indexes[1];
        delete[] sieve_start_indexes;
        delete[] power;
        delete[] log_Q;
        delete[] x;

        return factor;
}

void era_sieve(){
        char *sieve = new char[PRIME_BOUND + 1];
        memset(sieve, 0, PRIME_BOUND + 1);
        for(uint32_t i = 2 ; i < ceil(sqrt(PRIME_BOUND)) ; ++i){ 
                if(sieve[i]) 
                    continue;
                for(int j = i*i; j < PRIME_BOUND; j += i){
                        sieve[j] = 1;
                }
        }
        for(uint32_t i = 2 ; i < PRIME_BOUND ; ++i){
                if(!sieve[i]){
                        primes.push_back((mpz_class)i);
                }
        }
}


void handle_perfect_power(mpz_class& factor){
    {
        mpz_class root, rem;

        // Check root remainders up half of the amount of bits in factor.
        uint32_t max = mpz_sizeinbase(factor.get_mpz_t(), 2) / 2;
        for(uint32_t n = 2; n < max; ++n) {
            mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), factor.get_mpz_t(), n);
            if(rem == 0) {
                    // Push the n root factors.
                    for(uint32_t i = 0; i < n; ++i)
                            factors.push(root);
                    break;
            }
        }
    }
}

bool find_small_factors(mpz_class& factor){
    {
        for(uint32_t p = 0; p < primes.size(); ++p) {
                if(mpz_divisible_ui_p(factor.get_mpz_t(), primes[p].get_ui())) {
                        factors.push(primes[p]);
                        factors.push(factor / primes[p]);
                        return true;
                }
        }
        return false;
    }
}

