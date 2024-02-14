#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<strings.h>
#include<err.h>
#include<assert.h>

#define uint unsigned int
#define ulong unsigned long
#define ushort unsigned short

#define index(ngtypes,m,n) ((m) * ngtypes + n)

#define tri_index(ngtypes,m,n,p)\
    (ngtypes) * ((ngtypes) * ((ngtypes) + 1)/2 - ((ngtypes) - (m)) * ((ngtypes) - (m) + 1)/2 + (n) - (m)) + (p)

#define tri_check_bound(ngtypes,m,n,p) assert(tri_index(ngtypes,m,n,p) < ngtypes * ngtypes * (ngtypes + 1)/2)

#define for_each_pair(gtype1, gtype2, ngtypes) \
    for(gtype1 = 0; gtype1 < ngtypes; gtype1++)\
        for(gtype2 = gtype1; gtype2 < ngtypes; gtype2++)

#define For(gtype) for(gtype = 0; gtype < ngtypes; gtype++)

#define swapp(TYPE,a,b) {\
    TYPE *tmp = a;\
    a = b;\
    b = tmp;\
}

typedef struct _params{
    uint number_of_loci;
    uint number_of_alleles;
    double rate_of_mutation;
    double rate_of_recombination;
    /* Not implemented */
    // ulong population_size;
    uint number_of_generations;
} params;

void print_list(char* label, double *arr, uint length){
    uint i = 0;
    printf("%s: [%g",label,arr[i]);
    for(;i < length;i++)
        printf(" %g",arr[i]);
    printf("]\n");
}

double sum(double *state, uint length){
    double rtvl = 0.;
    do{
        rtvl += state[--length];
    } while(length);
}


uint power(uint n, uint p){
    uint rtvl= 1;
    while(p){
        rtvl *= n;
        --p;
    }
    return rtvl;
}

ushort hamming_dst(params p, uint m, uint n){
    uint nalleles = p.number_of_alleles;
    assert( nalleles > 0 );

    ushort rtvl = 0;
    while(m || n){
        ushort digit1 = m % nalleles;
        ushort digit2 = n % nalleles;
        rtvl += ( digit1 != digit2 );
        m /= nalleles;
        n /= nalleles;
    }
    return rtvl;
}

double *mk_mutation_tbl(params p){
    uint ngtypes = power(p.number_of_alleles, p.number_of_loci);
    double *rtvl = calloc(ngtypes * ngtypes, sizeof(double));

    uint gtype1, gtype2;

    for_each_pair(gtype1, gtype2, ngtypes){
        ushort hd = hamming_dst(p, gtype1, gtype2);
        rtvl[index(ngtypes, gtype1, gtype2)]\
            = rtvl[index(ngtypes, gtype2, gtype1)]\
            = power(p.rate_of_mutation, hd) * power(1. - p.rate_of_mutation, p.number_of_loci - hd);
    }

    return rtvl;
}

double *mk_recombination_tbl(params p){
    uint nloci = p.number_of_loci;
    uint nalleles = p.number_of_alleles;
    assert( nalleles != 0 );
    uint ngtypes = power(nalleles, nloci);
    double *rtvl = calloc(ngtypes * ngtypes * (ngtypes + 1)/2, sizeof(double));

    uint gtype, gtype1, gtype2;

    for_each_pair(gtype1, gtype2, ngtypes){
        for(gtype = 0; gtype < ngtypes; gtype++){
            uint x = gtype1, y = gtype2, z = gtype;
            double entry = 1.;
            /* Think of gtypes as base 'nalleles' multidigit numbers. */
            while(x || y || z){
                uint digit1 = x % nalleles;
                uint digit2 = y % nalleles;
                uint digit = z % nalleles;
                if(digit != digit1 && digit != digit2){
                    entry = 0.;
                    break;
                }
                else if(digit1 != digit2)
                    entry *= 0.5;
                x /= nalleles;
                y /= nalleles;
                z /= nalleles;
            }
            tri_check_bound(ngtypes, gtype1, gtype2, gtype);

            entry *= p.rate_of_recombination;
            entry += gtype1 == gtype || gtype2 == gtype?0.5 * (1 - p.rate_of_mutation):0;
            rtvl[tri_index(ngtypes, gtype1, gtype2, gtype)] = entry;
        }
    }
    return rtvl;
}

double *step(params p, double *state, double *scratch,
        double *fitness,
        double *recombination_tbl,
        double *mutation_tbl){
    uint nloci = p.number_of_loci;
    uint nalleles = p.number_of_alleles;
    uint ngtypes = power(nalleles, nloci);

    uint gtype1, gtype2, gtype3;

    /* Record "population size" */
    double start_sum = sum(state, ngtypes);

    /* First fitness */
    For(gtype1) state[gtype1] *= fitness[gtype1];

    /* Then recombination */
    bzero((void*)scratch, ngtypes * sizeof(double));
    for_each_pair(gtype1, gtype2, ngtypes){
        /* Take local segment of pair gtype1, gtype2 */
        double *ltbl = &recombination_tbl[tri_index(ngtypes,gtype1,gtype2,0)];
        For(gtype3) scratch[gtype3] += ltbl[gtype3] * state[gtype1] * state[gtype2];
    }

    swapp(double, state, scratch);

    /* Finally mutation */
    bzero((void*)scratch, ngtypes * sizeof(double));
    For(gtype1){
        /* See above comment in recombination block */
        double *ltbl = &mutation_tbl[index(ngtypes,gtype1,0)];
        For(gtype2)
            ltbl[gtype2] += ltbl[gtype2] * state[gtype1];
    }

    swapp(double, state, scratch);

    /* Normalize */
    double after_sum = sum(state, ngtypes);
    uint n = ngtypes;
    do{
        state[--n] *= start_sum / after_sum;
    } while(n);
            
    return state;
}

int main(){
    uint nread;
    params p;

    if((nread = read(0, &p, sizeof(p))) != sizeof(p)){
        errx(1, "Failed to read %lu bytes when reading parameters."
                "\nRead %u instead.",
                sizeof(p),
                nread);
    }

    printf("\nC program:\nGot %u loci\n"
            "Got %u alleles\n"
            "Got %g rate of mutation\n"
            "Got %g rate of recombination\n"
            // "Got %lu population size\n"
            "Got %u number of generations\n",
            p.number_of_loci,
            p.number_of_alleles,
            p.rate_of_mutation,
            p.rate_of_recombination,
            // ,p.population_size
            p.number_of_generations
            );

    uint nloci = p.number_of_loci;
    uint nalleles = p.number_of_alleles;
    uint ngtypes = power(nalleles, nloci);

    ulong state_size = ngtypes * sizeof(double);

    double *state = (double*)malloc(state_size);
    if((nread = read(0, state, state_size)) != state_size){
        errx(1, "Failed to read %lu bytes when reading: initial_state."
                "\nRead %u instead.",
                state_size,
                nread);
    }

    print_list("initial_state", state, ngtypes);

    double *fitness = (double*)malloc(state_size);
    if((nread = read(0, fitness, state_size)) != state_size){
        errx(1, "Failed to read %lu bytes when reading: fitness."
                "\nRead %u instead."
                ,state_size
                ,nread);
    }

    /* Piped from parser, so tell it that we're done by closing this
     * end. */
    if(close(0) != 0)
        warnx("Couldn't close pipe from parser.");

    print_list("fitness", fitness, ngtypes);

    double *recombination_tbl = mk_recombination_tbl(p);
    double *mutation_tbl = mk_mutation_tbl(p);
    
    double *scratch = calloc(ngtypes, sizeof(double));

    uint n = p.number_of_generations;
    do{
        state = step(p, state, scratch, fitness,
                recombination_tbl,
                mutation_tbl);
    } while(--n);
}
