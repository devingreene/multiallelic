#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>
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
    double threshold;
} params;

double sum(double *state, uint length){
    double rtvl = 0.;
    do{
        rtvl += state[--length];
    } while(length);
    return rtvl;
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
    /* TODO: Chase this restriction up the stack */
    assert( p.number_of_alleles >= 2);
    double nalts = (double)(p.number_of_alleles - 1);

    for_each_pair(gtype1, gtype2, ngtypes){
        ushort hd = hamming_dst(p, gtype1, gtype2);
        rtvl[index(ngtypes, gtype1, gtype2)] =
            rtvl[index(ngtypes, gtype2, gtype1)] =
            pow(p.rate_of_mutation/nalts, hd) *
            pow(1. - p.rate_of_mutation, p.number_of_loci - hd);
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
        For(gtype){
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
                if(digit1 != digit2)
                    entry *= 0.5;
                x /= nalleles;
                y /= nalleles;
                z /= nalleles;
            }
            tri_check_bound(ngtypes, gtype1, gtype2, gtype);

            entry *= p.rate_of_recombination;
            if(gtype1 == gtype2 && gtype2 == gtype)
                entry += (1 - p.rate_of_recombination);
            else if(gtype == gtype1 || gtype == gtype2)
                entry += 0.5 * (1 - p.rate_of_recombination);

            /* Double probability for crosses */
            if(gtype1 != gtype2) entry *= 2.;

            rtvl[tri_index(ngtypes, gtype1, gtype2, gtype)] = entry;
        }
    }
    return rtvl;
}

double *step(params p, double *state, double *scratch,
        double *fitness,
        double *recombination_tbl,
        double *mutation_tbl,
        double population_size){
    uint nloci = p.number_of_loci;
    uint nalleles = p.number_of_alleles;
    uint ngtypes = power(nalleles, nloci);

    uint gtype1, gtype2, gtype3;

    /* First fitness */
    For(gtype1) state[gtype1] *= fitness[gtype1];

    /* Then recombination */
    bzero((void*)scratch, ngtypes * sizeof(double));
    for_each_pair(gtype1, gtype2, ngtypes){
        /* Take local segment of pair gtype1, gtype2 */
        double *ltbl = &recombination_tbl[tri_index(ngtypes,gtype1,gtype2,0)];
        For(gtype3) scratch[gtype3] += ltbl[gtype3] *
            state[gtype1] *
            state[gtype2];
    }

    swapp(double, state, scratch);

    /* Finally mutation */
    bzero((void*)scratch, ngtypes * sizeof(double));
    For(gtype1){
        /* See above comment in recombination block */
        double *ltbl = &mutation_tbl[index(ngtypes,gtype1,0)];
        For(gtype2)
            scratch[gtype2] += ltbl[gtype2] * state[gtype1];
    }

    swapp(double, state, scratch);

    /* Normalize */
    double after_sum = sum(state, ngtypes);

    while(ngtypes--)
        state[ngtypes] *= population_size / after_sum;
            
    return state;
}

uint check_threshold(params p, double *state, uint *target_genotypes,
        uint ngtypes, double population_size){
    uint i = 0;
    for(; i < ngtypes; i++){
        if(target_genotypes[i] == (uint)-1)
            return -1;
        assert(target_genotypes[i] < ngtypes);

        if(state[target_genotypes[i]] > p.threshold * population_size){
            return target_genotypes[i];
        }
    }
    return -1;
}

char *base_converter(params p, uint x){
    assert( p.number_of_alleles <= 10 );

    char *rtvl = (char *)malloc(p.number_of_loci * sizeof(char));
    char s = '0';
    uint n = p.number_of_loci;

    memset((void*)rtvl, (int)s, n);
    while(n--){
        rtvl[n] = s + x % p.number_of_alleles;
        x /= p.number_of_alleles;
    }
    return rtvl;
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

    double *fitness = (double*)malloc(state_size);
    if((nread = read(0, fitness, state_size)) != state_size){
        errx(1, "Failed to read %lu bytes when reading: fitness."
                "\nRead %u instead."
                ,state_size
                ,nread);
    }

    uint target_genotypes_size = ngtypes * sizeof(uint);
    uint *target_genotypes = (uint*)malloc(target_genotypes_size);
    if((nread = read(0, target_genotypes, target_genotypes_size)) != target_genotypes_size){
        errx(1, "Failed to read %u bytes when reading: target_genotypes."
                "\nRead %u instead.",
                target_genotypes_size,
                nread);
    }

    /* Piped from parser, so tell it that we're done by closing this
     * end. */
    if(close(0) != 0)
        warnx("Couldn't close pipe from parser.");

    double *recombination_tbl = mk_recombination_tbl(p);
    double *mutation_tbl = mk_mutation_tbl(p);
    
    double *scratch = calloc(ngtypes, sizeof(double));

    uint n = p.number_of_generations;

    /* The simulation: We record the current sum and use it to normalize
     * after each step */
    double population_size = sum(state, ngtypes);
    while(n--){
        uint gtype;
        if((gtype = check_threshold(p, state, target_genotypes,
                        ngtypes, population_size)) != (uint)-1)
            goto passed_threshold;
        state = step(p, state, scratch, fitness, recombination_tbl,
                mutation_tbl, population_size);
        continue;
passed_threshold:
        printf("Genotype %s passed threshold at generation %d.\n",
                base_converter(p, gtype),
                p.number_of_generations - n - 1);
        return 0;
    }
    printf("%u generations run, no genotype passed threshold.\n",
            p.number_of_generations);
}
