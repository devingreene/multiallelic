#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/ioctl.h>
#include<termios.h>
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

struct params {
  uint number_of_loci;
  uint number_of_alleles;
  double rate_of_mutation;
  double rate_of_recombination;
  uint number_of_generations;
  double threshold;
  uint allow_double_mutations; /* 0 or 1; effect not implemented yet */
};

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

ushort hamming_dst(struct params p, uint m, uint n){
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

double *mk_mutation_tbl(struct params p){
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

double *mk_recombination_tbl(struct params p){
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

double *step(struct params p, double *state, double *scratch,
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

uint check_threshold(struct params p, double *state, uint *target_genotypes,
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

/* Write genotype x into buf as number_of_loci base-nalleles digits, most
 * significant first ("12" == 1*nalleles + 2).  buf must hold
 * number_of_loci + 1 chars.  Returns buf. */
char *genotype_to_str(struct params p, uint x, char *buf){
    assert( p.number_of_alleles <= 10 );

    uint n = p.number_of_loci;
    buf[n] = '\0';
    while(n--){
        buf[n] = '0' + x % p.number_of_alleles;
        x /= p.number_of_alleles;
    }
    return buf;
}

/* Heap-allocated version; the caller frees the result. */
char *base_converter(struct params p, uint x){
    return genotype_to_str(p, x, malloc(p.number_of_loci + 1));
}

void print_report(FILE *stream, struct params p, uint ngtypes, uint time, double *state){
    /* One reusable buffer instead of a malloc per genotype per report line */
    char gstr[p.number_of_loci + 1];

    /* Width of the largest generation number, computed exactly; the old
     * ceil(log10(n)) gave width 0 for n == 1 and is off at powers of 10
     * only by luck of rounding. */
    int width = snprintf(NULL, 0, "%u", p.number_of_generations - 1);

    fprintf(stream, "%*u|", width, time);
    for(uint i = 0; i < ngtypes; i++)
        fprintf(stream, " %s:%.0f", genotype_to_str(p, i, gstr), state[i]);
    fputc('\n', stream);
}

int main(int argc, char **argv){
    char progress = 1;

    /* Check for progress option */
    while(argc--)
      if(strcmp(argv[argc], "-s") == 0
	 || strcmp(argv[argc],"--silent") == 0){
	progress = 0;
	break;
      }

    /* Check if stderr is terminal */
    if(!isatty(fileno(stderr))) progress = 0;

    // TODO Adapt if terminal size changes?
    /* Get winsize */
    struct winsize w;
    if(!progress || ioctl(fileno(stderr), TIOCGWINSZ, &w))
       progress = 0;

    /* Comment tag // yaml means data must be parsed from yaml file */
    struct params p; // ymal

    uint nloci = p.number_of_loci;
    uint nalleles = p.number_of_alleles;
    uint ngtypes = power(nalleles, nloci);

    ulong state_size = ngtypes * sizeof(double);

    /* From initial_state */
    double *state = (double*)malloc(state_size); // yaml
    double *fitness = (double*)malloc(state_size); // yaml

    uint target_genotypes_size = ngtypes * sizeof(uint);
    uint *target_genotypes = (uint*)malloc(target_genotypes_size); // yaml
    uint report_points_size = 0; // TODO
    uint *report_points = (uint*)malloc(report_points_size); // yaml

    double *recombination_tbl = mk_recombination_tbl(p);
    double *mutation_tbl = mk_mutation_tbl(p);
    
    double *scratch = calloc(ngtypes, sizeof(double));

    uint n = p.number_of_generations;

    /* The simulation: We record the current sum and use it to normalize
     * after each step */
    double population_size = sum(state, ngtypes);
    uint i;
    ushort current = 0;

    // Parameters for memstream
    char *buf;
    size_t end;

    FILE *stream = open_memstream(&buf, &end);
    fprintf(stream, "\nReport points:\n");
    for(i = 0; i < n; i++){
        if(i == report_points[0]){
            print_report(stream, p, ngtypes, i, state);
            report_points++;
        }
        if(progress){
            uint bar = 0;
            if(current < w.ws_col*i/n){
                fprintf(stderr, "%c", '\r');
                for(; bar < w.ws_col*i/n; bar++)
                    fprintf(stderr, "%c", '-');
                for(; bar < w.ws_col; bar++)
                    fprintf(stderr, "%c", ' ');
                fflush(stderr);
                current = w.ws_col*i/n;
            }
        }
        uint gtype;
        if((gtype = check_threshold(p, state, target_genotypes,
                        ngtypes, population_size)) != (uint)-1)
            goto passed_threshold;
        state = step(p, state, scratch, fitness, recombination_tbl,
                mutation_tbl, population_size);
        continue;
passed_threshold:
        {
            char gstr[p.number_of_loci + 1];
            printf("%sGenotype %s passed threshold at generation %u.\n",
                    progress?"\n":"", genotype_to_str(p, gtype, gstr), i);
        }
        goto report;
    }
    printf("%s%u generations run, no genotype passed threshold.\n",
            progress?"\n":"",p.number_of_generations);
report:
    fclose(stream);
    printf("%s", buf);
}
