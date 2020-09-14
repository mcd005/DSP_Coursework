/*
DSP ENGI4151 COURSEWORK 19/20 (6_2_2020)
GROUP: 12
NAMES: EVAN SUTCLIFFE, WILL PANTON, JAMIE MCDONALD
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "fft.h"


#define std_bytelength 4
#define bufferSize 64 //[MegaBytes] e.g RPi has 2000MB of RAM in total
#define NUM_THREADS 1
#define OUTPUT_SIZE  (((bufferSize * 262144 / 4 + (87 - 1)) * 1)/52)

/* function declaration */
float * init_impulse(int, int, char *, float, int);
float * conv(float *, int, float *, int, int, int, int);
float * conv2(float *, int, float *, int);
float * convfft(float *, int, float *, int);
float * resample(float *, int, int, int, int, float *);
float * decimate(float *, int, int);
void separateToSave(FILE *, FILE *, float *, float *, int);
double Bessel(double);
void *thread_main(void *);
/* function declaration */

static float coeff100[87] = {
    5468287.0,
    16772956.0,
    -31999817.0,
    37472136.0,
    -30069695.0,
    10910445.0,
    14025207.0,
    -35801823.0,
    45498985.0,
    -37906106.0,
    14756344.0,
    15811541.0,
    -42703149.0,
    54241342.0,
    -45291330.0,
    17747549.0,
    18424191.0,
    -49486405.0,
    62745110.0,
    -51999677.0,
    20298023.0,
    20905098.0,
    -55830048.0,
    70381819.0,
    -58019333.0,
    22558201.0,
    23066919.0,
    -61212228.0,
    76820609.0,
    -62987436.0,
    24362852.0,
    24746013.0,
    -65407370.0,
    81632317.0,
    -66618558.0,
    25617255.0,
    25918657.0,
    -68160120.0,
    84652087.0,
    -68766293.0,
    26324539.0,
    26448505.0,
    -69315418.0,
    85682076.0,
    -69315418.0,
    26448505.0,
    26324539.0,
    -68766293.0,
    84652087.0,
    -68160120.0,
    25918657.0,
    25617255.0,
    -66618558.0,
    81632317.0,
    -65407370.0,
    24746013.0,
    24362852.0,
    -62987436.0,
    76820609.0,
    -61212228.0,
    23066919.0,
    22558201.0,
    -58019333.0,
    70381819.0,
    -55830048.0,
    20905098.0,
    20298023.0,
    -51999677.0,
    62745110.0,
    -49486405.0,
    18424191.0,
    17747549.0,
    -45291330.0,
    54241342.0,
    -42703149.0,
    15811541.0,
    14756344.0,
    -37906106.0,
    45498985.0,
    -35801823.0,
    14025207.0,
    10910445.0,
    -30069695.0,
    37472136.0,
    -31999817.0,
    16772956.0,
    5468287.0,
}; /*<-------------FIR filter for 1MHz SUM
sampling frequency: 2500000 Hz

fixed point precision: 32 bits

* 0 Hz - 970000 Hz
  gain = 0
  desired attenuation = -15 dB
  actual attenuation = n/a

* 980000 Hz - 1020000 Hz
  gain = 1
  desired ripple = 10 dB
  actual ripple = n/a

* 1030000 Hz - 1250000 Hz
  gain = 0
  desired attenuation = -20 dB
  actual attenuation = n/a

*/
static float coeff105[87] = {
    100315739.0,
    4353367.0,
    -12687487.0,
    23906977.0,
    -34703569.0,
    41553255.0,
    -41683176.0,
    33939762.0,
    -19259563.0,
    581051.0,
    17789981.0,
    -31223516.0,
    35986416.0,
    -30332197.0,
    15134207.0,
    6159516.0,
    -28274657.0,
    45389275.0,
    -52635114.0,
    47446741.0,
    -30383950.0,
    5155690.0,
    22207821.0,
    -44830485.0,
    56732910.0,
    -54416865.0,
    37902280.0,
    -10904997.0,
    -19968686.0,
    46914890.0,
    -63256407.0,
    64012736.0,
    -48806880.0,
    21134131.0,
    12312807.0,
    -43237266.0,
    63793064.0,
    -68614987.0,
    56265149.0,
    -29664362.0,
    -4619248.0,
    38024376.0,
    -62161637.0,
    70946992.0,
    -62161637.0,
    38024376.0,
    -4619248.0,
    -29664362.0,
    56265149.0,
    -68614987.0,
    63793064.0,
    -43237266.0,
    12312807.0,
    21134131.0,
    -48806880.0,
    64012736.0,
    -63256407.0,
    46914890.0,
    -19968686.0,
    -10904997.0,
    37902280.0,
    -54416865.0,
    56732910.0,
    -44830485.0,
    22207821.0,
    5155690.0,
    -30383950.0,
    47446741.0,
    -52635114.0,
    45389275.0,
    -28274657.0,
    6159516.0,
    15134207.0,
    -30332197.0,
    35986416.0,
    -31223516.0,
    17789981.0,
    581051.0,
    -19259563.0,
    33939762.0,
    -41683176.0,
    41553255.0,
    -34703569.0,
    23906977.0,
    -12687487.0,
    4353367.0,
    100315739.0
}; /*<--------FIR filter for 1.05MHz DIFF

sampling frequency: 2500000 Hz

fixed point precision: 32 bits

* 0 Hz - 1020000 Hz
  gain = 0
  desired attenuation = -20 dB
  actual attenuation = n/a

* 1030000 Hz - 1070000 Hz
  gain = 1
  desired ripple = 10 dB
  actual ripple = n/a

* 1080000 Hz - 1250000 Hz
  gain = 0
  desired attenuation = -20 dB
  actual attenuation = n/a

*/
static float coeffLP[37] = {
    15449059.0,
    16351747.0,
    19944009.0,
    19137775.0,
    12104021.0,
    -1697787.0,
    -20954854.0,
    -42160867.0,
    -60069213.0,
    -68489172.0,
    -61658885.0,
    -35706521.0,
    10088116.0,
    72565443.0,
    144846017.0,
    217286877.0,
    279196897.0,
    320880475.0,
    335578967.0,
    320880475.0,
    279196897.0,
    217286877.0,
    144846017.0,
    72565443.0,
    10088116.0,
    -35706521.0,
    -61658885.0,
    -68489172.0,
    -60069213.0,
    -42160867.0,
    -20954854.0,
    -1697787.0,
    12104021.0,
    19137775.0,
    19944009.0,
    16351747.0,
    15449059.0
};
pthread_mutex_t lock;

typedef struct Thread_cache {
    int N;
    bool extraIteration;
    long int r;
    int I;
    int iter;
    int offset;
    FILE* left;
    FILE* right;
} Thread_cache;

int main(int argc, char *argv[]) {

    // check if we have 3 arguments - remember that the first one is the file name
    if (argc != 4) {
        printf("ERROR: not enough input parameters\n");
        printf("USE: %s input right_out left_out\n", argv[0]);
        exit(1);
    }

    int threadCount = NUM_THREADS;
    pthread_t thread[NUM_THREADS - 1];
    struct Thread_cache thread_c[NUM_THREADS];

    FILE *dataIn;
    remove(argv[3]); //empties any existing output file
    remove(argv[2]);
    FILE *a = fopen(argv[3], "w"); //left
    FILE *b = fopen(argv[2], "w"); //right
    fclose(a);
    fclose(b);
    /******************/
    dataIn = fopen(argv[1], "rb");
    /* if (dataIn == NULL) {
        printf("ERROR: %s does not exist\n", argv[2]);
        exit(1);
    } */
    if (dataIn == NULL) {
        printf("File 'test_input.dat' not found - check directory"); // Program exits if the file pointer returns NULL.
        return (1);
    }
    fseek(dataIn, 0L, SEEK_END);
    long int sz = ftell(dataIn); // measure dataIn's length by seeking to end and reading
    fclose(dataIn);

    int I = sz / (bufferSize * 1048576); // no. of iterations for given buffer size

    long int r = sz % (bufferSize * 1048576); //calc. remainder buffer size for final iteration
    I = I + (r > 0);
    int last_thread = threadCount;
    int iter = (I / threadCount);

    if (I < threadCount) {
        iter = 1;
        last_thread = I;
    }

    int offset = 0;
    for (int i = 0; i < last_thread; i++) {
        int extra_count = (i < (I % last_thread)); // +1 if non zero divisor , distributed evenly 
        thread_c[i].I = I;
        thread_c[i].N = bufferSize * 262144;
        thread_c[i].extraIteration = false;
        thread_c[i].iter = iter + extra_count;
        thread_c[i].r = 0;
        thread_c[i].offset = offset;
        thread_c[i].left = fopen(argv[3], "w"); //left
        thread_c[i].right = fopen(argv[2], "w"); //right
        offset += iter + extra_count;
    }
    thread_c[last_thread - 1].extraIteration = (r > 0); //extra chunk not size of N
    thread_c[last_thread - 1].r = r / 4; //extra chunk not size of N
    for (int i = 0; i < last_thread - 1; i++) {

        if (pthread_create(&thread[i], NULL, thread_main, (void*) &thread_c[i])) {
            fprintf(stderr, "Error creating thread\n");
            return 1;

        }
    }
    thread_main(&thread_c[last_thread - 1]);
    for (int i = 0; i < last_thread - 1; i++) {
        pthread_join(thread[i], NULL);
    }
    return 0;
}

void *thread_main(void *x_void_ptr) {

    struct Thread_cache cache = *((struct Thread_cache*) x_void_ptr);
    int N = cache.N; // size of chunk to be read and processed 
    int iter = cache.iter; // number of iteration done by this thread
    int staticN = N;
    int offset = cache.offset;

    FILE *dataIn; // each thread has a separate pointer to the file (allows concurrent access for r/w)
    FILE *left_out_file;
    FILE *right_out_file;
    int L = 1;
    int M = 52;
    int dec = 4;
    int FIRorder = 87;
    int lengh_of_output = (N / dec + (FIRorder - 1)) * L / M; //Length of convolution result (and therefore SUM/DIFF I assume)
    int write_offset = offset * lengh_of_output * sizeof (float);

    left_out_file = cache.left;
    right_out_file = cache.right;
    fseek(left_out_file, write_offset, SEEK_SET);
    fseek(right_out_file, write_offset, SEEK_SET);

    float *input;
    dataIn = fopen("test_input.dat", "r");
    //seek(DataIn, sizeof(float)*offset, SEEK_SET);

    for (int j = 0; j < iter; j++) {
        //printf("%i,%i || Begin \n", i, j);
        fprintf(stderr, "□□");
        long int offset_byte = (offset + j) * N; // No. Samples offset from file beginning (n.b or i * bufferSize * 262144 is declaring after N = r/4;)
        if ((j == iter - 1) && cache.extraIteration) {
            input = (float *) malloc(cache.r);
            N = cache.r;
        }

        input = (float *) malloc(N * 4); //malloc is stored on heap (slow but large)
        memset(input, 0, (N * 4));
        fseek(dataIn, std_bytelength*offset_byte, SEEK_SET);
        fread(input, sizeof (float), N, dataIn);

        //RESAMPLE fs = 10MHz ---> 2.5MHz (decimate by: 4)
        float *input_dec;
        //ignore pre decimation filtering action as next step will filter out these freqs. anyway
        input_dec = decimate(input, N, dec);
        free(input);

        //Demodulating SUM component
        FIRorder = 87;
        float *SUM;
        SUM = conv(input_dec, N / dec, coeff100, FIRorder, 1000000, (offset + j), staticN / dec);

        //Demodulating DIFF component
        float *DIFF;
        DIFF = conv(input_dec, N / dec, coeff105, FIRorder, 1050000, (offset + j), staticN / dec);
        free(input_dec);


        //INSERT DOWN-SAMPLING CODE HERE
        FIRorder = 48;
        int L = 1; //interpolation factor
        int M = 52; //decimation factor
        SUM = resample(SUM, N / dec, L, M, FIRorder, coeffLP);
        DIFF = resample(DIFF, N / dec, L, M, FIRorder, coeffLP);

        //Left.dat/Right.dat SEPARATION AND SAVE CODE
        FIRorder = 87;
        lengh_of_output = ((N / dec + (FIRorder - 1)) * L) / M; //Length of convolution result (and therefore SUM/DIFF I assume)
        separateToSave(left_out_file, right_out_file, SUM, DIFF, lengh_of_output);

        free(SUM);
        free(DIFF);
    }
    fclose(dataIn);
    fclose(left_out_file);
    fclose(right_out_file);

    fprintf(stderr, "■");
    return NULL;

}

float * init_impulse(int theta2, int theta1, char *window, float beta, int FIRorder) {
    int i;
    int theta1n = theta1 * 2 * M_PI / 48000;
    int theta2n = theta2 * 2 * M_PI / 48000;
    float *coeff;
    int M = FIRorder;
    coeff = (float *) malloc((M + 1) * sizeof (float));
    if (theta1 != 0) {
        for (i = 0; i < M; i++) {
            coeff[i] = (theta2n / M_PI)*(sin(theta2n * ((float) i - (float) M / 2))) / (theta2n * ((float) i - (float) M / 2))-
                    (theta1n / M_PI)*(sin(theta1n * ((float) i - (float) M / 2))) / (theta1n * ((float) i - (float) M / 2));
        }
    } else {
        for (i = 0; i < M; i++)
            coeff[i] = (theta2n / M_PI)*(sin(theta2n * ((float) i - (float) M / 2))) / (theta2n * ((float) i - (float) M / 2));
    }

    // replace the NaN with the correct value
    coeff[M / 2] = theta2n / M_PI - theta1n / M_PI;

    //  apply selected window
    if (strcmp(window, "hamming") == 0) {
        for (i = 0; i < M; i++)
            coeff[i] = coeff[i]*(0.54 - 0.46 * cos(2 * M_PI * i / (M - 1)));
    } else if (strcmp(window, "rectangular") == 0) {
    } else if (strcmp(window, "hanning") == 0) {
        for (i = 0; i < M; i++)
            coeff[i] = coeff[i]*(0.5 + 0.5 * cos(2 * M_PI * i / (M - 1)));
    } else if (strcmp(window, "blackman") == 0) {
        for (i = 0; i < M; i++)
            coeff[i] = coeff[i]*(0.42 + 0.5 * cos(2 * M_PI * i / (M - 2)) + 0.08 * cos(4 * M_PI * i / (M - 2)));
    } else if (strcmp(window, "kaiser") == 0) {
        double Arg;
        for (i = 0; i < M; i++) {
            Arg = beta * sqrt(1 - pow((i / (M - 1) / 2), 2));
            coeff[i] = coeff[i]*(Bessel(Arg) / Bessel(beta));
        }
    } else if (strcmp(window, "none") == 0) {
    } else printf("error: Window type: \"%s\" could not be found - no window was applied\n", window);

    return coeff;
}

float * conv(float *input, int N, float *firCoefficients, int FIRorder, int Fnyquist, int iteration, int staticN) {

    int Lr = N + (FIRorder - 1); //Length of convolution result
    int Li = N + 2 * (FIRorder - 1); //Length of zero padded input

    float *x;
    x = (float *) malloc(Li * sizeof (float));

    float *padding;
    padding = (float *) malloc(FIRorder * sizeof (float));
    memset(padding, 0, (FIRorder * sizeof (float)));

    memcpy(x, padding, (FIRorder - 1) * sizeof (float)); // zero padding at beginning
    memcpy(x + (FIRorder - 1), input, N * sizeof (float)); // add input signal
    memcpy(x + (N + (FIRorder - 1)), padding, (FIRorder - 1) * sizeof (float)); // zero padding at end

    float *y;
    y = (float *) malloc(Lr * sizeof (float));

    //perform convolution
    int i, j;
    int phaseOffset = iteration*staticN;
    for (i = 0; i < Lr; i++) {
        for (j = 0; j < FIRorder; j++) {
            y[i] += x[i + j] * firCoefficients[j]; //multiply and accumulate
        }
        y[i] = y[i] * cos(Fnyquist * 8 * M_PI * (i + phaseOffset) / 10000000); //local oscillator nested inside convolution to speed up memory access
    }
    free(x);
    return y;
}

float * conv2(float *input, int N, float *firCoefficients, int FIRorder) {

    int Lr = N + (FIRorder - 1); //Length of convolution result
    int Li = N + 2 * (FIRorder - 1); //Length of zero padded input

    float *x;
    x = (float *) malloc(Li * sizeof (float));

    float *padding;
    padding = (float *) malloc(FIRorder * sizeof (float));
    memset(padding, 0, (FIRorder * sizeof (float)));

    memcpy(x, padding, (FIRorder - 1) * sizeof (float)); // zero padding at beginning
    memcpy(x + (FIRorder - 1), input, N * sizeof (float)); // add input signal
    memcpy(x + (N + (FIRorder - 1)), padding, (FIRorder - 1) * sizeof (float)); // zero padding at end

    float *y;
    y = (float *) malloc(Lr * sizeof (float));

    //perform convolution
    int i, j;
    for (i = 0; i < Lr; i++) {
        for (j = 0; j < FIRorder; j++) {
            y[i] += x[i + j] * firCoefficients[j]; //multiply and accumulate
        }
    }
    free(x);
    return y;

}

double Bessel(double x)// This would be used with the Kaiser window.
{
    double Sum = 0.0, XtoIpower;
    int i, j, Factorial;
    for (i = 1; i < 10; i++) {
        XtoIpower = pow(x / 2.0, (double) i);
        Factorial = 1;
        for (j = 1; j <= i; j++)Factorial *= j;
        Sum += pow(XtoIpower / (double) Factorial, 2.0);
    }
    return (1.0 + Sum);
}

float * resample(float *input, int N, int p, int q, int FIRorder, float *LPFcoefficients) {
    float *resampled;
    int resampledLength = (int) ((float) N * ((float) p / (float) q)); //Ints insteasd??
    resampled = (float *) malloc(resampledLength * sizeof (float));
    memset(resampled, 0, (resampledLength * sizeof (float)));

    float *delayLine;
    int taps_per_subfilter = FIRorder / p;
    delayLine = (float *) malloc(taps_per_subfilter * sizeof (float));
    memset(delayLine, 0, (taps_per_subfilter * sizeof (float)));

    const float *ptr_subcoeffs;

    int phase_num = p;

    int index = 0;
    int index2 = 0;
    while (index < N) {
        while (phase_num >= p) { //This exist so samples are kept shifting through delay line until you get to one you won't discard by decimation
            phase_num -= p;
            for (int tap = taps_per_subfilter - 1; tap >= 1; tap--) {
                delayLine[tap] = delayLine[tap - 1];
            }
            delayLine[0] = input[index]; //Do this with indices?
            index++;
        }
        while (phase_num < p) {
            ptr_subcoeffs = LPFcoefficients + phase_num;

            float sum = 0.0;
            for (int tap = 0; tap < taps_per_subfilter; tap++) {
                sum += *ptr_subcoeffs * delayLine[tap];
                ptr_subcoeffs += p;
            }
            resampled[index2] = sum;
            index2++;
            phase_num += q;
        }
    }

    free(delayLine);
    return resampled;
}

void separateToSave(FILE *left_out_file, FILE *right_out_file, float *SUM, float *DIFF, int len) {
    float *L;
    float *R;
    L = (float *) malloc(len * sizeof (float));
    R = (float *) malloc(len * sizeof (float));
    for (int i = 0; i < len; i++) {
        R[i] = (SUM[i] + DIFF[i]) / 2.0;
        L[i] = SUM[i] - R[i];
        R[i] = R[i]/1320000000;
        L[i] = L[i]/1320000000;
    }
    pthread_mutex_lock(&lock);
    fwrite(R, sizeof (float), len, right_out_file);
    fwrite(L, sizeof (float), len, left_out_file);
    pthread_mutex_unlock(&lock);
    free(L);
    free(R);
}

float * convfft(float *input, int N, float *firCoefficients, int FIRorder) {

    float *coefficients;
    coefficients = (float *) malloc(N * sizeof (float));

    int *padding;
    padding = (int *) malloc((N - FIRorder) * sizeof (int));
    memset(padding, 0, ((N - FIRorder) * sizeof (int))); //we cant rely on malloc initilising values to zero

    memcpy(coefficients, padding, (N - FIRorder) * sizeof (float)); // zero padding at beginning
    memcpy(coefficients, firCoefficients, FIRorder * sizeof (float)); // add firCoeffciencts signal

    float *y;
    y = (float *) malloc(N * sizeof (float));

    Fft_convolveReal(input, coefficients, y, N);

    free(coefficients);
    free(padding);
    return y;
}

float * decimate(float *input, int N, int M) {
    float *decimated;
    int decimatedLength = (int) ((float) N / (float) M);
    decimated = (float *) malloc(decimatedLength * sizeof (float));
    memset(decimated, 0, (decimatedLength * sizeof (float))); //set the whole new array to values of zero

    int index = 0;
    for (int i = 0; i < decimatedLength; i++) {
        decimated[i] = input[index];
        index = index + M;
    }
    return decimated;
}