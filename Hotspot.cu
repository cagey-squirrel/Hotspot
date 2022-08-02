#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* maximum power density possible (say 300W for a 10mm x 10mm chip)	*/
#define MAX_PD (3.0e6)
/* required precision in degrees	*/
#define PRECISION 0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
/* capacitance fitting factor	*/
#define FACTOR_CHIP 0.5

#define BLOCK_DIM_y 16
#define BLOCK_DIM_x 16

/* chip parameters	*/
const float t_chip = 0.0005;
const float chip_height = 0.016;
const float chip_width = 0.016;

/* ambient temperature, outside of box range*/
const float amb_temp = 80.0;

int num_omp_threads;

void single_iteration(float *result, float *temp, float *power, int row, int col, float Cap_1, float Rx_1, float Ry_1, float Rz_1, float step) {

    for (int r = 0; r < row; r++) {
        for (int c = 0; c < col; c++) {
            float delta;
            // corner cases
            if ((r == 0) && (c == 0)) {
                /* Corner 1 */
                delta = (Cap_1) * (power[0] + (temp[1] - temp[0]) * Rx_1 + (temp[col] - temp[0]) * Ry_1 + (amb_temp - temp[0]) * Rz_1);
                //printf("%d: Corner 1 case delta = %f \n", r*col+c, delta);
            } else if ((r == 0) && (c == col - 1)) {
                /* Corner 2 */
                delta = (Cap_1) * (power[c] + (temp[c - 1] - temp[c]) * Rx_1 + (temp[c + col] - temp[c]) * Ry_1 + (amb_temp - temp[c]) * Rz_1);
                //printf("%d: Corner 2 case delta = %f \n", r*col+c, delta);
            } else if ((r == row - 1) && (c == col - 1)) {
                /* Corner 3 */
                delta = (Cap_1) * (power[r * col + c] + (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 + (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 +
                                   (amb_temp - temp[r * col + c]) * Rz_1);
                //printf("%d: Corner 3 case delta = %f \n", r*col+c, delta);
            } else if ((r == row - 1) && (c == 0)) {
                /* Corner 4	*/
                delta = (Cap_1) * (power[r * col] + (temp[r * col + 1] - temp[r * col]) * Rx_1 + (temp[(r - 1) * col] - temp[r * col]) * Ry_1 + (amb_temp - temp[r * col]) * Rz_1);
                //printf("%d: Corner 4 case delta = %f \n", r*col+c, delta);
            } else if (r == 0) {
                /* Edge 1 */
                delta = (Cap_1) * (power[c] + (temp[c + 1] + temp[c - 1] - 2.0 * temp[c]) * Rx_1 + (temp[col + c] - temp[c]) * Ry_1 + (amb_temp - temp[c]) * Rz_1);
                //printf("%d: Edge 1 case delta = %f \n", r*col+c, delta);
            } else if (c == col - 1) {
                /* Edge 2 */
                delta = (Cap_1) * (power[r * col + c] + (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] - 2.0 * temp[r * col + c]) * Ry_1 +
                                   (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 + (amb_temp - temp[r * col + c]) * Rz_1);
                //printf("%d: Edge 2 case delta = %f \n", r*col+c, delta);                   
            } else if (r == row - 1) {
                /* Edge 3 */
                delta = (Cap_1) * (power[r * col + c] + (temp[r * col + c + 1] + temp[r * col + c - 1] - 2.0 * temp[r * col + c]) * Rx_1 +
                                   (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 + (amb_temp - temp[r * col + c]) * Rz_1);
                //printf("%d: Edge 3 case delta = %f \n", r*col+c, delta);
            } else if (c == 0) {
                /* Edge 4 */
                delta = (Cap_1) * (power[r * col] + (temp[(r + 1) * col] + temp[(r - 1) * col] - 2.0 * temp[r * col]) * Ry_1 + (temp[r * col + 1] - temp[r * col]) * Rx_1 +
                                   (amb_temp - temp[r * col]) * Rz_1);
                //printf("%d: Edge 4 case delta = %f \n", r*col+c, delta);
            } else {
                // base case
                delta = (Cap_1 * (power[r * col + c] + (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] - 2.f * temp[r * col + c]) * Ry_1 +
                                  (temp[r * col + c + 1] + temp[r * col + c - 1] - 2.f * temp[r * col + c]) * Rx_1 + (amb_temp - temp[r * col + c]) * Rz_1));
                //printf("delta temp = %f \n", power[r * col + c]);
                //printf("%d: Base case delta = %f \n", r*col+c, delta);
            }
            result[r * col + c] += delta;
        }
    }
}



void my_single_iteration(float *result, float *temp, float *power, int rows, int cols, float Cap_1, float Rx_1, float Ry_1, float Rz_1, float step) {

    for (int r = 1; r < rows-1; r++) {
        for (int c = 1; c < cols-1; c++) {
            float delta;

            delta = (Cap_1 * 
                        (power[(r-1) * (cols-2) + c-1] +                                          // conditional - 
                        (temp[(r + 1) * cols + c] + temp[(r - 1) * cols + c] - temp[r * cols + c] - temp[r * cols + c] * (r > 1 && r < rows-2)) * Ry_1 +
                        (temp[r * cols + c + 1] + temp[r * cols + c - 1] - temp[r * cols + c]     - temp[r * cols + c] * (c > 1 && c < cols-2)) * Rx_1 + 
                        (amb_temp - temp[r * cols + c]) * Rz_1)
                    );

            //printf("%d: Delta = %f \n", (r-1)*(cols-2)+(c-1), delta);
            result[r * cols + c] += delta;
        }
    }

   
}

__global__ void hotspot_cuda_kernell(float* result, float* temp, float* power, int w, int h, float Cap_1, float Rx_1, float Ry_1, float Rz_1){

    int index_y = blockIdx.y * (blockDim.y-2) + threadIdx.y;
    int index_x = blockIdx.x * (blockDim.x-2) + threadIdx.x;
    int middle_square_pos = index_y * w + index_x;

    int shared_index_y = threadIdx.y;
    int shared_index_x = threadIdx.x;

    __shared__ float shared_temp[BLOCK_DIM_y][BLOCK_DIM_x];

    if((index_x) < (w) && index_y < (h)) {
        shared_temp[threadIdx.y][threadIdx.x] = temp[middle_square_pos];
    }

    __syncthreads();

    if(index_x < (w-1) && index_y < (h-1)) {
        if((shared_index_x > 0) && (shared_index_x < (blockDim.x - 1)) && (shared_index_y > 0) && (shared_index_y < (blockDim.y - 1))) {
            int middle_square_pos = index_y * w + index_x;
            float delta = (Cap_1 * 
                        (power[(index_y - 1) * (w-2) + index_x-1] +                                                          // conditional - 
                        (shared_temp[shared_index_y+1][shared_index_x] + shared_temp[shared_index_y-1][shared_index_x] - shared_temp[shared_index_y][shared_index_x] - shared_temp[shared_index_y][shared_index_x] * (index_y > 1 && index_y < h-2)) * Ry_1 +
                        (shared_temp[shared_index_y][shared_index_x+1] + shared_temp[shared_index_y][shared_index_x-1] - shared_temp[shared_index_y][shared_index_x] - shared_temp[shared_index_y][shared_index_x] * (index_x > 1 && index_x < w-2)) * Rx_1 + 
                        (amb_temp - shared_temp[shared_index_y][shared_index_x]) * Rz_1)
                    );
            result[middle_square_pos] += delta;
        }
    }

     __syncthreads();
    

}

void compute_tran_temp(float **result, int num_iterations, float **temp, float *power, int row, int col, int run_cuda) {

    float grid_height = chip_height / row;
    float grid_width = chip_width / col;

    float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * grid_width * grid_height;
    float Rx = grid_width / (2.0 * K_SI * t_chip * grid_height);
    float Ry = grid_height / (2.0 * K_SI * t_chip * grid_width);
    float Rz = t_chip / (K_SI * grid_height * grid_width);

    float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    float step = PRECISION / max_slope / 1000.0;

    float Rx_1 = 1.f / Rx;
    float Ry_1 = 1.f / Ry;
    float Rz_1 = 1.f / Rz;
    float Cap_1 = step / Cap;

    if(!run_cuda)
        for (int i = 0; i < num_iterations; i++) {
            single_iteration(*result, *temp, power, row, col, Cap_1, Rx_1, Ry_1, Rz_1, step);
            float *tmp = *temp;
            *temp = *result;
            *result = tmp;
        }
    
    if(run_cuda) {
        int grid_height, grid_width, size, power_size;
        float *r_gpu, *t_gpu, *power_gpu;

        grid_height = (row+BLOCK_DIM_y-3) / (BLOCK_DIM_y-2);
        grid_width = (col+BLOCK_DIM_x-3) / (BLOCK_DIM_x-2);
        dim3 grid_dim(grid_height, grid_width);
        dim3 block_dim(BLOCK_DIM_y, BLOCK_DIM_x);

        size = (row+2) * (col+2) * sizeof(float);
        power_size = (row) * (col) * sizeof(float);
        cudaMalloc(&r_gpu, size);
        cudaMalloc(&t_gpu, size);
        cudaMalloc(&power_gpu, power_size);
        cudaMemcpy(power_gpu, power, power_size, cudaMemcpyHostToDevice);
        cudaMemcpy(t_gpu, *temp, size , cudaMemcpyHostToDevice);

        cudaEvent_t start = cudaEvent_t();
	    cudaEvent_t stop = cudaEvent_t();
	    cudaEventCreate(&start);
	    cudaEventCreate(&stop);
	    cudaEventRecord(start, 0);

        for (int i = 0; i < num_iterations; i++) {
            //my_single_iteration(r, t, power, row+2, col+2, Cap_1, Rx_1, Ry_1, Rz_1, step);
            hotspot_cuda_kernell<<<grid_dim, block_dim>>>(r_gpu, t_gpu, power_gpu, row+2, col+2, Cap_1, Rx_1, Ry_1, Rz_1);
            float *tmp = r_gpu;
            r_gpu = t_gpu;
            t_gpu = tmp;
        }

        cudaEventRecord(stop, 0);
	    cudaEventSynchronize(stop);
	    float elapsed = 0.f;
	    cudaEventElapsedTime(&elapsed, start, stop);
	    printf("Parallel implementation execution time = %f \n", elapsed);
        
        cudaEventDestroy(start);
	    cudaEventDestroy(stop);

        cudaMemcpy(*temp, t_gpu, size, cudaMemcpyDeviceToHost);
        cudaFree(r_gpu);
        cudaFree(t_gpu);
        cudaFree(power_gpu);
    }
}

void fatal(char *s) {
    fprintf(stderr, "error: %s\n", s);
    exit(1);
}

void writeoutput(float *vect, int grid_rows, int grid_cols, char *file) {
    int i, j;
    FILE *fp;
    char str[256];
    if ((fp = fopen(file, "w")) == 0) printf("The file was not opened\n");
    for (i = 0; i < grid_rows; i++) {
        for (j = 0; j < grid_cols; j++) {

            sprintf(str, "%g\n", vect[i * grid_cols + j]);
            fputs(str, fp);
        }
    }
    fclose(fp);
}

void my_writeoutput(float *vect, int grid_rows, int grid_cols, char *file) {
    int i, j;
    FILE *fp;
    char str[256];
    if ((fp = fopen(file, "w")) == 0) printf("The file was not opened\n");
    for (i = 1; i < grid_rows+1; i++) {
        for (j = 1; j < grid_cols+1; j++) {

            sprintf(str, "%g\n", vect[i * (grid_cols+2) + j]);
            fputs(str, fp);
        }
    }
    fclose(fp);
}



void read_input(float *vect, int grid_rows, int grid_cols, char *file) {
    int i;
    FILE *fp;
    char str[256];
    float val;

    fp = fopen(file, "r");
    if (!fp) fatal("file could not be opened for reading");

    for (i = 0; i < grid_rows * grid_cols; i++) {
        fgets(str, 256, fp);
        if (feof(fp)) fatal("not enough lines in file");
        if ((sscanf(str, "%f", &val) != 1)) fatal("invalid file format");
        vect[i] = val;
    }

    fclose(fp);
}

void usage(int argc, char **argv) {
    fprintf(stderr, "Usage: %s <grid_rows> <grid_cols> <sim_time> <no. of threads><temp_file> <power_file>\n", argv[0]);
    fprintf(stderr, "\t<grid_rows>  - number of rows in the grid (positive integer)\n");
    fprintf(stderr, "\t<grid_cols>  - number of columns in the grid (positive integer)\n");
    fprintf(stderr, "\t<sim_time>   - number of iterations\n");
    fprintf(stderr, "\t<no. of threads>   - number of threads\n");
    fprintf(stderr, "\t<temp_file>  - name of the file containing the initial temperature values of each cell\n");
    fprintf(stderr, "\t<power_file> - name of the file containing the dissipated power values of each cell\n");
    fprintf(stderr, "\t<output_file> - name of the output file\n");
    exit(1);
}

int main(int argc, char **argv) {

    int grid_rows, grid_cols, sim_time;
    float *temp, *power, *result, *my_result, *my_temp;
    char *tfile, *pfile, *ofile;
    const double ACCURACY = 0.01;

    /* check validity of inputs	*/
    if (argc != 8) usage(argc, argv);
    if ((grid_rows = atoi(argv[1])) <= 0 || (grid_cols = atoi(argv[2])) <= 0 || (sim_time = atoi(argv[3])) <= 0 || (num_omp_threads = atoi(argv[4])) <= 0) usage(argc, argv);

    /* allocate memory for the temperature and power arrays	*/
    temp = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    power = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    result = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    my_temp = (float *)calloc((grid_rows+2) * (grid_cols+2), sizeof(float));
    
    if (!temp || !power || !result) fatal("unable to allocate memory");

    /* read initial temperatures and input power	*/
    tfile = argv[5];
    pfile = argv[6];
    ofile = argv[7];

    read_input(temp, grid_rows, grid_cols, tfile);
    read_input(power, grid_rows, grid_cols, pfile);

    // Copy temp into temp for cuda 
    for(int i = 0; i < grid_rows; i++) {
        for(int j = 0; j < grid_cols; j++) {
            my_temp[(grid_cols+2)*(i+1) + j+1] = temp[i*grid_cols+j];
        }
    }
    
    cudaEvent_t start = cudaEvent_t();
	cudaEvent_t stop = cudaEvent_t();
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

    compute_tran_temp(&result, sim_time, &temp, power, grid_rows, grid_cols, 0);
    
    cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsed = 0.f;
	cudaEventElapsedTime(&elapsed, start, stop);
	printf("Gold implementation execution time = %f \n", elapsed);

    writeoutput(temp, grid_rows, grid_cols, ofile);

    my_result = (float *)calloc((grid_rows+2) * (grid_cols+2), sizeof(float));
    
    start = cudaEvent_t();
	stop = cudaEvent_t();
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

    compute_tran_temp(&my_result, sim_time, &my_temp, power, grid_rows, grid_cols, 1);

    cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
    elapsed = 0.f;
	cudaEventElapsedTime(&elapsed, start, stop);
	printf("Full parallel implementation execution time = %f \n", elapsed);

    char* my_o_file = (char*) malloc((strlen(ofile) + 3) * sizeof(char));
    *my_o_file = '\0';
    strcat(my_o_file, ofile);
    strcat(my_o_file, "_my");

    my_writeoutput(my_temp, grid_rows, grid_cols, my_o_file);

    int largest_mistake = 0;
    for(int i = 0; i < grid_rows; i++) {
        for(int j = 0; j < grid_cols; j++) {
            if(fabs(temp[i*grid_cols + j] - my_temp[(grid_cols+2)*(i+1) + j+1]) > ACCURACY) {
                printf("mistake at %d %d expected %f got %f \n", i, j, temp[i*grid_cols + j], my_temp[(grid_cols+2)*(i+1) + j+1]);
                int mistake = fabs(temp[i*grid_cols + j] - my_temp[(grid_cols+2)*(i+1) + j+1]);
                if(mistake > largest_mistake) {
                    largest_mistake = mistake;
                }
                printf("\n\n Test FAILED \n\n");
                exit(-1);
            }
        }
    }

    //printf("Largest mistake = %d \n", largest_mistake);

    printf("\n\n Test PASSED \n\n");


    /* cleanup	*/
    free(temp);
    free(my_temp);
    free(power);
    free(result);
    free(my_result);

    return 0;
}