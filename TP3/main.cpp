#include <stdio.h>
#include <math.h>
#include <omp.h>

#define PI 4*atan(1.)

double f(double x){
    // return 4./(1 + x*x);
    // return sin(x);
    // return exp(x);
    return x*x;
}

/**
     * # Mid-Point integral computing
     * This function computes the integral of `f` in the range `xMin,xMax` using the
     * mid-point rectangles method.
     * ## Parameters
     * - double f(double x)
     *      Callable function to integrate.
     * - double xMin
     *      Starting point of the integral
     * - double xMax
     *      Finishing point of the integral
     * - double nPoints
     *      Number of points to use in the integral
     * ## Returns
     * - double integral
     *      The integral of f(x) between xMin and xMax: \int_xMin^xMax f(x)dx
     */
double midpointIntegral(double f(double x), double xMin, double xMax, long long int nPoints,int N_THREADS){
    double dx = (xMax - xMin)/nPoints;
    // printf("dx = %E...\n",dx);
    double integral = 0;
    omp_set_num_threads(N_THREADS);
    #pragma omp parallel
    {
        double local_integral = 0;
        #pragma omp for
        for(long long int i = 0; i<nPoints; i++){
            local_integral += f(xMin + i*dx + 0.5*dx);
        }
        #pragma omp critical
        integral += local_integral*dx;
        // printf("nThreads = %d\tt_id = %d\tlocal_integral = %lf\n",nThreads,t_id,local_integral);
    }
    return integral;
}

double midpointIntegralOMP(double f(double x), double xMin, double xMax, long long int nPoints){
    double dx = (xMax - xMin)/nPoints;
    // printf("dx = %E...\n",dx);
    double integral = 0;
    // omp_set_num_threads(N_THREADS);
    #pragma omp parallel for reduction(+:integral) schedule(static)
        for(long long int i = 0; i<nPoints; i++){
            integral += f(xMin + i*dx + 0.5*dx)*dx;
        }
    return integral;
}

double accuracy(double numericalValue, double expectedValue){
    double err = fabs(numericalValue - expectedValue);
    // double rel_err = err/expectedValue;
    return err;
}

int main(){
    FILE *file,*timeFile;
    const int average = 200;
    file = fopen("data.txt","w");
    timeFile = fopen("time_data.txt","w");
    long long int number_of_points = pow(10,6);
    int i = 4;
    for(int i = 1; i < 66; i += 1){
        double total = 0;
        double sigma = 0;
        double totalOMP = 0;
        double sigmaOMP = 0;
        double numValue = midpointIntegral(f,0,1,number_of_points,i);
        double expValue = 1./3;
        // printf("integral = %.20lf\t%.20lf\n",numValue,expValue);
        double err_magnitude = accuracy(numValue,expValue);
        printf("(%lf,%lf)\terr = (%lf)\n",numValue,expValue,err_magnitude);
        fprintf(file,"%d\t%.32lf\n",i,err_magnitude);
        for(int j = 0; j<average; j++){
            double clock = omp_get_wtime();
            numValue = midpointIntegral(f,0,1,number_of_points,i);
            clock = omp_get_wtime() - clock;
            total += clock;
            sigma += clock*clock;
        }
        total /= average;
        sigma = sqrt((sigma/average - total*total)/average);
        // for(int j = 0; j<average; j++){
        //     double clock = omp_get_wtime();
        //     midpointIntegralOMP(f,0,1,number_of_points);
        //     clock = omp_get_wtime() - clock;
        //     totalOMP += clock;
        //     sigmaOMP += clock*clock;
        // }
        // totalOMP /= average;
        // sigmaOMP = sqrt((sigmaOMP/average - totalOMP*totalOMP)/average);
        fprintf(timeFile,"%d\t%lf\t%lf\t%lf\t%lf\n",i,total,sigma,totalOMP,sigmaOMP);
    }
    fclose(file);fclose(timeFile);
    printf("\n\n\n");
    return 0;
}