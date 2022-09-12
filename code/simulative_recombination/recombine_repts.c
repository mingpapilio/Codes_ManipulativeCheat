/**********************************************************************************************************
 * The is the source code for the individual-based simulation in 'The evolution of manipulative cheating'.
 * In this file, we simulate population with coevolving traits and log the trait values each 5 generations
 * Notes:
 *      1. Please make sure random number generator is located at the right directory (dSFMT)
 *      2. Switches (s1, s2, s3, s4) are implemented to make changing simulation settings easier. In particular
 *         , s1-3 is to make traits able to mutate, s4 is controlling the number of mechanisms to manipulate 
 *         and suppress manipulation.
 *      3. Change rec_rate for modifying the rate of recombination, where two haploids exchange one random trait
 *         value with each other with probability of rec_rate. When rec-rate is 0, it is equivalent to the oter 
 *         simulation code file.
 *      4. Feel free to email 'ming.liu.ac[at]gmail.com' if encounter problems with execution
 * Last edit: 12 Sep, 2022 (Ming Liu)
 **********************************************************************************************************/
// Including crucial packages and MersenneTwister random number generator
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
// Global functions
    double Mean_array(double p[], int length);
    double SD_array(double p[], int length);
    double sum(double p[], int length);
    int RandFromProb(double p[], int length, double RandNum);
    double normal_dist_BM (double mean, double sd, double u1, double u2);
    double **d_array2d(long size_1, long size_2);
    void free_d_array2d(double **x);
    double ****d_array4d(long size_1, long size_2, long size_3, long size_4);
    void free_d_array4d(double ****x);
    double safe_add(double a, double b);
// Main code
int main(){
    // Switches
        int s1= 1;                      // Whether selfishness is mutable in the model
        int s2= 1;                      // Whether manipulation is mutable in the model
        int s3= 1;                      // Whether suppression of manipulation is mutable in the model
        int s4= 1;                      // 0: single-mechanism model, 1: two-mechanisms model
    // General parameters
        int T= 3E3;                     // Number of generations 
        int N_rep= 10;                  // Number of repetitions
        int num_ini_idvl= 2;            // Number of initial group members (Relatedness= its inverse)
        int n= 100;                     // Number of individual in each social group
        int N_samp= 1E3;                // Number of sampling in calculating R
        int N_total= 1E5;               // Total number of individual in the population
        int G= ceil((double)N_total/n); // Number of groups
        double mut_rate= 0.01;          // Mutation rate
        double rec_rate= 0.0;           // Recombination rate
        double MutStep= 0.1;            // Size of mutation
        double a= 0.7;                  // Nonlinearity parameter fo x
        double b= 0.7;                  // Nonlinearity parameter fo u
        int N_trait= 5;

    // Output settings
        FILE *out_z, *out_x1, *out_u1, *out_x2, *out_u2;
        out_z= fopen("zlog.txt", "w");
        fprintf(out_z, "T\trep_1\trep_2\trep_3\trep_4\trep_5\trep_6\trep_7\trep_8\trep_9\trep_10\n");
        out_x1= fopen("x1log.txt", "w");
        fprintf(out_x1, "T\trep_1\trep_2\trep_3\trep_4\trep_5\trep_6\trep_7\trep_8\trep_9\trep_10\n");
        out_u1= fopen("u1log.txt", "w");
        fprintf(out_u1, "T\trep_1\trep_2\trep_3\trep_4\trep_5\trep_6\trep_7\trep_8\trep_9\trep_10\n");
        out_x2= fopen("x2log.txt", "w");
        fprintf(out_x2, "T\trep_1\trep_2\trep_3\trep_4\trep_5\trep_6\trep_7\trep_8\trep_9\trep_10\n");
        out_u2= fopen("u2log.txt", "w");
        fprintf(out_u2, "T\trep_1\trep_2\trep_3\trep_4\trep_5\trep_6\trep_7\trep_8\trep_9\trep_10\n");
        time_t time_start, time_end;
        time_start= time(NULL);

    // Random number genertor
        int seed;
        dsfmt_t dsfmt;
        seed= time(NULL);
        if(seed==0)seed= 1;
        dsfmt_init_gen_rand(&dsfmt,seed);
    // Variables
        int i,j,k,l;
        int t, rep, idx, rand_samp, g_samp, n_samp, num_div_idvl, new_geno, geno_counter, idvl_counter, compare_counter, rand_samp2, g_samp2, n_samp2;
        double pp, pp2, pp3, u1, u2, social_temp, rho, z_temp, x1_temp, u1_temp, x2_temp, u2_temp, cost, idvl_sum, nebr_sum, p_global, exgroup_freq, pop_sum, grup_sum;
        double avg_P, asd_P, avg_porp, asd_porp, sd_avg_P, sd_asd_P, sd_avg_porp, sd_asd_porp, mut_size;
        double ftemp1, z_grup, z_othr, x1_grup, x1_othr, u1_grup, u1_othr, x1_nsum, u1_nsum, x2_grup, x2_othr, u2_grup, u2_othr, x2_nsum, u2_nsum;

    // Temporary spaces
        double ****Popn= d_array4d(N_rep, G, n, N_trait);     // The population matrix
        double ****PopnNext= d_array4d(N_rep, G, n, N_trait); // The population in next generation, sampled from the offsprings
        double ****Fitness= d_array4d(N_rep, G, n, 1);
        double **GrpFit= d_array2d(N_rep, G);
        double **z_tmp= d_array2d(N_rep, G*n);
        double **x1_tmp= d_array2d(N_rep, G*n);
        double **u1_tmp= d_array2d(N_rep, G*n);
        double **x2_tmp= d_array2d(N_rep, G*n);
        double **u2_tmp= d_array2d(N_rep, G*n);
        //
        double globalZ[N_rep];
        double globalX1[N_rep];
        double globalU1[N_rep];
        double globalX2[N_rep];
        double globalU2[N_rep];
        double globalW[N_rep];
        double log_ini_idvl[num_ini_idvl][N_trait]; // p
        double Founder[G*num_ini_idvl][N_trait];
        double GrpTemp[n];

    // Index of the individual
        int id_z=   0; // Probability of becoming division 1
        int id_x1=  1;
        int id_u1=  2;
        int id_x2=  3;
        int id_u2=  4;

        int id_geno=0; // Genotype
        int id_freq=1; // Frequency
    // Initialization
        for(rep=0; rep<N_rep; ++rep){
            for(i=0; i<G; ++i){
                for(j=0; j<n; ++j) {
                    Popn[rep][i][j][id_z]= 0.0;
                    Popn[rep][i][j][id_x1]= 0.0;
                    Popn[rep][i][j][id_u1]= 0.0;
                    Popn[rep][i][j][id_x2]= 0.0;
                    Popn[rep][i][j][id_u2]= 0.0;
        }}}
    // Start of simulation
        for(t=1; t<= T; ++t){
            for(rep=0; rep<N_rep; ++rep){
                // Determine the fitness
                for(i=0; i<G; ++i){
                    // Population growth? Fixed at N, growth happens at the end of each generation (after drawing individuals)
                    z_grup= x1_grup= u1_grup= x1_nsum= u1_nsum= x2_grup= u2_grup= x2_nsum= u2_nsum= 0.0;
                    for(j=0; j<n; ++j){
                        z_grup+= Popn[rep][i][j][id_z]/n;
                        x1_grup+= Popn[rep][i][j][id_x1]/n;
                        u1_grup+= Popn[rep][i][j][id_u1]/n;
                        x1_nsum+= pow(Popn[rep][i][j][id_x1],a)/n;
                        u1_nsum+= pow(Popn[rep][i][j][id_u1],b)/n;
                        x2_grup+= Popn[rep][i][j][id_x2]/n;
                        u2_grup+= Popn[rep][i][j][id_u2]/n;
                        x2_nsum+= pow(Popn[rep][i][j][id_x2],a)/n;
                        u2_nsum+= pow(Popn[rep][i][j][id_u2],b)/n;
                    }
                    // NAN handling
                    if(z_grup!= z_grup) z_grup=0.0;
                    if(x1_grup!= x1_grup) x1_grup=0.0;
                    if(u1_grup!= u1_grup) u1_grup=0.0;
                    if(x2_grup!= x2_grup) x2_grup=0.0;
                    if(u2_grup!= u2_grup) u2_grup=0.0;
                    // Getting group fitness
                    GrpFit[rep][i]= 1-z_grup-x1_grup-u1_grup-x2_grup-u2_grup;
                    if(GrpFit[rep][i]<0.0) GrpFit[rep][i]= 0.0;
                    //
                    for(j=0; j<n; ++j){
                        //
                        z_othr= (z_grup*n- Popn[rep][i][j][id_z])/(n-1);
                        x1_othr= (x1_nsum*n- pow(Popn[rep][i][j][id_x1],a))/(n-1);
                        u1_othr= (u1_nsum*n- pow(Popn[rep][i][j][id_u1],b))/(n-1);
                        x2_othr= (x2_nsum*n- pow(Popn[rep][i][j][id_x2],a))/(n-1);
                        u2_othr= (u2_nsum*n- pow(Popn[rep][i][j][id_u2],b))/(n-1);
                        // Getting individual fitness
                        if(z_grup!= 0.0){
                            ftemp1= Popn[rep][i][j][id_z];
                            // ftemp1
                            if(Popn[rep][i][j][id_x1]!= 0.0){
                                ftemp1+= (1-u1_othr)*z_othr*pow(Popn[rep][i][j][id_x1],a);
                            }
                            if(x1_othr!= 0.0){
                                ftemp1+= -1*(1-pow(Popn[rep][i][j][id_u1],b))*Popn[rep][i][j][id_z]*x1_othr;
                            }
                            if(Popn[rep][i][j][id_x2]!= 0.0){
                                ftemp1+= (1-u2_othr)*z_othr*pow(Popn[rep][i][j][id_x2],a);
                            }
                            if(x2_othr!= 0.0){
                                ftemp1+= -1*(1-pow(Popn[rep][i][j][id_u2],b))*Popn[rep][i][j][id_z]*x2_othr;
                            }
                            //
                            if (ftemp1< 0.0) ftemp1=0.0;
                        } 
                        else ftemp1= 0.0;
                        //
                        Fitness[rep][i][j][0]= ftemp1;
                    }
                }
                // Mutation& reproduction
                for(i=0; i<G; ++i){
                    for(j=0; j< num_ini_idvl; ++j){
                    // Sample an individual to form a group
                        // Fitness-based selection
                        g_samp= RandFromProb(GrpFit[rep], G, dsfmt_genrand_open_close(&dsfmt));
                        for(k=0; k<n; ++k) GrpTemp[k]= Fitness[rep][g_samp][k][0];
                        n_samp= RandFromProb(GrpTemp, n, dsfmt_genrand_open_close(&dsfmt));

                        // Mutation of z, selfishness
                        if(s1==1 && t< T){
                            if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                            }
                            else mut_size= 0.0;
                            z_temp= Popn[rep][g_samp][n_samp][id_z]+ mut_size;
                        }
                        else z_temp= Popn[rep][g_samp][n_samp][id_z];
                        // Mutation of x, manipulation
                        if(s2==1 && t> 100){
                            if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                            }
                            else mut_size= 0.0;
                            x1_temp= Popn[rep][g_samp][n_samp][id_x1]+ mut_size;
                            //
                            if(s4==1){
                                if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                    mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                                }
                                else mut_size= 0.0;
                                x2_temp= Popn[rep][g_samp][n_samp][id_x2]+ mut_size;
                            }
                        }
                        else{
                            x1_temp= Popn[rep][g_samp][n_samp][id_x1];
                            x2_temp= Popn[rep][g_samp][n_samp][id_x2];
                        }

                        // Mutation of u, suppression
                        if(s3==1 && t> 200){
                            if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                            }
                            else mut_size= 0.0;
                            u1_temp= Popn[rep][g_samp][n_samp][id_u1]+ mut_size;
                            //
                            if(s4==1){
                                if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                    mut_size= MutStep*normal_dist_BM(0,1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                                }
                                else mut_size= 0.0;
                                u2_temp= Popn[rep][g_samp][n_samp][id_u2]+ mut_size;
                            }
                        }
                        else{
                            u1_temp= Popn[rep][g_samp][n_samp][id_u1];
                            u2_temp= Popn[rep][g_samp][n_samp][id_u2];
                        }
                    // Dealing with out-of-boundary values
                        if(z_temp> 1.0) z_temp= 1.0-1E-6;
                        if(z_temp< 0.0) z_temp= 1E-6;
                        if(x1_temp> 1.0) x1_temp= 1.0-1E-6;
                        if(x1_temp< 0.0) x1_temp= 1E-6;
                        if(u1_temp> 1.0) u1_temp= 1.0-1E-6;
                        if(u1_temp< 0.0) u1_temp= 1E-6;
                        if(x2_temp> 1.0) x2_temp= 1.0-1E-6;
                        if(x2_temp< 0.0) x2_temp= 1E-6;
                        if(u2_temp> 1.0) u2_temp= 1.0-1E-6;
                        if(u2_temp< 0.0) u2_temp= 1E-6;
                    // Logging trait values
                        log_ini_idvl[j][id_z]= z_temp;
                        log_ini_idvl[j][id_x1]= x1_temp;
                        log_ini_idvl[j][id_u1]= u1_temp;
                        log_ini_idvl[j][id_x2]= x2_temp;
                        log_ini_idvl[j][id_u2]= u2_temp;
                    }
                    // Recombination
                    for(j=0; j<(num_ini_idvl/2); ++j){
                        if(dsfmt_genrand_open_close(&dsfmt)< rec_rate){
                            if(s4==0){
                                rand_samp= floor(dsfmt_genrand_open_close(&dsfmt)*3);
                                z_temp= log_ini_idvl[j][rand_samp];
                                log_ini_idvl[j][rand_samp]= log_ini_idvl[j+1][rand_samp];
                                log_ini_idvl[j+1][rand_samp]= z_temp;
                            }
                            if(s4==1){
                                rand_samp= floor(dsfmt_genrand_open_close(&dsfmt)*5);
                                z_temp= log_ini_idvl[j][rand_samp];
                                log_ini_idvl[j][rand_samp]= log_ini_idvl[j+1][rand_samp];
                                log_ini_idvl[j+1][rand_samp]= z_temp;
                            }
                        }
                    }
                    // Put all founders to the founder pool
                    for(j=0; j< num_ini_idvl; ++j){
                        idx= i*num_ini_idvl+j;
                        for(k=0; k<N_trait; ++k) Founder[idx][k]= log_ini_idvl[j][k];
                    }
                }
                // Forming the new population by sampling founders from founder pool
                for(i=0; i<G; ++i){
                    //
                    for(j=0; j< num_ini_idvl; ++j){
                        rand_samp= floor(dsfmt_genrand_open_close(&dsfmt)*num_ini_idvl*G);
                        for(k=0; k<N_trait; ++k) log_ini_idvl[j][k]= Founder[rand_samp][k];
                    }
                    //
                    for(j=0; j<n; ++j){
                        rand_samp= floor(dsfmt_genrand_open_close(&dsfmt)*num_ini_idvl);
                        // 
                        PopnNext[rep][i][j][id_z]= log_ini_idvl[rand_samp][id_z];
                        PopnNext[rep][i][j][id_x1]= log_ini_idvl[rand_samp][id_x1];
                        PopnNext[rep][i][j][id_u1]= log_ini_idvl[rand_samp][id_u1];
                        PopnNext[rep][i][j][id_x2]= log_ini_idvl[rand_samp][id_x2];
                        PopnNext[rep][i][j][id_u2]= log_ini_idvl[rand_samp][id_u2];
                    }
                // End of group
                }
            // End of repetition
            }
            
            if(t%5==0) {
                fprintf(out_z, "%d\t", t);
                fprintf(out_x1, "%d\t", t);
                fprintf(out_u1, "%d\t", t);
                fprintf(out_x2, "%d\t", t);
                fprintf(out_u2, "%d\t", t);
                //
                for(rep=0; rep<N_rep; ++rep){
                    // Calculate average P
                    for(i=0; i<G; ++i){
                        for(j=0; j<n; ++j) {
                            z_tmp[rep][i*n+j]= Popn[rep][i][j][id_z];
                            x1_tmp[rep][i*n+j]= Popn[rep][i][j][id_x1];
                            u1_tmp[rep][i*n+j]= Popn[rep][i][j][id_u1];
                            x2_tmp[rep][i*n+j]= Popn[rep][i][j][id_x2];
                            u2_tmp[rep][i*n+j]= Popn[rep][i][j][id_u2];
                        }
                    }
                    globalZ[rep]= Mean_array(z_tmp[rep], G*n);
                    globalX1[rep]= Mean_array(x1_tmp[rep], G*n);
                    globalU1[rep]= Mean_array(u1_tmp[rep], G*n);
                    globalX2[rep]= Mean_array(x2_tmp[rep], G*n);
                    globalU2[rep]= Mean_array(u2_tmp[rep], G*n);
                    // 
                    fprintf(out_z, "%lf\t", globalZ[rep]);
                    fprintf(out_x1, "%lf\t", globalX1[rep]);
                    fprintf(out_u1, "%lf\t", globalU1[rep]);
                    fprintf(out_x2, "%lf\t", globalX2[rep]);
                    fprintf(out_u2, "%lf\t", globalU2[rep]);
                }
                fprintf(out_z, "\n");
                fprintf(out_x1, "\n");
                fprintf(out_u1, "\n");
                fprintf(out_x2, "\n");
                fprintf(out_u2, "\n");
            }
            // Producing the next generation
            for(rep=0; rep<N_rep; ++rep){
                for(i=0; i<G; ++i){
                    for(j=0; j<n; ++j){
                        for(k=0; k<N_trait; ++k) Popn[rep][i][j][k]= PopnNext[rep][i][j][k];
            }}}
        }
    // Erasing all the matrices
        time_end= time(NULL);
        printf("This simulation lasted %ld seconds.\n",time_end-time_start);
        free_d_array2d(z_tmp);
        free_d_array2d(x1_tmp);
        free_d_array2d(u1_tmp);
        free_d_array2d(x2_tmp);
        free_d_array2d(u2_tmp);
        free_d_array4d(Fitness);
        free_d_array2d(GrpFit);
        free_d_array4d(Popn);
        free_d_array4d(PopnNext);
        fclose(out_z);
        fclose(out_x1);
        fclose(out_u1);
        fclose(out_x2);
        fclose(out_u2);
    // End
    return EXIT_SUCCESS;
}
////////////////////// Functions are defined in below //////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp;
    temp= 0.0;
    if(length> 0){
        for (i=0; i<length; ++i) temp+= p[i]/length;
        return temp;
    }
    else return NAN;
}

double SD_array(double p[], int length){
    if(length> 1){
        int i;
        double avg, temp;
        avg= temp= 0.0;
        for (i=0; i<length; ++i) avg+= p[i]/length;
        for (i=0; i<length; ++i) temp+= pow(p[i]-avg, 2);
        return sqrt(temp/(length-1));
    }
    else return NAN;
}

double sum(double p[], int length){
    int i;
    double sum=0.0;
    for(i=0; i< length; ++i) sum= safe_add(sum, p[i]);
    return sum;
}

int RandFromProb(double p[], int length, double RandNum){
    int i,temp,check;
    double sum, pp;
    sum= pp= 0.0;
    check= temp= 0;
    for(i=0; i< length; ++i){
        if(p[i]> 0.0) sum= safe_add(sum, p[i]/length);
    }
    if(sum> 0.0){
        for(i=0; i< length; ++i){
            if(RandNum> pp/sum) check=1;
            if(p[i]> 0.0) pp+= p[i]/length;
            if(RandNum< pp/sum&& check==1){
                temp= i;
                break;
    }}}
    else{
        temp= floor(length*RandNum);
    }
    return temp;
}

double normal_dist_BM (double mean, double sd, double u1, double u2){
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}

double **d_array2d(long size_1, long size_2){
    double **x;
    long i;
    x= (double **) malloc((size_t)(size_1*sizeof(double)));
    for(i=0; i< size_1; i++) x[i]= (double *) malloc((size_t)(size_2*sizeof(double)));
    return x;
}
void free_d_array2d(double **x){
    free(x[0]);
    free(x);
}

double ****d_array4d(long size_1, long size_2, long size_3, long size_4){
    double ****x;
    long i,j,k;
    x= (double ****) malloc((size_t)(size_1*sizeof(double ***)));
    for(i=0; i< size_1; i++) {
        x[i]= (double ***) malloc((size_t)(size_2*sizeof(double **)));
        for(j=0; j< size_2; j++){
            x[i][j]= (double **) malloc((size_t)(size_3*sizeof(double *)));
            for(k=0; k< size_3; k++) x[i][j][k]= (double *) malloc((size_t)(size_4*sizeof(double)));
        }
    }
    return x;
}

void free_d_array4d(double ****x){
    free(x[0][0][0]);
    free(x[0][0]);
    free(x[0]);
    free(x);
}

double safe_add(double a, double b){
    if (a> 0 && b> DBL_MAX- a) printf("Overflow in safe_add.\n");
    else if (a< 0 && b< - DBL_MAX- a) printf("Underflow in safe_add.\n");
    // DBL_MIN is the smallest (closest value to 0) floating number
    return a+b;
}

