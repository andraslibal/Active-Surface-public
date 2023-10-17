
//Active Mattter generic code to be used for future active matter projects
//one file, modular functions and adaptable to different new projects no hard coding
//nothing fancy, only considerations are simplicity, understandability and speed of execution
//Copyright @Libal Andras 2020 May - Simulations and Modeling Group

//fix color fix radius
//transmit theta place a small arrow


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "random_numrec.h"

#define PI 3.14159265358979323846

// pinningsite distribution type
// -1 - uniform distribution
// 0 - density growing with y (/ distribution)
// 1 - ^ distribution
// 2 - same as 1, but doubled
// 3 - same as 2, but doubled
// 4 - placing pinningsites in ^ shapes
#define dist_type -1

// #define non_periodic_boundary

// uncomment this #define to include drive force
// 0,3 - perpendicular
// 1,4 - gradient
// 2,5 - opposite
// #define drive_forces_type 4

struct particle_struct
{
    // index of particle
    int i;
    int clusterID;
    struct particle_struct* next;
};

struct cluster_struct
{
    // number of particles in cluster
    int N;
    struct particle_struct* first;
    struct particle_struct* last;
};

struct global_struct
    {
        char baseName[200];
        uint64_t seedToSet;

        double SX,SY;
        double SX2,SY2;
        int N_particles;
        int N_pinningsites;    

        //disc shaped particle
        double *particle_x;
        double *particle_y;
        double *particle_R;
        int *particle_color;
        //forces on the particle
        double *particle_fx;
        double *particle_fy;
        //distance since last Verlet build
        double *particle_dx_so_far;
        double *particle_dy_so_far;
        //active particle
        double *particle_motor_force;
        double *particle_angle_rad;
        double *particle_cosfi;
        double *particle_sinfi;
        int *particle_motor_ellapsed_time;
        int *particle_motor_total_time;
        //passive particle
        int *particle_is_passive;

        //for clustering the particles
        struct particle_struct* particles;
        struct cluster_struct* particle_cluster;

        double *pinningsite_x;
        double *pinningsite_y;
        double *pinningsite_R;
        double *pinningsite_fmax;

        double drive_force;
        
        //grid for the pinning sites
        //2D specifying the grid x,y position
        //thiord dimension: 4 because max 4 pinning sites fit
        int ***pinningsite_grid;
        double pinning_grid_cell_size;
        int pinningsite_grid_Nx,pinningsite_grid_Ny;

        // grid for active sites
        double active_grid_division;
        double active_grid_division_half;
        int active_grid_neighbours;
        double **active_grid;
        int active_grid_Nx, active_grid_Ny;

        // mask matrices
        int **matrix_xp_yp;
        int **matrix_xp_ym;
        int **matrix_xm_yp;
        int **matrix_xm_ym;

        // patern for absorption part
        int *pattern_x;
        int *pattern_y;

        // active grid changes
        double active_grid_absorption;
        double active_grid_recovery;

        // stat
        int nr_non_passive_particles;
        int nr_moving_partices;
        double largest_x_com, largest_y_com;
        double largest_average_dr;
        double average_x_moving_in_passive_particle;
        double average_y_moving_in_passive_particle;
        double average_x_moving_in_moving_particle;
        double average_y_moving_in_moving_particle;
        double average_x_moving_in_standing_particle;
        double average_y_moving_in_standing_particle;

        double passive_particle_percentage;

        double threshold_l;
        double average_neighbours_inside_threshold;
        int **neighbours_inside_threshold;
        int *neighbours_inside_threshold_N;
        int **neighbours_inside_threshold_tmp;
        int *neighbours_inside_threshold_tmp_N;

        //generic particle properties
        double generic_particle_R;
        double generic_particle_k_spring;
        int generic_particle_motor_minimum_time;
        int generic_particle_motor_maximum_time;
        //generic pinningsite properties
        double generic_pinningsite_R;
        double generic_pinningsite_fxmax;

        //simulation time step
        double dt;

        int time;
        int total_runtime;

        int echo_time;
        int movie_time;
        int stat_time;

        FILE *movie_file;
        FILE *stat_file;

        int *Verlet_i;
        int *Verlet_j;
        int N_Verlet;
        char flag_Verlet_needs_rebuild;
        double Verlet_distance_to_rebuild;
        double Verlet_R;

    } global;


void initialize_from_file(int argc,char *argv[]);
void read_line_from_file(FILE *parameterfile,char descript[],char type[],void *value);
void initialize_particles();
void initialize_pinningsites();
void initialize_pinningsites_not_randomly();
void initialize_pinning_grid();

void initialize_matrices();

void initialize_active_grid();

void calculate_PBC_folded_distance(double *dx, double *dy,double x1,double y1,double x2,double y2);
void write_particles();
void write_pinningsites();

void calculate_external_forces();

void calculate_active_grid_forces();
void calculate_active_grid_recovery();

void calculate_internal_motor_forces();
void calculate_interparticle_forces();
void calculate_pinning_forces();
void calculate_interaction_with_wall();
void move_particles();
void clusterize_particles();
int largest_clusterID();
int number_of_clusters();
void centroid_of_cluster(double* x_com, double* y_com, int clusterID);
void calculate_average_distance_from_centroid_of_cluster();

void check_Verlet_need_to_rebuild(int i);
void rebuild_Verlet();

void check_affine_displacement();

void run_simulation();
void write_cmovie_frame();
void write_stat_data();

void density_distribution_to_file(int time_step, int N_bands);


void initialize_from_file(int argc,char *argv[])
{
FILE *parameterfile;
char mvi_fname[300];
char stat_fname[300];

    if (argc<2) 
        {
            printf("Usage: ./active <parameter_file>\n");
            fflush(stdout);
            exit(1);
        }

    parameterfile = fopen(argv[1],"rt");
    if (parameterfile==NULL)
        {
            printf("Cannot open parameter file %s\n",argv[1]);
            fflush(stdout);
            exit(1);
        }
    printf("Opened parameter file %s\n",argv[1]);
    read_line_from_file(parameterfile,"SX","double",&global.SX);
    read_line_from_file(parameterfile,"SY","double",&global.SY);
    global.SX2 = global.SX / 2.0;
    global.SY2 = global.SY / 2.0;
    read_line_from_file(parameterfile,"N_particles","int",&global.N_particles);
    read_line_from_file(parameterfile,"N_pinningsites","int",&global.N_pinningsites);   
    read_line_from_file(parameterfile,"dt","double",&global.dt);
    read_line_from_file(parameterfile,"total_runtime","int",&global.total_runtime);
    read_line_from_file(parameterfile,"echo_time","int",&global.echo_time);
    read_line_from_file(parameterfile,"movie_time","int",&global.movie_time);
    read_line_from_file(parameterfile,"stat_time","int",&global.stat_time);

    read_line_from_file(parameterfile,"generic_particle_R","double",&global.generic_particle_R);
    read_line_from_file(parameterfile,"generic_particle_k_spring","double",&global.generic_particle_k_spring);
    read_line_from_file(parameterfile,"generic_particle_motor_minimum_time","int",&global.generic_particle_motor_minimum_time);
    //read_line_from_file(parameterfile,"generic_particle_motor_maximum_time","int",&global.generic_particle_motor_maximum_time);
    global.generic_particle_motor_maximum_time = 2 * global.generic_particle_motor_minimum_time;
    printf("Calculated: generic_particle_motor_maximum_time = %d\n",global.generic_particle_motor_maximum_time);
    read_line_from_file(parameterfile,"generic_pinningsite_R","double",&global.generic_pinningsite_R);
    read_line_from_file(parameterfile,"generic_pinningsite_fxmax","double",&global.generic_pinningsite_fxmax);

    read_line_from_file(parameterfile,"active_grid_division","double",&global.active_grid_division);
    read_line_from_file(parameterfile,"active_grid_absorption","double",&global.active_grid_absorption);
    read_line_from_file(parameterfile,"active_grid_recovery","double",&global.active_grid_recovery);

    read_line_from_file(parameterfile,"seedToSet","uint64_t",&global.seedToSet);

    read_line_from_file(parameterfile,"passive_particle_percentage","double",&global.passive_particle_percentage);

    if (global.seedToSet != 0) setseed(global.seedToSet);

    //Verlet list pre-init
    global.Verlet_distance_to_rebuild = global.generic_particle_R;
    printf("Calculated: Verlet_distance_to_rebuild = %lf\n",global.Verlet_distance_to_rebuild);
    global.Verlet_i = NULL;
    global.Verlet_j = NULL;
    global.N_Verlet = 0;
    global.Verlet_R = 4.0 * global.generic_particle_R;
    printf("Calculated: Verlet_R = %lf\n",global.Verlet_R);
    printf("Particle density is = %lf\n",global.N_particles * global.generic_particle_R * global.generic_particle_R * PI / global.SX / global.SY);

    //open files
    #ifdef drive_forces_type
        global.drive_force = 1.;

        #if drive_forces_type == 0 || drive_forces_type == 3
            if (global.seedToSet == 0)
                sprintf(global.baseName, "perpendicular/%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
            else
                sprintf(global.baseName, "seeds/perpendicular/%ld_%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", global.seedToSet, dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
        #elif drive_forces_type == 1 || drive_forces_type == 4
            if (global.seedToSet == 0)
                sprintf(global.baseName, "gradient/%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
            else
                sprintf(global.baseName, "seeds/gradient/%ld_%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", global.seedToSet, dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
        #elif drive_forces_type == 2 || drive_forces_type == 5
            if (global.seedToSet == 0)
                sprintf(global.baseName, "opposite/%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
            else
                sprintf(global.baseName, "seeds/opposite/%ld_%d%dN%d_Npin%d_fpin%.2lf_fD%.2lf_minRuntime%d_Tsteps%d_dt%lf", global.seedToSet, dist_type, drive_forces_type, global.N_particles, global.N_pinningsites, global.generic_pinningsite_fxmax, global.drive_force, global.generic_particle_motor_minimum_time, global.total_runtime, global.dt);
        #endif
    #else
        if (global.seedToSet == 0)
            sprintf(global.baseName, "s%dN%d_p_perc%.1lf_Tsteps%d_dt%lf_actDiv%lf_actAbs%lf_actRec%lf", dist_type, global.N_particles, global.passive_particle_percentage, global.total_runtime, global.dt, global.active_grid_division, global.active_grid_absorption, global.active_grid_recovery);
        else
            sprintf(global.baseName, "seeds/seed_%ld_s%dN%d_p_perc%.1lf_Tsteps%d_dt%lf_actDiv%lf_actAbs%lf_actRec%lf", global.seedToSet, dist_type, global.N_particles, global.passive_particle_percentage, global.total_runtime, global.dt, global.active_grid_division, global.active_grid_absorption, global.active_grid_recovery);
    #endif

    printf("Base name: %s\n", global.baseName);

    sprintf(mvi_fname, "movies/%s.mvi", global.baseName);
    global.movie_file = fopen(mvi_fname,"wb");

    sprintf(stat_fname, "stat/%s.txt", global.baseName);
    global.stat_file = fopen(stat_fname,"w");
}

void read_line_from_file(FILE *parameterfile,char descript[],char type[],void *value)
{
    double temp_double;
    int temp_int;
    uint64_t temp_uint64_t;
    char temp_descript[100];
    int read_items;
    
    if (strcmp(type,"double")==0)
        read_items = fscanf(parameterfile,"%s %lf",temp_descript,&temp_double);
    else if (strcmp(type,"int")==0)
        read_items = fscanf(parameterfile,"%s %d",temp_descript,&temp_int);
    else if (strcmp(type,"uint64_t")==0)
        read_items = fscanf(parameterfile,"%s %ld",temp_descript,&temp_uint64_t);

    if (read_items!=2)
        {
            printf("Error on parameter file line %s\n",descript);
            exit(1);
        }

    if (strcmp(temp_descript,descript)==0) 
        {

        if (strcmp(type,"double")==0)
            {
                *(double *)value = temp_double;
                printf("Parameter : %s = %lf\n",descript,*(double *)value);
            }
        else if (strcmp(type,"int")==0)
            {
                *(int *)value = temp_int;
                printf("Parameter : %s = %d\n",descript,*(int *)value);
            }
        else if (strcmp(type,"uint64_t")==0)
            {
                *(uint64_t *)value = temp_uint64_t;
                printf("Parameter : %s = %ld\n",descript,*(uint64_t *)value);
            }


        }
    else 
        {
            printf("Expected %s in parameterfile\n",descript);
            exit(1);
        }
}

void initialize_particles()
{
    int i,j,N;
    double min_dist;
    double density;
    double x_temp,y_temp;
    double dx,dy,dr;
    char overlap;

    //this means the discs cover about half of the system (0.4)
    //with random deposition we can go up to around 0.5
    density = 0.4;
    min_dist = sqrt(global.SX * global.SY / (double)global.N_particles / PI * density);
    printf("Min dist = %lf Max radius = %lf\n",min_dist,min_dist/2.0);

    global.particle_x = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_y = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_R = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_color = (int *)malloc(global.N_particles*sizeof(int));

    global.particle_fx = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_fy = (double *)malloc(global.N_particles*sizeof(double));

    global.particle_dx_so_far = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_dy_so_far = (double *)malloc(global.N_particles*sizeof(double));

    global.particle_motor_force = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_angle_rad = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_cosfi = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_sinfi = (double *)malloc(global.N_particles*sizeof(double));
    global.particle_motor_ellapsed_time = (int *)malloc(global.N_particles*sizeof(int));
    global.particle_motor_total_time = (int *)malloc(global.N_particles*sizeof(int));

    global.particle_cluster = (struct cluster_struct*)malloc(global.N_particles*sizeof(struct cluster_struct));
    global.particles = (struct particle_struct*)malloc(global.N_particles*sizeof(struct particle_struct));

    global.particle_is_passive = (int *)malloc(global.N_particles*sizeof(int));

    global.nr_non_passive_particles = 0;

    global.average_x_moving_in_passive_particle = 0.0;
    global.average_y_moving_in_passive_particle = 0.0;
    global.average_x_moving_in_moving_particle = 0.0;
    global.average_y_moving_in_moving_particle = 0.0;
    global.average_x_moving_in_standing_particle = 0.0;
    global.average_y_moving_in_standing_particle = 0.0;

    global.average_neighbours_inside_threshold = 0.0;

    global.neighbours_inside_threshold = (int **)malloc(global.N_particles*sizeof(int*));
    global.neighbours_inside_threshold_N = (int *)malloc(global.N_particles*sizeof(int));

    global.neighbours_inside_threshold_tmp = (int **)malloc(global.N_particles*sizeof(int*));
    global.neighbours_inside_threshold_tmp_N = (int *)malloc(global.N_particles*sizeof(int));

    for(i=0;i<global.N_particles;i++)
        {

            overlap = 1;

            //shoould have an exit condition for very dense systems
            //if the density was set right it will not lock into this
            while (overlap) 
                {

                x_temp = global.SX * Rand();
                y_temp = global.SY * Rand();
                overlap = 0;

                for(j=0;j<i;j++)
                    {
                        calculate_PBC_folded_distance(&dx,&dy,x_temp,y_temp,global.particle_x[j],global.particle_y[j]);
                        dr =  sqrt(dx*dx+dy*dy);
                        if (dr<min_dist)
                        {
                            overlap = 1;
                            break;
                        }
                        #ifdef non_periodic_boundary
                        else if (global.SX-x_temp<=global.generic_particle_R || x_temp<=global.generic_particle_R || global.SY-y_temp<=global.generic_particle_R || y_temp<=global.generic_particle_R)
                        {
                            overlap = 1;
                            break;
                        }
                        #endif
                    }
                } 

            //found a good position that does not overlap
            global.particle_x[i] = x_temp;
            global.particle_y[i] = y_temp;

            //radius of the particle
            global.particle_R[i] = global.generic_particle_R;

            //randomly half will go right or left this is for testing
            // if (Rand() < 0.5)
            //     {
            //         global.particle_color[i] = 2;
            //         global.particle_angle_rad[i] = 0.0;
            //         global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
            //         global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);
            //     }
            // else
            //     {
            //         global.particle_color[i] = 3;
            //         global.particle_angle_rad[i] = PI;
            //         global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
            //         global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);
            //     }
            

            global.particle_angle_rad[i] = 2.0 * PI * Rand();
            global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
            global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);
            global.particle_color[i] = 4;
            // if (global.particle_cosfi[i]<0) global.particle_color[i] = 2;
            // else global.particle_color[i] = 3;

            global.particle_motor_force[i] = 1.0;
            global.particle_motor_ellapsed_time[i] = 0;
            global.particle_motor_total_time[i] = global.generic_particle_motor_minimum_time;
            global.particle_motor_total_time[i] += (int)floor(Rand()*(global.generic_particle_motor_maximum_time - global.generic_particle_motor_minimum_time));

            global.particles[i].i = i;

            //printf("%d %d\n",i,global.particle_motor_total_time[i]);

            global.particle_fx[i] = 0.0;
            global.particle_fy[i] = 0.0;

            global.particle_dx_so_far[i] = 0.0;
            global.particle_dy_so_far[i] = 0.0;

            if (Rand() < global.passive_particle_percentage) {
                global.particle_is_passive[i] = 1;
            }
            else {
                global.particle_is_passive[i] = 0;
                global.nr_non_passive_particles += 1;
            }

            global.neighbours_inside_threshold_N[i] = 0;
            global.neighbours_inside_threshold[i] = (int *)malloc(global.neighbours_inside_threshold_N[i]*sizeof(int));

            global.neighbours_inside_threshold_tmp_N[i] = 0;
            global.neighbours_inside_threshold_tmp[i] = (int *)malloc(global.neighbours_inside_threshold_tmp_N[i]*sizeof(int));
            
        }

    printf("%d Particles initialized successfully\n",global.N_particles);
}

#ifdef dist_type
    #if dist_type == 0 || dist_type == 2
        // inverse transform sampling
        // current CDF (cumulative distribution function): f(x) = x^2 / a^2
        // because probability density function is g(x) = 2*x / a^2
        double arbitrary_distribution(double x)
        {
            return sqrt(x);
        }
    #elif dist_type == 1 || dist_type == 3
        // inverse transform sampling
        // current CDF (cumulative distribution function): f(x) = 2 / a^2 * x^2 (if x in [0, a/2]), and f(x) = 2 / a * x * (2 - x / a) - 1 (if x in (a/2, a])
        // because probability density function is g(x) = 4 / a^2 * x (if x in [0, a/2]), and g(x) = -4 / a^2 * (x - a) (if x in (a/2, a])
        double arbitrary_distribution(double x)
        {
            if (x <= 0.5) return sqrt(x / 2.);
            return (2. - sqrt(2. * (1. - x))) / 2.;
        }
    #elif dist_type == 4
        double y_V_shape(double x, double x0, double y0, double width, double height)
        {
            if (x >= x0 && x <= x0+width/2) return y0 + 2 * height / width * (x - x0);
            else if (x > x0+width/2 && x <= x0+width) return y0 - 2 * height / width * (x - x0 - width);
        }

        void place_in_V_shape(double x, double y, double width, double height, double dist_between_sites, int i)
        {
            double dx, dy, cos_alpha, x_temp;
            cos_alpha = cos(atan(2 * height / width));

            x_temp = x;
            while (x_temp < x+width && i < global.N_pinningsites)
            {
                global.pinningsite_x[i] = x_temp;
                global.pinningsite_y[i] = y_V_shape(x_temp, x, y, width, height);
                global.pinningsite_R[i] = global.generic_pinningsite_R;
                global.pinningsite_fmax[i] = global.generic_pinningsite_fxmax;

                i++;
                x_temp += (2 * global.generic_pinningsite_R + dist_between_sites) * cos_alpha;
            }
        }
    #endif
#endif

void initialize_pinningsites()
{
    int i,j;
    double min_dist;
    double density;
    double x_temp,y_temp;
    double dx,dy,dr;
    char overlap;
    
    min_dist = 2 * global.generic_pinningsite_R;
    //covering only half the area with pinning sites
    #if dist_type == -1
        density = (double)(global.N_pinningsites * PI * global.generic_pinningsite_R * global.generic_pinningsite_R) / (global.SX * global.SY *0.5);
    #elif dist_type >= 0 && dist_type <= 3
        density = (double)(global.N_pinningsites * PI * global.generic_pinningsite_R * global.generic_pinningsite_R) / (global.SX * global.SY);
    #endif
    printf("Pinning site density is = %lf\n",density);

    global.pinningsite_x = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_y = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_R = (double *)malloc(global.N_pinningsites*sizeof(double));
    global.pinningsite_fmax = (double *)malloc(global.N_pinningsites*sizeof(double));


    for(i=0;i<global.N_pinningsites;i++)
        {

            overlap = 1;

            //shoould have an exit condition for very dense systems
            //if the density was set right it will not lock into this
            while (overlap) 
                {

                x_temp = global.SX * Rand();
                #ifdef dist_type
                    #if dist_type == 0 || dist_type == 1
                        y_temp = global.SY * arbitrary_distribution(Rand());
                    #elif dist_type == 2 || dist_type == 3
                        if (Rand() < 0.5) y_temp = global.SY / 2. * arbitrary_distribution(Rand());
                        else y_temp = global.SY / 2. * arbitrary_distribution(Rand()) + global.SY / 2.;
                    #elif dist_type == -1
                        y_temp = global.SY * Rand() * 0.5;
                    #endif
                #endif
                overlap = 0;

                for(j=0;j<i;j++)
                    {
                        calculate_PBC_folded_distance(&dx,&dy,x_temp,y_temp,global.pinningsite_x[j],global.pinningsite_y[j]);
                        dr =  sqrt(dx*dx+dy*dy);
                        if (dr<min_dist)
                        {
                            overlap = 1;
                            break;
                        }
                        #ifdef non_periodic_boundary
                        else if (global.SX-x_temp<=global.generic_pinningsite_R || x_temp<=global.generic_pinningsite_R || global.SY-y_temp<=global.generic_pinningsite_R || y_temp<=global.generic_pinningsite_R)
                        {
                            overlap = 1;
                            break;
                        }
                        #endif
                    }
                } 

            //found a good position that does not overlap
            global.pinningsite_x[i] = x_temp;
            global.pinningsite_y[i] = y_temp;

            //radius of the pinningsite
            global.pinningsite_R[i] = global.generic_pinningsite_R;
            global.pinningsite_fmax[i] = global.generic_pinningsite_fxmax;
        }

    printf("%d Pinningsites initialized successfully\n",global.N_pinningsites);
}

void initialize_pinningsites_not_randomly()
{

}

void add_pinningsite_to_grid(int i, int j, int l)
{
    int k;

    global.pinningsite_grid[i][j][0] ++;
    k = global.pinningsite_grid[i][j][0];
    global.pinningsite_grid[i][j] = (int*) realloc(global.pinningsite_grid[i][j], (k+1) * sizeof(int));
    global.pinningsite_grid[i][j][k] = l;
}

void initialize_pinning_grid()
{
    int i, j, l;

    global.pinning_grid_cell_size = 2.0 * global.generic_pinningsite_R;
    global.pinningsite_grid_Nx = (int)floor(global.SX /global.pinning_grid_cell_size) + 1;
    global.pinningsite_grid_Ny = (int)floor(global.SY /global.pinning_grid_cell_size) + 1;
    
    printf("Pinning grid size = %lf\n",global.pinning_grid_cell_size);
    printf("Pinning grid size %d x %d\n",global.pinningsite_grid_Nx,global.pinningsite_grid_Ny);

    global.pinningsite_grid = (int***) malloc(global.pinningsite_grid_Nx * sizeof(int**));
    for (i=0; i<global.pinningsite_grid_Nx; i++)
    {
        global.pinningsite_grid[i] = (int**) malloc(global.pinningsite_grid_Ny * sizeof(int*));
        for (j=0; j<global.pinningsite_grid_Ny; j++)
        {
            //first number determines the number of pinningsites in a grid cell
            global.pinningsite_grid[i][j] = (int*) malloc(sizeof(int));
            global.pinningsite_grid[i][j][0] = 0;
        }
    }

    for (l=0; l<global.N_pinningsites; l++)
    {
        i = (int) floor(global.pinningsite_x[l] / global.pinning_grid_cell_size);
        j = (int) floor(global.pinningsite_y[l] / global.pinning_grid_cell_size);

        add_pinningsite_to_grid(i, j, l);
    }
}

void initialize_matrices() {
	
	int i, j, k;
	int act_matrix_index;

	global.active_grid_neighbours = pow((int)(2 * global.generic_particle_R / global.active_grid_division) + 1, 2) / 2;

	global.matrix_xm_yp = (int**) malloc(4 * sizeof(int*));
	global.matrix_xp_ym = (int**) malloc(4 * sizeof(int*));
	global.matrix_xp_yp = (int**) malloc(4 * sizeof(int*));
	global.matrix_xm_ym = (int**) malloc(4 * sizeof(int*));
    for (i=0; i<4; i++)
    {
        global.matrix_xm_yp[i] = (int*) malloc(global.active_grid_neighbours * sizeof(int));
        global.matrix_xp_ym[i] = (int*) malloc(global.active_grid_neighbours * sizeof(int));
        global.matrix_xp_yp[i] = (int*) malloc(global.active_grid_neighbours * sizeof(int));
        global.matrix_xm_ym[i] = (int*) malloc(global.active_grid_neighbours * sizeof(int));
        for (j=0; j<global.active_grid_neighbours; j++)
        {
            global.matrix_xm_yp[i][j] = 0;
            global.matrix_xp_ym[i][j] = 0;
            global.matrix_xp_yp[i][j] = 0;
            global.matrix_xm_ym[i][j] = 0;
        }
    }

    // xm_yp

    act_matrix_index = 0;

	for (i=-2; i <= 1; i++) {
		for (j=0; j <= 1; j++) {
			global.matrix_xm_yp[0][act_matrix_index] = i;
			global.matrix_xm_yp[1][act_matrix_index] = j;
			if (i < 0) {
				global.matrix_xm_yp[2][act_matrix_index] = -1;
				if (j <= 0) {
					global.matrix_xm_yp[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xm_yp[3][act_matrix_index] = 1;	
				}
			}
			else {
				global.matrix_xm_yp[2][act_matrix_index] = 1;
				if (j <= 0) {
					global.matrix_xm_yp[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xm_yp[3][act_matrix_index] = 1;	
				}
			}
			act_matrix_index++;
		}
	}

	for (i=-1; i <= 0; i++) {
		for (j=-1; j <= 2; j++) {
			for (k=0; k<act_matrix_index; k++) {
				if ((global.matrix_xm_yp[0][k] == i) && (global.matrix_xm_yp[1][k] == j)) {
					break;
				}
			}
			if (k == act_matrix_index) {
				global.matrix_xm_yp[0][act_matrix_index] = i;
				global.matrix_xm_yp[1][act_matrix_index] = j;
				if (i < 0) {
					global.matrix_xm_yp[2][act_matrix_index] = -1;
					if (j <= 0) {
						global.matrix_xm_yp[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xm_yp[3][act_matrix_index] = 1;	
					}
				}
				else {
					global.matrix_xm_yp[2][act_matrix_index] = 1;
					if (j <= 0) {
						global.matrix_xm_yp[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xm_yp[3][act_matrix_index] = 1;	
					}
				}
				act_matrix_index++;
			}
		}
	}

    // xp_ym

    act_matrix_index = 0;

	for (i=-1; i <= 2; i++) {
		for (j=-1; j <= 0; j++) {
			global.matrix_xp_ym[0][act_matrix_index] = i;
			global.matrix_xp_ym[1][act_matrix_index] = j;
			if (i <= 0) {
				global.matrix_xp_ym[2][act_matrix_index] = -1;
				if (j < 0) {
					global.matrix_xp_ym[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xp_ym[3][act_matrix_index] = 1;	
				}
			}
			else {
				global.matrix_xp_ym[2][act_matrix_index] = 1;
				if (j < 0) {
					global.matrix_xp_ym[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xp_ym[3][act_matrix_index] = 1;	
				}
			}
			act_matrix_index++;
		}
	}

	for (i=0; i <= 1; i++) {
		for (j=-2; j <= 1; j++) {
			for (k=0; k<act_matrix_index; k++) {
				if ((global.matrix_xp_ym[0][k] == i) && (global.matrix_xp_ym[1][k] == j)) {
					break;
				}
			}
			if (k == act_matrix_index) {
				global.matrix_xp_ym[0][act_matrix_index] = i;
				global.matrix_xp_ym[1][act_matrix_index] = j;
				if (i <= 0) {
					global.matrix_xp_ym[2][act_matrix_index] = -1;
					if (j < 0) {
						global.matrix_xp_ym[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xp_ym[3][act_matrix_index] = 1;	
					}
				}
				else {
					global.matrix_xp_ym[2][act_matrix_index] = 1;
					if (j < 0) {
						global.matrix_xp_ym[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xp_ym[3][act_matrix_index] = 1;	
					}
				}
				act_matrix_index++;
			}
		}
	}

    // xp_yp

    act_matrix_index = 0;

	for (i=-1; i <= 2; i++) {
		for (j=0; j <= 1; j++) {
			global.matrix_xp_yp[0][act_matrix_index] = i;
			global.matrix_xp_yp[1][act_matrix_index] = j;
			if (i <= 0) {
				global.matrix_xp_yp[2][act_matrix_index] = -1;
				if (j <= 0) {
					global.matrix_xp_yp[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xp_yp[3][act_matrix_index] = 1;	
				}
			}
			else {
				global.matrix_xp_yp[2][act_matrix_index] = 1;
				if (j <= 0) {
					global.matrix_xp_yp[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xp_yp[3][act_matrix_index] = 1;	
				}
			}
			act_matrix_index++;
		}
	}

	for (i=0; i <= 1; i++) {
		for (j=-1; j <= 2; j++) {
			for (k=0; k<act_matrix_index; k++) {
				if ((global.matrix_xp_yp[0][k] == i) && (global.matrix_xp_yp[1][k] == j)) {
					break;
				}
			}
			if (k == act_matrix_index) {
				global.matrix_xp_yp[0][act_matrix_index] = i;
				global.matrix_xp_yp[1][act_matrix_index] = j;
				if (i <= 0) {
					global.matrix_xp_yp[2][act_matrix_index] = -1;
					if (j <= 0) {
						global.matrix_xp_yp[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xp_yp[3][act_matrix_index] = 1;	
					}
				}
				else {
					global.matrix_xp_yp[2][act_matrix_index] = 1;
					if (j <= 0) {
						global.matrix_xp_yp[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xp_yp[3][act_matrix_index] = 1;	
					}
				}
				act_matrix_index++;
			}
		}
	}

    // xm_ym

    act_matrix_index = 0;

	for (i=-2; i <= 1; i++) {
		for (j=-1; j <= 0; j++) {
			global.matrix_xm_ym[0][act_matrix_index] = i;
			global.matrix_xm_ym[1][act_matrix_index] = j;
			if (i < 0) {
				global.matrix_xm_ym[2][act_matrix_index] = -1;
				if (j < 0) {
					global.matrix_xm_ym[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xm_ym[3][act_matrix_index] = 1;	
				}
			}
			else {
				global.matrix_xm_ym[2][act_matrix_index] = 1;
				if (j < 0) {
					global.matrix_xm_ym[3][act_matrix_index] = -1;
				}
				else {
					global.matrix_xm_ym[3][act_matrix_index] = 1;	
				}
			}
			act_matrix_index++;
		}
	}

	for (i=-1; i <= 0; i++) {
		for (j=-2; j <= 1; j++) {
			for (k=0; k<act_matrix_index; k++) {
				if ((global.matrix_xm_ym[0][k] == i) && (global.matrix_xm_ym[1][k] == j)) {
					break;
				}
			}
			if (k == act_matrix_index) {
				global.matrix_xm_ym[0][act_matrix_index] = i;
				global.matrix_xm_ym[1][act_matrix_index] = j;
				if (i < 0) {
					global.matrix_xm_ym[2][act_matrix_index] = -1;
					if (j < 0) {
						global.matrix_xm_ym[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xm_ym[3][act_matrix_index] = 1;	
					}
				}
				else {
					global.matrix_xm_ym[2][act_matrix_index] = 1;
					if (j < 0) {
						global.matrix_xm_ym[3][act_matrix_index] = -1;
					}
					else {
						global.matrix_xm_ym[3][act_matrix_index] = 1;	
					}
				}
				act_matrix_index++;
			}
		}
	}

	// for (i=0; i<global.active_grid_neighbours; i++) {
	// 	printf("%d %d\t%d %d\n", global.matrix_xm_yp[0][i], global.matrix_xm_yp[1][i], global.matrix_xm_yp[2][i], global.matrix_xm_yp[3][i]);
	// }
	// printf("xm_yp\n");

	// for (i=0; i<global.active_grid_neighbours; i++) {
	// 	printf("%d %d\t%d %d\n", global.matrix_xp_ym[0][i], global.matrix_xp_ym[1][i], global.matrix_xp_ym[2][i], global.matrix_xp_ym[3][i]);
	// }
	// printf("xp_ym\n");

	// for (i=0; i<global.active_grid_neighbours; i++) {
	// 	printf("%d %d\t%d %d\n", global.matrix_xp_yp[0][i], global.matrix_xp_yp[1][i], global.matrix_xp_yp[2][i], global.matrix_xp_yp[3][i]);
	// }
	// printf("xp_yp\n");

	// for (i=0; i<global.active_grid_neighbours; i++) {
	// 	printf("%d %d\t%d %d\n", global.matrix_xm_ym[0][i], global.matrix_xm_ym[1][i], global.matrix_xm_ym[2][i], global.matrix_xm_ym[3][i]);
	// }
	// printf("xm_ym\n");

}

void initialize_active_grid() {

	int i, j;

	global.active_grid_Nx = (int) floor(global.SX / global.active_grid_division) + 1;
	global.active_grid_Ny = (int) floor(global.SY / global.active_grid_division) + 1;

	global.active_grid_division_half = (global.active_grid_division / 2.0);

	printf("Active grid size %d x %d\n", global.active_grid_Nx, global.active_grid_Ny);

	global.active_grid = (double**) malloc(global.active_grid_Nx * sizeof(double*));

    for (i=0; i<global.active_grid_Nx; i++)
    {
        global.active_grid[i] = (double*) malloc(global.active_grid_Ny * sizeof(double));
        for (j=0; j<global.active_grid_Ny; j++)
        {
            global.active_grid[i][j] = Rand();
        }
    }
}

void calculate_PBC_folded_distance(double *dx, double *dy,double x1,double y1,double x2,double y2)
{
    double dxtemp,dytemp;

    dxtemp = x2 - x1;
    dytemp = y2 - y1;

    #ifndef non_periodic_boundary
    //Periodic Boundary Conditions Fold-Back
    if (dxtemp>global.SX2) dxtemp -= global.SX;
    if (dxtemp<-global.SX2) dxtemp += global.SX;
    if (dytemp>global.SY2) dytemp -= global.SY;
    if (dytemp<-global.SY2) dytemp += global.SY;
    #endif
    
    *dx = dxtemp;
    *dy = dytemp;
}

void calculate_external_forces()
{
    int i;

    for(i=0;i<global.N_particles;i++)
        {
        	if (global.particle_color[i] == 2) {
        		global.particle_fx[i] += 0.5;
        	}
        	else {
        		global.particle_fx[i] -= 0.5;
        	}
        }
}

void calculate_active_grid_forces() {

	int i, k;
	int center_grid_i, center_grid_j;
	double center_grid_x, center_grid_y;
	double act_particle_fx, act_particle_fy;
    int act_i, act_j;

    global.nr_moving_partices = 0;

	for(i=0;i<global.N_particles;i++) {

        center_grid_i = (int) floor(global.particle_x[i] / global.active_grid_division);
        center_grid_j = (int) floor(global.particle_y[i] / global.active_grid_division);

        center_grid_x = center_grid_i * global.active_grid_division + global.active_grid_division_half;
        center_grid_y = center_grid_j * global.active_grid_division + global.active_grid_division_half;

        if (global.particle_is_passive[i] == 0) {
		
    		// printf("%f %f\n", global.particle_x[i], global.particle_y[i]);

    		act_particle_fx = 0;
    		act_particle_fy = 0;

    		// printf("%d - %d | %f - %f\n", center_grid_i, center_grid_j, center_grid_x, center_grid_y);
    		// printf("-----\n");

    		if (global.particle_x[i] < center_grid_x) {
    			if (global.particle_y[i] < center_grid_y) {
    				for (k=0; k<global.active_grid_neighbours; k++) {

                        act_i = center_grid_i + global.matrix_xm_ym[0][k];
                        act_j = center_grid_j + global.matrix_xm_ym[1][k];

                        // PBC foldback
                        if (act_i < 0) act_i = global.active_grid_Nx + act_i;
                        else if (act_i >= global.active_grid_Nx) act_i = act_i - global.active_grid_Nx;
                        if (act_j < 0) act_j = global.active_grid_Ny + act_j;
                        else if (act_j >= global.active_grid_Ny) act_j = act_j - global.active_grid_Ny;

    					act_particle_fx += global.active_grid[act_i][act_j] * global.matrix_xm_ym[2][k];
    					act_particle_fy += global.active_grid[act_i][act_j] * global.matrix_xm_ym[3][k];

    					// absorption
    					global.active_grid[act_i][act_j] -= global.active_grid_absorption;
    					if (global.active_grid[act_i][act_j] < 0) {
                            global.active_grid[act_i][act_j] = 0;
                        }
    				}
    			}
    			else {
    				for (k=0; k<global.active_grid_neighbours; k++) {

                        act_i = center_grid_i + global.matrix_xm_yp[0][k];
                        act_j = center_grid_j + global.matrix_xm_yp[1][k];

                        // PBC foldback
                        if (act_i < 0) act_i = global.active_grid_Nx + act_i;
                        else if (act_i >= global.active_grid_Nx) act_i = act_i - global.active_grid_Nx;
                        if (act_j < 0) act_j = global.active_grid_Ny + act_j;
                        else if (act_j >= global.active_grid_Ny) act_j = act_j - global.active_grid_Ny;

    					act_particle_fx += global.active_grid[act_i][act_j] * global.matrix_xm_yp[2][k];
    					act_particle_fy += global.active_grid[act_i][act_j] * global.matrix_xm_yp[3][k];

    					// absorption
                        global.active_grid[act_i][act_j] -= global.active_grid_absorption;
                        if (global.active_grid[act_i][act_j] < 0) {
                            global.active_grid[act_i][act_j] = 0;
                        }
    				}
    			}
    		}
    		else {
    			if (global.particle_y[i] < center_grid_y) {
    				for (k=0; k<global.active_grid_neighbours; k++) {

                        act_i = center_grid_i + global.matrix_xp_ym[0][k];
                        act_j = center_grid_j + global.matrix_xp_ym[1][k];

                        // PBC foldback
                        if (act_i < 0) act_i = global.active_grid_Nx + act_i;
                        else if (act_i >= global.active_grid_Nx) act_i = act_i - global.active_grid_Nx;
                        if (act_j < 0) act_j = global.active_grid_Ny + act_j;
                        else if (act_j >= global.active_grid_Ny) act_j = act_j - global.active_grid_Ny;

    					act_particle_fx += global.active_grid[act_i][act_j] * global.matrix_xp_ym[2][k];
    					act_particle_fy += global.active_grid[act_i][act_j] * global.matrix_xp_ym[3][k];

    					// absorption
                        global.active_grid[act_i][act_j] -= global.active_grid_absorption;
                        if (global.active_grid[act_i][act_j] < 0) {
                            global.active_grid[act_i][act_j] = 0;
                        }
    				}

    			}
    			else {
    				for (k=0; k<global.active_grid_neighbours; k++) {

                        act_i = center_grid_i + global.matrix_xp_yp[0][k];
                        act_j = center_grid_j + global.matrix_xp_yp[1][k];

                        // PBC foldback
                        if (act_i < 0) act_i = global.active_grid_Nx + act_i;
                        else if (act_i >= global.active_grid_Nx) act_i = act_i - global.active_grid_Nx;
                        if (act_j < 0) act_j = global.active_grid_Ny + act_j;
                        else if (act_j >= global.active_grid_Ny) act_j = act_j - global.active_grid_Ny;

    					act_particle_fx += global.active_grid[act_i][act_j] * global.matrix_xp_yp[2][k];
    					act_particle_fy += global.active_grid[act_i][act_j] * global.matrix_xp_yp[3][k];

    					// absorption
                        global.active_grid[act_i][act_j] -= global.active_grid_absorption;
                        if (global.active_grid[act_i][act_j] < 0) {
                            global.active_grid[act_i][act_j] = 0;
                        }
    				}
    			}
    		}

    		// 	printf("time = %d, i = %d, x = %f, y = %f, fx = %f, fy = %f\n", global.time, i, global.particle_x[i], global.particle_y[i], act_particle_fx, act_particle_fy);

    		global.particle_fx[i] += act_particle_fx;
    		global.particle_fy[i] += act_particle_fy;

            if (( (global.particle_fx[i]*global.particle_fx[i]) + (global.particle_fy[i]*global.particle_fy[i]) ) < 0.0001) {
                global.particle_color[i] = 3;
            }
            else {
                global.particle_color[i] = 2;
                global.nr_moving_partices += 1;
            }

    		global.particle_angle_rad[i] = atan2(global.particle_fy[i],global.particle_fx[i]);

        }

		// exit(0);
	}
}

void calculate_active_grid_recovery() {

	int i, j;

	for (i=0; i<global.active_grid_Nx; i++)
    {
        for (j=0; j<global.active_grid_Ny; j++)
        {
            
            // recovery
            global.active_grid[i][j] += global.active_grid_recovery;
            if (global.active_grid[i][j] > 1) {
                global.active_grid[i][j] = 1;
            }
        }
    }

}

void calculate_internal_motor_forces()
{
    int i;
    double theta;

    for(i=0;i<global.N_particles;i++)
        {
            global.particle_fx[i] += global.particle_motor_force[i] * global.particle_cosfi[i];
            global.particle_fy[i] += global.particle_motor_force[i] * global.particle_sinfi[i];

            global.particle_motor_ellapsed_time[i]++;
            
            if (global.particle_motor_ellapsed_time[i] >= global.particle_motor_total_time[i])
                {
                    //generate a new angle
                    global.particle_angle_rad[i] = 2.0 * PI * Rand();
                    global.particle_cosfi[i] = cos(global.particle_angle_rad[i]);
                    global.particle_sinfi[i] = sin(global.particle_angle_rad[i]);
                    global.particle_motor_ellapsed_time[i] = 0;
                    global.particle_motor_total_time[i] = global.generic_particle_motor_minimum_time;
                    global.particle_motor_total_time[i] += (int)floor(Rand()*(global.generic_particle_motor_maximum_time - global.generic_particle_motor_minimum_time));
                    if (global.particle_cosfi[i]<0) global.particle_color[i] = 2;
                    else global.particle_color[i] = 3;
                }
            
        }
}

void calculate_interparticle_forces()
{
    int i,j,ii;
    double dx,dy,dr,delta_compressed;
    double fx,fy;

    for(ii=0;ii<global.N_Verlet;ii++)
        {
            i = global.Verlet_i[ii];
            j = global.Verlet_j[ii];
            calculate_PBC_folded_distance(&dx,&dy,global.particle_x[i],global.particle_y[i],global.particle_x[j],global.particle_y[j]);
            dr =  sqrt(dx*dx+dy*dy);
            delta_compressed = dr - global.particle_R[i] - global.particle_R[j];
            if (delta_compressed<0)
                {
                    fx = delta_compressed * global.generic_particle_k_spring * dx/dr;
                    fy = delta_compressed * global.generic_particle_k_spring * dy/dr;
                    
                    global.particle_fx[i] += fx;
                    global.particle_fy[i] += fy;
                    global.particle_fx[j] -= fx;
                    global.particle_fy[j] -= fy;
                }
        }

}


void interaction_particle_pinningsiteGridCell(int k, int i, int j)
{
    double dx, dy;
    int l, i_pinningsite;

    for (l=1; l<=global.pinningsite_grid[i][j][0]; l++)
    {
        i_pinningsite = global.pinningsite_grid[i][j][l];
        calculate_PBC_folded_distance(&dx, &dy, global.particle_x[k], global.particle_y[k], global.pinningsite_x[i_pinningsite], global.pinningsite_y[i_pinningsite]);

        if ((dx*dx + dy*dy) <= global.pinningsite_R[i_pinningsite] * global.pinningsite_R[i_pinningsite])
        {
            global.particle_fx[k] += global.pinningsite_fmax[i_pinningsite] * dx / global.pinningsite_R[i_pinningsite];
            global.particle_fy[k] += global.pinningsite_fmax[i_pinningsite] * dy / global.pinningsite_R[i_pinningsite];
        }
    }
}

void calculate_pinning_forces()
{
    int i, j, i_left, i_right, j_down, j_up, k;

    for (k=0; k<global.N_particles; k++)
    {
        i = (int) floor(global.particle_x[k] / global.pinning_grid_cell_size);
        j = (int) floor(global.particle_y[k] / global.pinning_grid_cell_size);

        i_left = i - 1;
        i_right = i + 1;
        j_down = j - 1;
        j_up = j + 1;

        // PBC foldback
        if (i_left < 0) i_left = global.pinningsite_grid_Nx - 1;
        else if (i_right >= global.pinningsite_grid_Nx) i_right = 0;
        if (j_down < 0) j_down = global.pinningsite_grid_Ny - 1;
        else if(j_up >= global.pinningsite_grid_Ny) j_up = 0;

        interaction_particle_pinningsiteGridCell(k, i, j);
        interaction_particle_pinningsiteGridCell(k, i, j_down);
        interaction_particle_pinningsiteGridCell(k, i, j_up);
        interaction_particle_pinningsiteGridCell(k, i_left, j);
        interaction_particle_pinningsiteGridCell(k, i_left, j_down);
        interaction_particle_pinningsiteGridCell(k, i_left, j_up);
        interaction_particle_pinningsiteGridCell(k, i_right, j);
        interaction_particle_pinningsiteGridCell(k, i_right, j_down);
        interaction_particle_pinningsiteGridCell(k, i_right, j_up);
    }
}

#ifdef drive_forces_type
void calculate_drive_force()
{
    int i;

    for (i=0; i<global.N_particles; i++)
    {
        #if drive_forces_type == 0 || drive_forces_type == 3
        // perpendicular
        global.particle_fx[i] += global.drive_force;
        #elif drive_forces_type == 1 || drive_forces_type == 4
        // same as gradient
        global.particle_fy[i] += global.drive_force;
        #elif drive_forces_type == 2 || drive_forces_type == 5
        // opposite
        global.particle_fy[i] -= global.drive_force;
        #endif
    }
}
#endif

#ifdef non_periodic_boundary
    void calculate_interaction_with_wall()
    {
        int i;

        for (i=0; i<global.N_particles; i++)
        {
            if (global.SX - global.particle_x[i] <= global.particle_R[i])
            {
                global.particle_fx[i] += 2*(global.SX - global.particle_x[i] - global.particle_R[i]) * global.generic_particle_k_spring;
            }
            else if (global.particle_x[i] <= global.particle_R[i])
            {
                global.particle_fx[i] += 2*(global.particle_R[i] - global.particle_x[i]) * global.generic_particle_k_spring;
            }
            
            if (global.SY - global.particle_y[i] <= global.particle_R[i])
            {
                global.particle_fy[i] += 2*(global.SY - global.particle_y[i] - global.particle_R[i]) * global.generic_particle_k_spring;
            }
            else if (global.particle_y[i] <= global.particle_R[i])
            {
                global.particle_fy[i] += 2*(global.particle_R[i] - global.particle_y[i]) * global.generic_particle_k_spring;
            }
        }
    }
#endif

void move_particles()
{
    int i;
    double dx,dy;
    double act_x_moving_in_passive_particle = 0.0;
    double act_y_moving_in_passive_particle = 0.0;
    double act_x_moving_in_moving_particle = 0.0;
    double act_y_moving_in_moving_particle = 0.0;
    double act_x_moving_in_standing_particle = 0.0;
    double act_y_moving_in_standing_particle = 0.0;

    double number_of_passive = 0.0;
    double number_of_moving = 0.0;
    double number_of_standing = 0.0;

    for(i=0;i<global.N_particles;i++)
    {

            //Brownian (overdamped) motion
            dx = global.particle_fx[i] * global.dt;
            dy = global.particle_fy[i] * global.dt;

            //calculate Verlet criteria displacements
            global.particle_dx_so_far[i] += dx;
            global.particle_dy_so_far[i] += dy; 
            check_Verlet_need_to_rebuild(i);

            global.particle_x[i] += dx;
            global.particle_y[i] += dy;

            if (global.particle_is_passive[i]) {
                act_x_moving_in_passive_particle += dx;
                act_y_moving_in_passive_particle += dy;
                number_of_passive += 1;
            }
            else {
                if (global.particle_color[i] == 2) {
                    act_x_moving_in_moving_particle += dx;
                    act_y_moving_in_moving_particle += dy;
                    number_of_moving += 1;
                }
                else {
                    act_x_moving_in_standing_particle += dx;
                    act_y_moving_in_standing_particle += dy;
                    number_of_standing += 1;
                }
            }

            //place back into PBC box
            if (global.particle_x[i]>global.SX) global.particle_x[i] -= global.SX;
            if (global.particle_y[i]>global.SY) global.particle_y[i] -= global.SY;
            if (global.particle_x[i]<= 0.0) global.particle_x[i] += global.SX;
            if (global.particle_y[i]<= 0.0) global.particle_y[i] += global.SY;

    }

    if (number_of_passive != 0) {
        global.average_x_moving_in_passive_particle += act_x_moving_in_passive_particle/number_of_passive;
        global.average_y_moving_in_passive_particle += act_y_moving_in_passive_particle/number_of_passive;
    }
    else {
        global.average_x_moving_in_passive_particle += 0;
        global.average_y_moving_in_passive_particle += 0;

    }

    if (number_of_moving != 0) {        
        global.average_x_moving_in_moving_particle += act_x_moving_in_moving_particle/number_of_moving;
        global.average_y_moving_in_moving_particle += act_y_moving_in_moving_particle/number_of_moving;
    }
    else {
        global.average_x_moving_in_moving_particle += 0;
        global.average_y_moving_in_moving_particle += 0;
    }

    if (number_of_standing != 0) {
        global.average_x_moving_in_standing_particle += act_x_moving_in_standing_particle/number_of_standing;
        global.average_y_moving_in_standing_particle += act_y_moving_in_standing_particle/number_of_standing;
    }
    else {
        global.average_x_moving_in_standing_particle += 0;
        global.average_y_moving_in_standing_particle += 0;

    }
    

    //zero forces
    for(i=0;i<global.N_particles;i++)
        {
             global.particle_fx[i] = 0.0;
             global.particle_fy[i] = 0.0;
        }

}

void check_Verlet_need_to_rebuild(int i)
{
    double dr2;

    dr2 =  global.particle_dx_so_far[i] * global.particle_dx_so_far[i] + global.particle_dy_so_far[i] * global.particle_dy_so_far[i];
    if (dr2 > global.Verlet_distance_to_rebuild * global.Verlet_distance_to_rebuild) global.flag_Verlet_needs_rebuild = 1;
}



void rebuild_Verlet()
{
    int i,j;
    double dx,dy;
    double dr2;
    double Verlet_R2;

    //discs are in contact at a distance less than 2R
    //the distance to which the Verlet list goes out is 4R 
    Verlet_R2 = global.Verlet_R * global.Verlet_R;

    global.N_Verlet = 0;
    global.Verlet_i = (int *) realloc(global.Verlet_i,global.N_Verlet*sizeof(int));
    global.Verlet_j = (int *) realloc(global.Verlet_j,global.N_Verlet*sizeof(int));


    for(i=0;i<global.N_particles;i++)
        for(j=i+1;j<global.N_particles;j++)
            {
                calculate_PBC_folded_distance(&dx,&dy,global.particle_x[i],global.particle_y[i],global.particle_x[j],global.particle_y[j]);
                dr2 = dx*dx + dy*dy;
                if (dr2<Verlet_R2)
                    {
                        global.N_Verlet++;
                        global.Verlet_i = (int *) realloc(global.Verlet_i,global.N_Verlet*sizeof(int));
                        global.Verlet_j = (int *) realloc(global.Verlet_j,global.N_Verlet*sizeof(int));
                        global.Verlet_i[ global.N_Verlet-1] = i;
                        global.Verlet_j[ global.N_Verlet-1] = j;
                    }
            }

    //test write
    /*
    for(i=0;i<global.N_Verlet;i++)
        {
            printf("%d %d\n",global.Verlet_i[i],global.Verlet_j[i]);
        }
    */
    //printf("Verlet rebuilt at time %d list lenght = %d\n",global.time,global.N_Verlet);

    //refresh the distance couters
    for(i=0;i<global.N_particles;i++)
        {
            global.particle_dx_so_far[i] = 0.0;
            global.particle_dy_so_far[i] = 0.0;
        }
    global.flag_Verlet_needs_rebuild = 0;
}

void check_affine_displacement()
{
    int i,j,ii,jj;
    double dx,dy;
    double dr2;
    double threshold_l2;
    int counter;
    double act_average_neighbours_inside_threshold = 0.0;

    threshold_l2 = global.threshold_l * global.threshold_l;


    for(i=0;i<global.N_particles;i++)
    {
        global.neighbours_inside_threshold_tmp_N[i] = 0;
        
        for(j=0;j<global.N_particles;j++)
        {
            calculate_PBC_folded_distance(&dx,&dy,global.particle_x[i],global.particle_y[i],global.particle_x[j],global.particle_y[j]);
            dr2 = dx*dx + dy*dy;
            if (dr2<threshold_l2)
            {
                global.neighbours_inside_threshold_tmp_N[i]++;
                global.neighbours_inside_threshold_tmp[i] = (int *) realloc(global.neighbours_inside_threshold_tmp[i],global.neighbours_inside_threshold_tmp_N[i]*sizeof(int));
                global.neighbours_inside_threshold_tmp[i][global.neighbours_inside_threshold_tmp_N[i]-1] = j;
            }
        }

        counter = 0;

        for(ii=0;ii<global.neighbours_inside_threshold_N[i];ii++)
        {
            for(jj=0;jj<global.neighbours_inside_threshold_tmp_N[i];jj++)
            {
                if(global.neighbours_inside_threshold[i][ii]==global.neighbours_inside_threshold_tmp[i][jj])
                {
                    counter++;
                }
            }
        }

        act_average_neighbours_inside_threshold += counter;

        global.neighbours_inside_threshold_N[i] = global.neighbours_inside_threshold_tmp_N[i];
        global.neighbours_inside_threshold[i] = (int *) realloc(global.neighbours_inside_threshold[i],global.neighbours_inside_threshold_N[i]*sizeof(int));
        
        for(ii=0;ii<global.neighbours_inside_threshold_N[i];ii++)
        {
            global.neighbours_inside_threshold[i][ii] = global.neighbours_inside_threshold_tmp[i][ii];
        }

    }

    global.average_neighbours_inside_threshold += (act_average_neighbours_inside_threshold/global.N_particles);
}

void clusterize_particles()
{
    int i, j, k, cluster_i, cluster_j;
    double dx, dy;
    struct particle_struct* temp;

    for (i=0; i<global.N_particles; i++)
    {
        global.particles[i].next = NULL;
        global.particles[i].clusterID = i;
        global.particle_cluster[i].N = 1;
        global.particle_cluster[i].first = global.particles+i;
        global.particle_cluster[i].last = global.particles+i;
    }

    for (k=0; k<global.N_Verlet; k++)
    {
        i = global.Verlet_i[k];
        j = global.Verlet_j[k];

        if (global.particle_is_passive[i] && global.particle_is_passive[j]) {

            calculate_PBC_folded_distance(&dx, &dy, global.particle_x[i], global.particle_y[i], global.particle_x[j], global.particle_y[j]);
            if ((global.particles[i].clusterID != global.particles[j].clusterID) && (dx*dx+dy*dy <= (global.particle_R[i]+global.particle_R[j])*(global.particle_R[i]+global.particle_R[j])))
            {
                if (global.particle_cluster[global.particles[i].clusterID].N > global.particle_cluster[global.particles[j].clusterID].N)
                {
                    cluster_i = global.particles[i].clusterID;
                    cluster_j = global.particles[j].clusterID;
                }
                else
                {
                    cluster_i = global.particles[j].clusterID;
                    cluster_j = global.particles[i].clusterID;
                }

                temp = global.particle_cluster[cluster_j].first;
                while (temp)
                {
                    temp->clusterID = cluster_i;
                    temp = temp->next;
                }

                global.particle_cluster[cluster_i].N += global.particle_cluster[cluster_j].N;
                global.particle_cluster[cluster_i].last->next = global.particle_cluster[cluster_j].first;
                global.particle_cluster[cluster_i].last = global.particle_cluster[cluster_j].last;

                global.particle_cluster[cluster_j].N = 0;
                global.particle_cluster[cluster_j].first = NULL;
                global.particle_cluster[cluster_j].last = NULL;
            }
        }
    }
}

int largest_clusterID()
{
    int i, largest_cluster;

    largest_cluster = 0;

    while (!global.particle_is_passive[largest_cluster]) {
        largest_cluster++;
    }

    for (i=largest_cluster+1; i<global.N_particles; i++)
    {
        if (global.particle_is_passive[i] && (global.particle_cluster[largest_cluster].N < global.particle_cluster[i].N))
        {
            largest_cluster = i;
        }
    }
    
    return largest_cluster;
}

int number_of_clusters() {

    int i, counter;

    counter = 0;

    for (i=0; i<global.N_particles; i++)
    {
        if (global.particle_is_passive[i] && (global.particle_cluster[i].N > 0)) {
            counter++;
        }
    }
    
    return counter;
}

void centroid_of_cluster(double* x_com, double* y_com, int clusterID)
{
    int i;
    double *x_rad, *y_rad, x_xi, x_zeta, y_xi, y_zeta;
    struct particle_struct* temp;

    x_rad = (double*) malloc(global.particle_cluster[clusterID].N * sizeof(double));
    y_rad = (double*) malloc(global.particle_cluster[clusterID].N * sizeof(double));

    temp = global.particle_cluster[clusterID].first;
    i = 0;
    while (temp)
    {
        x_rad[i] = global.particle_x[temp->i] / global.SX2 * PI;
        y_rad[i] = global.particle_y[temp->i] / global.SY2 * PI;
        temp = temp->next;
        i++;
    }

    x_xi = 0.;
    x_zeta = 0.;
    y_xi = 0.;
    y_zeta = 0.;
    for (i=0; i<global.particle_cluster[clusterID].N; i++)
    {
        x_xi += cos(x_rad[i]);
        x_zeta += sin(x_rad[i]);
        y_xi += cos(y_rad[i]);
        y_zeta += sin(y_rad[i]);
    }
    x_xi /= global.particle_cluster[clusterID].N;
    x_zeta /= global.particle_cluster[clusterID].N;
    y_xi /= global.particle_cluster[clusterID].N;
    y_zeta /= global.particle_cluster[clusterID].N;

    *x_com = global.SX * (atan2(-x_zeta, -x_xi) / 2. / PI + 0.5);
    *y_com = global.SY * (atan2(-y_zeta, -y_xi) / 2. / PI + 0.5);

    free(x_rad);
    free(y_rad);
}

void calculate_average_distance_from_centroid_of_cluster() {
    
    int largest_cluster;
    double x_com, y_com;
    int i;
    double dx, dy, sum_dr;
    struct particle_struct* temp;

    largest_cluster = largest_clusterID();

    centroid_of_cluster(&x_com, &y_com, largest_cluster);

    temp = global.particle_cluster[largest_cluster].first;
    i = 0;
    while (temp)
    {
        calculate_PBC_folded_distance(&dx, &dy, x_com, y_com, global.particle_x[temp->i], global.particle_y[temp->i]);
        
        sum_dr += sqrt(dx*dx+dy*dy);

        temp = temp->next;
        i++;
    }

    global.largest_x_com = x_com;
    global.largest_y_com = y_com;
    global.largest_average_dr = sum_dr/i;
}

void write_particles()
{
    FILE *testout;
    int i;
    char fname[300];

    sprintf(fname, "particles/%s.txt", global.baseName);
    testout = fopen(fname,"wt");
    for(i=0;i<global.N_particles;i++)
        fprintf(testout,"%lf %lf\n",global.particle_x[i],global.particle_y[i]);
    fclose(testout);
     printf("%d Particles written to testfile\n",global.N_particles);
}

void write_pinningsites()
{
    FILE *testout;
    int i;
    char fname[300];

    sprintf(fname, "pinningsites/%s.txt", global.baseName);
    testout = fopen(fname,"wt");
    for(i=0;i<global.N_pinningsites;i++)
        fprintf(testout,"%lf %lf\n",global.pinningsite_x[i],global.pinningsite_y[i]);
    fclose(testout);
     printf("%d Pinning sites written to testfile\n",global.N_pinningsites);
}

void density_distribution_to_file(int time_step, int N_bands)
{
    double y_width;
    int i, j, *bands;
    char fname[300];
    FILE* f;

    sprintf(fname, "density_distribution/%s_tstep_%d.dat", global.baseName, time_step);
    f = fopen(fname, "wt");

    y_width = global.SY / N_bands;
    bands = (int*) malloc(N_bands * sizeof(int));
    for (i=0; i<N_bands; i++)
    {
        bands[i] = 0.;
    }

    for (i=0; i<global.N_particles; i++)
    {
        bands[(int) floor(global.particle_y[i] / y_width)] ++; 
    }

    for (i=0; i<N_bands; i++)
    {
        fprintf(f, "%lf %lf\n", i * y_width, (double) bands[i] / global.N_particles);
    }

    free(bands);
    fclose(f);
}

void average_density_distribution(double* bands, int N_bands)
{
    double y_width;
    int i, j;

    y_width = global.SY / N_bands;

    for (i=0; i<global.N_particles; i++)
    {
        bands[(int) floor(global.particle_y[i] / y_width)] ++;
    }
}

void average_density_distribution_to_file(double* bands, int N_bands)
{
    double N_total, y_width = global.SY / N_bands;
    int i;
    char fname[300];
    FILE* f;

    sprintf(fname, "density_distribution/avg_%s.dat", global.baseName);
    f = fopen(fname, "wt");

    N_total = 0;
    for (i=0; i<N_bands; i++)
    {
        N_total += bands[i];
    }

    for (i=0; i<N_bands; i++)
    {
        fprintf(f, "%lf %lf\n", i * y_width, bands[i] / N_total);
    }

    fclose(f);
}

// centroid of largest cluster in time, to file
void centroid_in_time(FILE* f, FILE* g)
{
    double x_com, y_com;

    centroid_of_cluster(&x_com, &y_com, largest_clusterID());

    fprintf(f, "%d %lf\n", global.time, x_com);fflush(f);
    fprintf(g, "%d %lf\n", global.time, y_com);fflush(g);
}

void print_cluster(int clusterID)
{
    FILE* f;
    char fname[300];
    struct particle_struct* temp;

    sprintf(fname, "clusters/%s_ID%d.dat", global.baseName, clusterID);
    f = fopen(fname, "w");

    temp = global.particle_cluster[clusterID].first;
    while (temp)
    {
        fprintf(f, "%lf %lf\n", global.particle_x[temp->i], global.particle_y[temp->i]);
        temp = temp->next;
    }

    fclose(f);
}

void put_every_particle_in_same_cluster()
{
    int i;

    global.particle_cluster[0].N = global.N_particles;
    global.particle_cluster[0].first = global.particles;
    global.particle_cluster[0].last = global.particles+global.N_particles-1;
    global.particles[0].clusterID = 0;

    for (i=1; i<global.N_particles; i++)
    {
        global.particles[i-1].next = global.particles+i;
        global.particles[i].clusterID = 0;
        global.particle_cluster[i].N = 0;
        global.particle_cluster[i].first = NULL;
        global.particle_cluster[i].last = NULL;
    }

    global.particles[i-1].next = NULL;
}

void run_simulation()
{
    // FILE *f, *g;
    // char centroid_fname[300];
    clock_t t_start = clock();

    // int N_bands = 50;
    // double* bands = (double*) calloc(N_bands, sizeof(double));

    // sprintf(centroid_fname, "clusters/com_in_time/%s_x.dat", global.baseName);
    // f = fopen(centroid_fname, "w");
    // sprintf(centroid_fname, "clusters/com_in_time/%s_y.dat", global.baseName);
    // g = fopen(centroid_fname, "w");
    
    for(global.time=0;global.time < global.total_runtime;global.time++)
        {
            // calculate_internal_motor_forces();
            calculate_interparticle_forces();
            // calculate_pinning_forces();
            // calculate_external_forces();

            calculate_active_grid_forces();
            calculate_active_grid_recovery();

            // #ifdef drive_forces_type
            //     #if drive_forces_type >= 3 && drive_forces_type <= 5
            //         if (global.time%80000 == 0) global.drive_force = -global.drive_force;
            //     #endif
            //     calculate_drive_force();
            // #endif

            // #ifdef non_periodic_boundary
            //     calculate_interaction_with_wall();
            // #endif

            if (global.time%global.echo_time==0)
            {
                double percent = (double) global.time / global.total_runtime * 100.0;
                printf("\r%d/%d (%.2lf%%; ETA: %.2lf s)               ", global.time,global.total_runtime, percent, (100.0-percent)*((double)(clock()-t_start)/CLOCKS_PER_SEC)/percent);
                fflush(stdout);
            }
            if (global.time%global.movie_time==0) write_cmovie_frame();
            if (global.time%global.stat_time==0) {
                clusterize_particles();
                // calculate_average_distance_from_centroid_of_cluster();
                // check_affine_displacement();
                write_stat_data();
                global.average_x_moving_in_passive_particle = 0.0;
                global.average_y_moving_in_passive_particle = 0.0;
                global.average_x_moving_in_moving_particle = 0.0;
                global.average_y_moving_in_moving_particle = 0.0;
                global.average_x_moving_in_standing_particle = 0.0;
                global.average_y_moving_in_standing_particle = 0.0;

                // global.average_neighbours_inside_threshold = 0.0;
            }
            // if (global.time%100000 == 0) density_distribution_to_file(global.time, N_bands);
            // if (global.time>=300000)
            // {
            //     average_density_distribution(bands, N_bands);
            //     if (global.time%1000 == 0)
            //     {
            //         clusterize_particles();
            //         // put_every_particle_in_same_cluster();
            //         centroid_in_time(f, g);
            //     }
            // }

            move_particles();
            if (global.flag_Verlet_needs_rebuild==1) rebuild_Verlet();
        }
    printf("\r%d/%d (100.00%%; ETA: 0.00 s)\n",global.time,global.total_runtime);
    printf("                     in %.2lf s\r", (double)(clock()-t_start)/CLOCKS_PER_SEC);

    // average_density_distribution_to_file(bands, N_bands);
    // free(bands);
    // fclose(f);
    // fclose(g);
}

void write_cmovie_frame()
{
    int i,j;
    float floatholder;
    int intholder;
    int counter;

    intholder = global.N_particles + ((global.active_grid_Nx) * (global.active_grid_Ny));
    fwrite(&intholder,sizeof(int),1,global.movie_file);

    intholder = global.time;
    fwrite(&intholder,sizeof(int),1,global.movie_file);

    for (i=0;i<global.N_particles;i++)
	{
        intholder = global.particle_color[i];
        //intholder = global.particle_cluster[i]+2;
	    fwrite(&intholder,sizeof(int),1,global.movie_file);
	    intholder = i; //ID
        //(int)floor(10*global.particle_R[i]);
        //ID replaced by radius integer
	    fwrite(&intholder,sizeof(int),1,global.movie_file);
	    floatholder = (float)global.particle_x[i];
	    fwrite(&floatholder,sizeof(float),1, global.movie_file);
	    floatholder = (float)global.particle_y[i];
	    fwrite(&floatholder,sizeof(float),1, global.movie_file);
	    floatholder = global.particle_angle_rad[i]; 
        //former cump disp, now using angle
	    fwrite(&floatholder,sizeof(float),1,global.movie_file);
	}

	counter = global.N_particles;

    for (i=0;i<global.active_grid_Nx;i++)
	{
		for (j=0;j<global.active_grid_Ny;j++)
		{
	        intholder = 5 + (int)(global.active_grid[i][j]*10);//special color for pins
	        //intholder = global.particle_cluster[i]+2;
		    fwrite(&intholder,sizeof(int),1,global.movie_file);
		    intholder = counter++; //ID
	        //(int)floor(10*global.particle_R[i]);
	        //ID replaced by radius integer
		    fwrite(&intholder,sizeof(int),1,global.movie_file);
		    floatholder = (float) (i * global.active_grid_division + global.active_grid_division_half);
		    fwrite(&floatholder,sizeof(float),1, global.movie_file);
		    floatholder = (float) (j * global.active_grid_division + global.active_grid_division_half);
		    fwrite(&floatholder,sizeof(float),1, global.movie_file);
		    floatholder = 0.0;
	        //former cump disp, now using angle
		    fwrite(&floatholder,sizeof(float),1,global.movie_file);
		}
	}
}

void write_stat_data() {
    fprintf(global.stat_file, "%d %d %d %d %d %f %f %f %f %f %f\n", global.time, global.nr_non_passive_particles, global.nr_moving_partices, 
        number_of_clusters(), global.particle_cluster[largest_clusterID()].N, 
        global.average_x_moving_in_passive_particle, global.average_y_moving_in_passive_particle,
        global.average_x_moving_in_moving_particle, global.average_y_moving_in_moving_particle,
        global.average_x_moving_in_standing_particle, global.average_y_moving_in_standing_particle);
    fflush(global.stat_file);
}



int main(int argc,char *argv[])
{
    char cmd[500];

    printf("Active Matter Simulation\n");
    initialize_from_file(argc,argv);
    // #if dist_type >= -1 && dist_type <= 3
    //     initialize_pinningsites();
    // #elif dist_type == 4
    //     initialize_pinningsites_not_randomly();
    // #endif
    // initialize_pinning_grid();
    
	initialize_matrices();

    initialize_active_grid();

    initialize_particles();
    rebuild_Verlet();

    // write_particles();
    // write_pinningsites();

    run_simulation();
    printf("Simulation completed\n");

    write_particles(); fflush(stdout);
    // write_pinningsites();

    fclose(global.movie_file);
    fclose(global.stat_file);

    // sprintf(cmd, "povray_anim/create_frames %s %lf %lf 0", global.baseName, global.SX, global.SY);
    // if (system(cmd) == -1)
    // {
    //     printf("Could not run command: %s\n", cmd);
    // }

    return 0;
}
