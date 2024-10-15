#include <iostream>
#include <fstream> // file
#include "ceres/ceres.h"
#include "ceres/problem.h"
#include "BA.h"
#include "ba_L1_smoothing.h"
#include <ctime>   // clock
#include "matrix/matrix.h"
#include <string>
#include "sys/time.h"

using namespace std;

Eigen::Vector3f FitterLeastSquareMethod(std::vector<std::vector<double>> &X, std::vector<double> &Y, uint8_t orders)
{

    // abnormal input verification
//    if (X.size() < 2 || Y.size() < 2 || X.size() != Y.size() || orders < 1)
//        exit(EXIT_FAILURE);

//    std::cout<<X.size()<<" "<<Y.size()<<std::endl;
    Eigen::VectorXf result;
    std::vector<float> X_l;
    std::vector<float> Y_l;
    Eigen::MatrixXf A(X.size(), 3);
    Eigen::MatrixXf b(Y.size(), 1);
//    A<<1,2,3 ,2,4,5,6,7,8,4,5,6;
    // map sample data from STL vector to eigen vector
    for(int i=0;i<X.size();i++)
    {
        A.row(i)<<X[i][0],X[i][1],X[i][2];
        for(int j=0;j<3;j++)
        {
            X_l.push_back(X[i][j]);
//            A<<X[i][j];
        }
    }
    for(int i=0;i<Y.size();i++)
    {
        b.row(i)<<Y[i];
        Y_l.push_back(Y[i]);
    }

    Eigen::Vector3f x = A.colPivHouseholderQr().solve(b);
//    std::cout << "The solution is:\n" << x << std::endl;
//    std::cout<<" x:"<<x[0]<<" y:"<<x[1]<<" z:"<<x[2]<<std::endl;
//    std::cout << "Here is the matrix A:\n" << A << std::endl;
    Eigen::Map<Eigen::VectorXf> sampleX(X_l.data(), X_l.size());
    Eigen::Map<Eigen::VectorXf> sampleY(Y_l.data(), Y_l.size());
//        std::cout<<sampleX.size()<<" "<<sampleY.size()<<std::endl;
//    Eigen::Map<Eigen::VectorXf> sampleY(Y.data(), Y.size());
//
//    Eigen::MatrixXf mtxVandermonde(X.size(), orders + 1);  // Vandermonde matrix of X-axis coordinate vector of sample data
//    Eigen::VectorXf colVandermonde = sampleX;              // Vandermonde column
//
//    // construct Vandermonde matrix column by column
//    for (size_t i = 0; i < orders + 1; ++i)
//    {
//        if (0 == i)
//        {
//            mtxVandermonde.col(0) = Eigen::VectorXf::Constant(X.size(), 1, 1);
//            continue;
//        }
//        if (1 == i)
//        {
//            mtxVandermonde.col(1) = colVandermonde;
//            continue;
//        }
//        colVandermonde = colVandermonde.array()*sampleX.array();
//        mtxVandermonde.col(i) = colVandermonde;
//    }
//
//    // calculate coefficients vector of fitted polynomial
//    Eigen::VectorXf result = (mtxVandermonde.transpose()*mtxVandermonde).inverse()*(mtxVandermonde.transpose())*sampleY;
//
    return x;
}
util::point3d triangulate(std::vector<CameraParameters> cameras, std::vector<util::point2d> points)
{
//    for()
//    std::cout<<"camera_i:"<<camera_i.size()<<std::endl;
//    std::cout<<"point_i:"<<point_i.size()<<std::endl;
    int num_eqs = 2*points.size();
    int num_vars = 3;

    double *As = new double[num_eqs*num_vars];
    double *bs = new double[num_eqs];

    std::vector<std::vector<double>> A_l;
//    std::vector<std::vector<double>> b_l;
    std::vector<double> b_l;


    double *x = new double[num_vars];

    for(int i=0;i<points.size();i++)
    {
        double s, alpha, beta, gamma, t0, t1;

        std::vector<double > a1;
        std::vector<double > a2;

        s = cameras[i].s;
        alpha = cameras[i].alpha;
        beta = cameras[i].beta;
        gamma = cameras[i].gamma;
        t0 = cameras[i].t0;
        t1 = cameras[i].t1;

        double cos_alpha = cos(alpha);
        double sin_alpha = sin(alpha);
        double cos_beta = cos(beta);
        double sin_beta = sin(beta);
        double cos_gamma = cos(gamma);
        double sin_gamma = sin(gamma);

        double* A = As+6*i;
        double* b = bs+2*i;


        a1.push_back(cos_beta);
        a1.push_back(sin_alpha*sin_beta);
        a1.push_back(-cos_alpha*sin_beta);

        a2.push_back(0);
        a2.push_back(cos_alpha);
        a2.push_back(sin_alpha);

        A_l.push_back(a1);
        A_l.push_back(a2);


//        std::cout<<" "<<A[i]<<" ";

        b[0] = s*(cos_gamma*points[i].x+sin_gamma*points[i].y+t0);
        b[1] = s*(-sin_gamma*points[i].x+cos_gamma*points[i].y+t1);

        b_l.push_back(b[0]);
        b_l.push_back(b[1]);
    }

    std::vector<float> A_tmp;
    std::vector<float> B_tmp;

    for(int i=0;i<num_eqs*num_vars;i++)
    {
        A_tmp.push_back(As[i]);
    }
    for(int i=0;i<num_eqs;i++)
    {
        B_tmp.push_back(bs[i]);
    }


    Eigen::Vector3f result(FitterLeastSquareMethod(A_l, b_l, 3));
//
//    dgelsy_driver(As, bs, x, num_eqs, num_vars, 1);


    util::point3d tmp;



    tmp.x=result[0];
    tmp.y=result[1];
    tmp.z=result[2];
    return tmp;
}


void add_noise_uniform_cams(double *data,double ratio_noise,int ncams,double *data_with_noise){
    int i,j;
    double mean_temp;
    std::random_device rd;
    std::default_random_engine rng {rd()};
    for(i=0;i<6;i++){
        mean_temp = 0;
        for(j=0;j<ncams;j++){
            mean_temp += data[i+j*6];
        }
        mean_temp = mean_temp/ncams;
        if(mean_temp == 0) mean_temp = 0.1;
        //std::cout<<"Cam mean is "<<mean_temp<<std::endl;
        std::uniform_real_distribution<double> noise_data(-mean_temp*ratio_noise*0.01, mean_temp*ratio_noise*0.01); 
        for(j=0;j<ncams;j++) data_with_noise[i+j*6] = data[i+j*6]+noise_data(rng);
    }
}

void data_initial_new(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
    //Gauss noise
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();


    std::uniform_real_distribution<double> even_vmask(1, 100);
    std::uniform_real_distribution<double> range_vmask(5, 30);
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;

 
    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");


    
    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    add_noise_gauss_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");

 
    int count = 0;
    std::vector<util::point3d> point3d_vec;   
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }


 
    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++){
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


    //write_txt(p,ncams*6,"simulate_data_noOne/cams_noise.txt");

    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");


    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");

    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/p_real.txt");


    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}



void data_initial_exponential(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){

    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();

   
    std::uniform_real_distribution<double> even_vmask(1, 100);
    std::uniform_real_distribution<double> range_vmask(5, 30);
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;

 
    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");



    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    //add_noise_gauss_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    //add_noise_possion_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    add_noise_exponential_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");

 
    int count = 0;
    std::vector<util::point3d> point3d_vec;   
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }


   
    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++){
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


    //write_txt(p,ncams*6,"simulate_data_noOne/cams_noise.txt");

    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");


    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");

    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/p_real.txt");


    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}


void data_initial_uniform(int ncams,int n3Dpts,double value_min,double value_max,double ratio_outlier,double ratio_noise_cams,double range_img){
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();


    std::uniform_real_distribution<double> even_vmask(1, 100);
    std::uniform_real_distribution<double> range_vmask(5, 30);
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;


    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");



    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    add_noise_uniform_2dpts(real_2dpoints,value_min,value_max,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");


    int count = 0;
    std::vector<util::point3d> point3d_vec;   
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }



    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++)
    {
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");
    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");


    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");




    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}

void data_initial_kafang(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();


    std::uniform_real_distribution<double> even_vmask(1, 100);
    std::uniform_real_distribution<double> range_vmask(5, 30);
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;


    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");



    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    add_noise_kafang_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");


    int count = 0;
    std::vector<util::point3d> point3d_vec;   
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }



    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++)
    {
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");
    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");


    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");




    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}


void data_initial_possion(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();

 
    std::uniform_real_distribution<double> even_vmask(1, 100);
    std::uniform_real_distribution<double> range_vmask(5, 30);
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;

    
    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");


   
    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    //add_noise_gauss_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    add_noise_possion_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");

   
    int count = 0;
    std::vector<util::point3d> point3d_vec;   
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }


    
    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++){
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


    //write_txt(p,ncams*6,"simulate_data_noOne/cams_noise.txt");

    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");


    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");

    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/p_real.txt");


    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}





void data_initial_Laplace(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std:ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];
   
    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;
    
    fin.open("changedata_output/cams.txt", std::ios::in);
    for (i = 0; i < ncams*6; i++) fin >> cams[i];
    fin.close();
    fin.open("changedata_output/points3d.txt", std::ios::in);
    for (i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
    fin.close();
    fin.open("changedata_output/vmask.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();

    
    std::uniform_real_distribution<double> even_vmask(1, 100); 
    std::uniform_real_distribution<double> range_vmask(5, 30); 
    ratio_miss_pts = int(range_vmask(rng));
    n2dpts = 0;
        for(i=0;i<ncams*n3Dpts;i++) {
            vmask[i] = '1'-'1';
            if(even_vmask(rng)<100-ratio_miss_pts) {
                vmask[i] = '1'-'0';
                n2dpts++;
            }
        }
    n2dpts*=2;

    
    cams_noise = new double[ncams*6];
    add_noise_uniform_cams(cams,ratio_noise_cams,ncams,cams_noise);
    write_txt(cams_noise,ncams*6,"output_temp/cams_noise.txt");

    
    
    real_2dpoints = new double[n2dpts];
    pts2d_noise = new double[n2dpts];
    pts2d_noise_outlier = new double[n2dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    compute_real_2dpoints(p,vmask,ncams,n3Dpts,real_2dpoints);
    //std::cout<<"f is "<<f(p,0,real_2dpoints,vmask,weight_one,ncams,n3Dpts)<<std::endl;
    write_txt(real_2dpoints,n2dpts,"output_temp/real_2dpts.txt");
    //add_noise_gauss_2dpts(real_2dpoints,ratio_noise,n2dpts,pts2d_noise);
    add_noise_Laplace_2dpts(real_2dpoints,0,ratio_noise,n2dpts,pts2d_noise);
    write_txt(pts2d_noise,n2dpts,"output_temp/pts2d_noise.txt");

    num_outer = ratio_outlier*n2dpts*0.01/2;
    std::cout<<"add "<<num_outer<<" outliers"<<std::endl;
    index_outer = new int[num_outer];
    add_outer(pts2d_noise,range_img,n2dpts,num_outer,index_outer,pts2d_noise_outlier);
    write_txt(pts2d_noise_outlier,n2dpts,"output_temp/pts2d_outlier.txt");

    
    int count = 0;
    std::vector<util::point3d> point3d_vec;   //存储三角化后的3D点
    for(i=0;i<n3Dpts;i++){
        std::vector<CameraParameters> camera_i;
        std::vector<util::point2d> point_i;
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                camera_i.push_back(CameraParameters(cams[j*6],cams[j*6+1],cams[j*6+2],cams[j*6+3],cams[j*6+4],cams[j*6+5]));
                point_i.push_back(util::point2d(pts2d_noise_outlier[count*2],pts2d_noise_outlier[count*2+1]));
                count ++;
            }
        }
        util::point3d point3d_i = triangulate(camera_i, point_i);
        point3d_vec.push_back(point3d_i);
    }
    

    
    double one_avg_random,one_std_random;
    double *error_begin = new double[1];
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts; i++){
        pts3d_simulate[i*3] = point3d_vec[i].x;
        pts3d_simulate[i*3+1] = point3d_vec[i].y;
        pts3d_simulate[i*3+2] = point3d_vec[i].z;
    }
    for(i=0;i<ncams*6;i++) p[i] = cams[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


    //write_txt(p,ncams*6,"simulate_data_noOne/cams_noise.txt");
    
    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error no one is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"error_begin.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3]; 
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");

    
    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");
    
    for(i=0;i<ncams*6;i++) p[i] = cams_noise[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/p_real.txt");


    delete[] cams;
    delete[] pts3d;
    delete[] imgpts;
    delete[] vmask;
    delete[] cams_noise;
    delete[] real_2dpoints;
    delete[] pts2d_noise;
    delete[] pts2d_noise_outlier;
    delete[] p;
    delete[] pts3d_simulate;
    delete[] weight_one;
    delete[] ncams_vec;
    delete[] n3Dpts_vec;
    delete[] vec_gamma;
}

int main(int argc,char** argv){
    int method_index,ncams,n3Dpts,runing_times;
    double ratio_noise,ratio_outlier,ratio_noise_cams,range_img,mean_one,std_one,value_min,value_max;
    std:ifstream fin;
    
    ncams = 21;
    n3Dpts = 60;
    ratio_noise = 0.2*0.01;
    ratio_outlier = 5.0;
    ratio_noise_cams = 2.0;
    range_img = 1024;
    value_min = -8;
    value_max = 8;

    runing_times = std::atoi(argv[1]);
    
    
    data_initial_Laplace(ncams,n3Dpts,ratio_noise*range_img,ratio_outlier,ratio_noise_cams,range_img);
    
    int i,n2Dpts;
    double *mot,*pts2d,*mot_temp,*weight_one,*mot_new,*time_smoothing,*error_smoothing;
    char *vmask;
    std::clock_t start, end;

    time_smoothing = new double[1];
    error_smoothing = new double[1];
        
    mot = new double[ncams*6+n3Dpts*3];
    mot_temp = new double[ncams*6+n3Dpts*3];
    mot_new = new double[ncams*6+n3Dpts*3];
    fin.open("simulate_data/mot_random.txt", std::ios::in);
    for (i = 0; i < ncams*6+n3Dpts*3; i++) fin >> mot[i];
    fin.close();

    n2Dpts = 0;
    vmask = new char[n3Dpts*ncams];
    weight_one = new double[n3Dpts*ncams];
    fin.open("simulate_data/vmask_random.txt", std::ios::in);
    for (i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
        weight_one[i] = 1.0;
        if((int)vmask[i] == 1) n2Dpts++;
    }
    fin.close();
    n2Dpts *=2;

    pts2d = new double[n2Dpts];
    fin.open("simulate_data/pts2d_random.txt", std::ios::in);
    for (i = 0; i < n2Dpts; i++) fin >> pts2d[i];
    fin.close();

    std::cout<<"(*) Now run smoothing method:"<<std::endl;
    struct timeval t1,t2;

    double_copy(mot_temp,mot,ncams*6+n3Dpts*3);
    std::cout<<"initial error is "<<f(mot_temp,0,pts2d,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    gettimeofday(&t1,NULL);
    ba_L1_smoothing_grad(mot_temp,pts2d,vmask,ncams,n3Dpts,mot_new);
    gettimeofday(&t2,NULL);
    time_smoothing[0] = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;;
    error_smoothing[0] = f(mot_new,0,pts2d,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"Final L1 error of smoothing is "<<error_smoothing[0]<<", running time is "<<time_smoothing[0]<<"s"<<std::endl;
    write_txt_app<double>(time_smoothing,1,"output_temp/time_smoothing.txt");
    write_txt_app<double>(error_smoothing,1,"output_temp/error_smoothing.txt");

    double *p_recover,*pts2d_recover,*vec_gamma;
    p_recover = new double[ncams*6+n3Dpts*3];
    pts2d_recover = new double[n2Dpts];
    vec_gamma = new double[ncams];
    fin.open("simulate_data_noOne/one_avg.txt", std::ios::in);
    for (i = 0; i < runing_times+1; i++) fin >> mean_one;
    fin.close();
    fin.open("simulate_data_noOne/one_std.txt", std::ios::in);
    for (i = 0; i < runing_times+1; i++) fin >> std_one;
    fin.close();
    fin.open("simulate_data_noOne/vec_gamma.txt", std::ios::in);
    for (i = 0; i < ncams; i++) fin >> vec_gamma[i];
    fin.close();
    std::cout<<"avg is "<<mean_one<<", std is "<<std_one<<std::endl;

    Toone_recover_new(mot_new,pts2d,ncams,n3Dpts,n2Dpts,mean_one,std_one,vec_gamma,p_recover,pts2d_recover);
    error_smoothing[0] = f(p_recover,0,pts2d_recover,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"Final recovery L1 error is "<<error_smoothing[0]<<std::endl;
    write_txt_app<double>(error_smoothing,1,"error_Smoothing.txt");

    delete[] p_recover;
    delete[] pts2d_recover;
    delete[] vec_gamma;
    delete[] mot;
    delete[] pts2d;
    delete[] mot_temp;
    delete[] weight_one;
    delete[] mot_new;
    delete[] time_smoothing;
    delete[] error_smoothing;
    delete[] vmask;
}


