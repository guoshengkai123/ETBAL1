//

//

#include "BA.h"
template<>
bool ReprojectionError<6>::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{
    util::point2d p;
    CameraParameters cam;
    cam.s=parameters[0][0];
    cam.alpha=parameters[0][1];
    cam.beta=parameters[0][2];
    cam.gamma=parameters[0][3];
    cam.t0=parameters[0][4];
    cam.t1=parameters[0][5];

    double t0 = cam.t0;
    double t1 = cam.t1;

    Eigen::Map<const Eigen::Vector3d> point(parameters[1]);

    double X, Y, Z;
    X=parameters[1][0];
    Y=parameters[1][1];
    Z=parameters[1][2];

    p=cam.Project(point);

    residuals[0] = p.x  - _observation_x;
    residuals[1] = p.y  - _observation_y;

    double cos_alpha = cos(cam.alpha);
    double sin_alpha = sin(cam.alpha);
    double cos_beta = cos(cam.beta);
    double sin_beta = sin(cam.beta);
    double cos_gamma = cos(cam.gamma);
    double sin_gamma = sin(cam.gamma);
    double s_s = 1/cam.s;

    if(jacobians !=NULL)
    {
        if(jacobians[0] != NULL)
        {
//            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor> > J_se3(jacobians[0]);
            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> J(jacobians[0]);
            J(0,0)=(-cos_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)+sin_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;;
            J(0,1) = (cos_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)-sin_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(0,2) = cos_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(0,3) = -sin_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-cos_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(0,4) = -cos_gamma;
            J(0,5) = sin_gamma;
            J(1,0) = (-sin_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)-cos_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
            J(1,1) = (sin_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)+cos_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(1,2) = sin_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(1,3) = cos_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-sin_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(1,4) = -sin_gamma;
            J(1,5) = -cos_gamma;
        }

        if(jacobians[1] != NULL)
        {
            jacobians[1][0] = cos_gamma*cos_beta*s_s;
            jacobians[1][1] = (cos_gamma*sin_alpha*sin_beta-sin_gamma*cos_alpha)*s_s;
            jacobians[1][2] = (-cos_gamma*cos_alpha*sin_beta-sin_gamma*sin_alpha)*s_s;
            jacobians[1][3] = sin_gamma*cos_beta*s_s;
            jacobians[1][4] = (sin_gamma*sin_alpha*sin_beta+cos_gamma*cos_alpha)*s_s;
            jacobians[1][5] = (-sin_gamma*cos_alpha*sin_beta+cos_gamma*sin_alpha)*s_s;
        }
    }

    return true;

}

util::point2d CameraParameters::Project(const Eigen::Vector3d& Point3D) {

    double cos_alpha;
    double sin_alpha;
    double cos_beta;
    double sin_beta;
    double cos_gamma;
    double sin_gamma;

    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    cos_beta = cos(beta);
    sin_beta = sin(beta);
    cos_gamma = cos(gamma);
    sin_gamma = sin(gamma);

    double a1,a2;
    double u,v;

    a1 = (cos_beta * Point3D[0] + sin_alpha * sin_beta * Point3D[1] - cos_alpha * sin_beta * Point3D[2]) / s - t0;
    a2 = (cos_alpha * Point3D[1] + sin_alpha * Point3D[2]) / s - t1;

    u = cos_gamma * a1 - sin_gamma * a2;
    v = sin_gamma * a1 + cos_gamma * a2;

    util::point2d p;
    p.x=u;
    p.y=v;

    return p;
}

void LM_BA(double *Camera_noise,double *imgpts_state,char* vmask,int ncams,int n3Dpts,double *newp){
    int gap,i;
    gap=0;
    ceres::Problem problem;
    double tmp=0;
    double x_o,y_o,dij;

    PosePointParametersBlock init_states;
    init_states.create(ncams, n3Dpts);

    for(i=0;i<ncams;i++)
    {
        init_states.values[i*6]=Camera_noise[i*6];
        init_states.values[i*6+1]=Camera_noise[i*6+1];
        init_states.values[i*6+2]=Camera_noise[i*6+2];
        init_states.values[i*6+3]=Camera_noise[i*6+3];
        init_states.values[i*6+4]=Camera_noise[i*6+4];
        init_states.values[i*6+5]=Camera_noise[i*6+5];
    }
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++)
    {
        init_states.values[i]=Camera_noise[i];
    }

    for(int i=0; i<n3Dpts; i++){
        for(int j=0;j<ncams; j++){
            dij=1;
            if(vmask[i*ncams+j]){
                x_o=imgpts_state[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts_state[i*(ncams*2)+2*j+1-gap*2];
                ceres::CostFunction* cost_function;
                cost_function=new ReprojectionError<6>(x_o, y_o);
                ceres::LossFunction* lossFunc = new MADN_loss(dij);
                problem.AddResidualBlock(cost_function,lossFunc, init_states.pose(j), init_states.point(i));
            }
            else{
                gap++;
            }
        }
    }

    ceres::Solver::Options options;
   
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = false;
    options.max_num_iterations = 1000;

    
    ceres::Solver::Summary summary;

//    ceres::Solver(options,&problem,&summary);
    ceres::Solve(options, &problem, &summary);
    for(i=0;i<ncams*6+n3Dpts*3;i++) newp[i] = init_states.values[i];
    //std::cout << summary.BriefReport() << "\n";

}