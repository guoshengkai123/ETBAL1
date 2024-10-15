//

//

#ifndef BA_CERES_BA_H
#define BA_CERES_BA_H
#include "ceres/ceres.h"
#include <iostream>
#include <fstream> // file
#include "ceres/ceres.h"
#include "ceres/problem.h"

#include "matrix/matrix.h"
struct _point{
    union{
        float v[2];
        struct{float x;
            float y;};
    };
    _point(float _x, float _y):x(_x), y(_y){}
    _point(){}
    bool operator != (const _point& pt_){return !(x == pt_.x && y == pt_.y);}
};

struct _point3D{
    union{
        float v[3];
        struct{float x;
            float y;
            float z;};
    };
    _point3D(float _x, float _y, float _z):x(_x), y(_y), z(_z){}
    _point3D(){}
    bool operator != (const _point3D& pt_){return !(x == pt_.x && y == pt_.y && z == pt_.z);}
};

namespace util {
    typedef _point point2d;
    typedef _point3D point3d;
}

class CameraParameters
{
public:
    double s;
    double alpha;
    double beta;
    double gamma;
    double t0;
    double t1;
public:
    CameraParameters(){}
    CameraParameters(double s_, double alpha_, double beta_, double gamma_, double t0_, double t1_)
            : s(s_), alpha(alpha_), beta(beta_), gamma(gamma_), t0(t0_), t1(t1_){}
    util::point2d Project(const Eigen::Vector3d& Point3D);

};

class init_CameraParameters
{
public:
    std::vector<CameraParameters> cam;

public:
//    bool ReadCamera(const char* filename);
//    double* pose(int idx) {  return cam[idx] + idx * 6; }

};

class PosePointParametersBlock
{
public:
    PosePointParametersBlock(){}
    void create(int pose_num, int point_num)
    {
        poseNum = pose_num;
        pointNum = point_num;
        values = new double[pose_num * 6 + point_num * 3];
    }
    PosePointParametersBlock(int pose_num, int point_num): poseNum(pose_num), pointNum(point_num)
    {
        values = new double[pose_num * 6 + point_num * 3];
    }
//    ~PosePointParametersBlock() { delete[] values; }

    double* pose(int idx) {  return values + idx * 6; }

    double* point(int idx) { return values + poseNum * 6 + idx * 3; }


    int poseNum;
    int pointNum;
    double *values;    

};

template<int PoseBlockSize>
class ReprojectionError: public ceres::SizedCostFunction<2, PoseBlockSize, 3>
{
public:
    ReprojectionError(double observation_x,double observation_y):
            _observation_x(observation_x),
            _observation_y(observation_y){}

    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const;


private:
    double _observation_x;
    double _observation_y;

};

class Huber : public ceres::LossFunction {
public:
    explicit Huber(double delta) : delta_(delta) {}


    virtual void Evaluate(double s, double rho[3]) const {
//        std::cout<<"s:"<<s<<std::endl;
        const double delta2 = delta_ * delta_;
        if (s <= delta2) {
            rho[0] = s;         // loss
            rho[1] = 1.0;       // first derivative of loss
            rho[2] = 0.0;       // second derivative of loss
        } else {
            const double sqrt_s = sqrt(s);
            rho[0] = 2 * delta_ * sqrt_s - delta2;   // loss
            rho[1] = delta_ / sqrt_s;               // first derivative of loss
            rho[2] = -0.5 * delta_ / (s * sqrt_s);  // second derivative of loss
        }
    }
private:
    const double delta_;
};

class MADN_loss : public ceres::LossFunction {
public:
    explicit MADN_loss(double delta) : delta_(delta) {}

    virtual void Evaluate(double s, double rho[3]) const {
//        std::cout<<"delta_:"<<delta_<<std::endl;
        rho[0]=s*delta_;
        rho[1] = delta_;       // first derivative of loss
        rho[2] = 0.0;       // second derivative of loss
    }

private:
    const double delta_;
};
void LM_BA(double *Camera_noise,double *imgpts_state,char* vmask,int ncams,int n3Dpts,double *newp);
#endif //BA_CERES_BA_H
