#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>

void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {

    std::chrono::high_resolution_clock::time_point start, end;

    psi = 0.0;
    // svd decomposition for deformation gradient F
    start = std::chrono::high_resolution_clock::now();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> frameTime = end - start;
    // std::cout << "svd for F is " << frameTime.count() << std::endl;

    // get the singular values
    start = std::chrono::high_resolution_clock::now();
    Eigen::Vector3d singularValues = svd.singularValues();
    // compute energy density
    double I1 = singularValues.squaredNorm();
    double I2 = singularValues.array().pow(4).sum();
    double I3 = singularValues.array().pow(2).prod();
    end = std::chrono::high_resolution_clock::now();
    frameTime = end - start;
    // std::cout << "svd computing is " << frameTime.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    double log_I3 = std::log(I3);
    double log_I3_squared = std::pow(log_I3, 2);
    psi = 0.5 * C * (I1 - log_I3 - 3) + D / 8.0 * log_I3_squared;
    end = std::chrono::high_resolution_clock::now();
    frameTime = end - start;
    // std::cout << "final computing is " << frameTime.count() << std::endl;

    /*
    //------------------------
    double e;
    double J = F.determinant();
    double trace = (F.transpose() * F).trace();
    e = C * ((pow(J, -2.0 / 3.0) * trace) - 3) + D * pow(J - 1, 2);
    */
}