#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>
#include <iostream>
#include <chrono>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_q (temp_qdot) - scratch space for storing positions(or velocities)
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(
        Eigen::VectorXd &q,
        Eigen::VectorXd &qdot, double dt,
        const Eigen::SparseMatrixd &mass,
        ENERGY &energy, FORCE &force, STIFFNESS &stiffness,
        Eigen::VectorXd &tmp_q,
        Eigen::VectorXd &tmp_qdot,
        Eigen::VectorXd &tmp_force,
        const Eigen::SparseMatrixd &tmp_stiffness) {
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> frameTime = end - start;

    // set size
    tmp_force.resize(q.size());
    tmp_qdot.resize(qdot.size());
    tmp_q.resize(q.size());

    // for solving equation system
    Eigen::VectorXd tmp_g;
    Eigen::SparseMatrixd tmp_H;
    tmp_g.resize(q.size());
    tmp_H.resize(q.size(), q.size());

    // make initial guess
    // tmp_q = q;
    tmp_qdot = qdot;

    auto Jacobian = [&](Eigen::VectorXd& tmp_g, Eigen::VectorXd x) {
        tmp_force.setZero();
        force(tmp_force, q + dt * x, x);
        tmp_g = mass * (x - qdot) + (-tmp_force) * dt;
    };

    auto Hessian = [&](Eigen::SparseMatrixd& tmp_H, Eigen::VectorXd x) {
        std::chrono::high_resolution_clock::time_point lstart, lend;
        std::chrono::duration<double, std::milli> lTime = lend - lstart;

        lstart = std::chrono::high_resolution_clock::now();
        Eigen::SparseMatrixd K(tmp_stiffness);
        stiffness(K, q + dt * x, x);
        lend = std::chrono::high_resolution_clock::now();
        lTime = lend - lstart;
        // std::cout << "get full stiffness = " << lTime.count();

        lstart = std::chrono::high_resolution_clock::now();
        tmp_H = mass + K;
        lend = std::chrono::high_resolution_clock::now();
        lTime = lend - lstart;
        // std::cout << "stiffness -> Hessian = " << lTime.count();
    };

    // set iteration condition
    size_t maxSteps = 5;

    start = std::chrono::high_resolution_clock::now();
    // wrapper for implementing Newton method
    // double state = newtons_method(tmp_qdot, energy, Jacobian, Hessian, maxSteps, tmp_g, tmp_H);
    double state = projectedNewton(tmp_qdot, energy, Jacobian, Hessian, tmp_g, tmp_H);
    end = std::chrono::high_resolution_clock::now();
    frameTime = end - start;
    // std::cout << "the time of newton method = " << frameTime.count() << std::endl;

    // update state
    qdot = tmp_qdot;
    q = q + dt * qdot;
}

// for multi bodies simulation
template<typename ENERGY, typename FORCE, typename STIFFNESS, typename COLLISION>
inline void implicitEuler(
        Eigen::VectorXd &q,
        Eigen::VectorXd &qdot, double dt,
        const Eigen::SparseMatrixd &mass,
        ENERGY &energy, FORCE &force, STIFFNESS &stiffness, COLLISION &collision,
        Eigen::VectorXd &tmp_qdot,
        Eigen::VectorXd &tmp_force,
        const Eigen::SparseMatrixd &tmp_stiffness) {
    // set size
    tmp_force.resize(q.size());
    tmp_qdot.resize(q.size());

    // for solving equation system
    Eigen::VectorXd tmp_g;
    Eigen::SparseMatrixd tmp_H;
    tmp_g.resize(q.size());
    tmp_H.resize(q.size(), q.size());

    // make initial guess
    tmp_qdot = qdot;

    auto Jacobian = [&](Eigen::VectorXd& tmp_g, Eigen::VectorXd x) {
        tmp_force.setZero();
        force(tmp_force, q + dt * x, x);
        tmp_g = mass * (x - qdot) + (-tmp_force) * dt;
    };

    auto Hessian = [&](Eigen::SparseMatrixd& tmp_H, Eigen::VectorXd x) {
        Eigen::SparseMatrixd K(tmp_stiffness);
        stiffness(K, q + dt * x, x);
        tmp_H = mass + K;
    };

    // wrapper for implementing Newton method
    // double state = newtonsMethod(tmp_qdot, energy, Jacobian, Hessian, collision, maxSteps, tmp_g, tmp_H);
    double state = projectedNewtonMultiBodies(tmp_qdot, energy, Jacobian, Hessian, collision, tmp_g, tmp_H);

    // Eigen::VectorXd acc = (tmp_qdot - qdot) / dt;
    // std::cout << "acceleration: " << std::endl;
    // test::visualVector(acc);
    // std::cout << "total force: " << std::endl;
    // test::visualVector(mass * acc);

    // update state
    qdot = tmp_qdot;
    q = q + dt * qdot;
}