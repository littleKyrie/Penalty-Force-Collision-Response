#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>

#include <preconditioner.h>
#include <Jacobi_solver.h>
#include <conjugate_gradient_solver.h>
#include <direct_solver.h>
#include <Gauss_Seidel_Solver.h>
#include <colliders.h>
#include <spatial_hash.h>
#include <testTools.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(
        Eigen::VectorXd &x0,
        Objective &f,
        Jacobian &g,
        Hessian &H,
        unsigned int maxSteps,
        Eigen::VectorXd &tmp_g,
        Eigen::SparseMatrixd &tmp_H) {
    Eigen::VectorXd x = x0;

    // config parameters
    double alpha = 1.0;
    double p = 0.5;
    double c = 1e-6; // tolerance

    // set termination condition
    size_t i = 0;
    double gradNorm = 0;
    g(tmp_g, x);
    gradNorm = tmp_g.norm();
    double originalGradNorm = gradNorm;

    // iteration
    while (gradNorm > c * originalGradNorm) {
        std::cout << "iteration = " << i << "grad norm = " << gradNorm << std::endl;
        std::cout << "relative grad norm = " << gradNorm / originalGradNorm << std::endl;
        double solutionNorm = x.norm();

        // construct linear equations system
        H(tmp_H, x);

        // use Eigen solver
        // Eigen::VectorXd dir;
        // Eigen CG
        // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        // solver.compute(tmp_H);
        // dir = solver.solve(-1.0 * tmp_g);
        // Eigen direct solver
        // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        // solver.compute(tmp_H);
        // dir = solver.solve(-tmp_g);

        // use my solver
        Eigen::VectorXd dir;
        Eigen::VectorXd b = -1 * tmp_g;
        // Eigen::SparseMatrixd M;
        // bool isInverseM = jacobiPreconditioning(M, tmp_H);
        // bool isInverseM = choleskyPreconditioning(M, tmp_H);
        // bool solved = jacobiSolver(dir, tmp_H, b, true);
        // jacobi(dir, tmp_H, b);
        // gaussSeidel(dir, tmp_H, b, true);
        // bool solved = jacobi(dir, tmp_H, b, false);
        // bool solved = conjugateGradient(dir, tmp_H, b);
        // bool solved = preconditionedConjugateGradient(dir, tmp_H, b, M, isInverseM);
        bool solved = PCG(dir, tmp_H, b, "Cholesky");

        // backtracking line search
        alpha = 1.0;
        double energyNew = f(x + alpha * dir);
        // double energyUpRange = f(x) + c * alpha * dir.transpose() * tmp_g;
        double energyUpRange = f(x);
        while (energyNew > energyUpRange && alpha > c) {
            alpha *= p;
            energyNew = f(x + alpha * dir);
        }
        x = x + alpha * dir;

        double differenceNorm = alpha * dir.norm();
        std:: cout << "the relative difference of solution = " << differenceNorm / solutionNorm << std::endl;

        if ((differenceNorm / solutionNorm) < c) {
            break;
        }

        // update iteration condition
        i++;
        g(tmp_g, x);
        gradNorm = tmp_g.norm();
    }

    // start = std::chrono::high_resolution_clock::now();
    x0 = x;
    // end = std::chrono::high_resolution_clock::now();
    // frameTime = end - start;
    // std::cout << "the time of write back vectorXd " << frameTime.count() << std::endl;

    // std::cout << "iteration of Newton is " << i << std::endl;
    // std::cout << "the final gradient norm is " << gradNorm << std::endl;
   
    return 0.0;
}

// for multi bodies simulation
template<typename Objective, typename Jacobian, typename Hessian, typename Collision>
double newtonsMethod(
        Eigen::VectorXd &x0,
        Objective &f,
        Jacobian &g,
        Hessian &H,
        Collision &collision,
        unsigned int maxSteps,
        Eigen::VectorXd &tmp_g,
        Eigen::SparseMatrixd &tmp_H) {
    Eigen::VectorXd x = x0;

    // config parameters
    double alpha = 1.0;
    double p = 0.5;
    double c = 1e-8; // tolerance

    // set termination condition
    size_t i = 0;
    double gradNorm = 0;

    // iteration
    while (i < maxSteps) {
        // collision(x);

        g(tmp_g, x);
        gradNorm = tmp_g.norm();

        if (gradNorm <= c) {
            break;
        }

        // construct linear equations system
        H(tmp_H, x);

        // use Eigen solver
        Eigen::VectorXd dir;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.compute(tmp_H);
        dir = solver.solve(-1.0 * tmp_g);

        // backtracking line search
        alpha = 1.0;
        // std::cout << c * dir.transpose() * tmp_g << std::endl;
        // std::cout << "dir: " << std::endl;
        // test::visualVector(dir);
        double energyUpRange = f(x) + c * dir.transpose() * tmp_g;
        // collision(x + alpha * dir);
        double energyNew = f(x + alpha * dir);
        while (energyNew > energyUpRange && alpha > c) {
            alpha *= p;
            // collision(x + alpha * dir);
            energyNew = f(x + alpha * dir);
        }
        x = x + alpha * dir;

        // std::cout << "X at iteration " << i+1 << " with alpha = " << alpha << std::endl;
        // test::visualVector(x);

        // update iteration condition
        i++;
    }

    x0 = x;

    return 0.0;
}

// test another kind of convergence condition mentioned in <IPC> based on infinity norm
template<typename Objective, typename Jacobian, typename Hessian>
double projectedNewton(
        Eigen::VectorXd &x0,
        Objective &f,
        Jacobian &g,
        Hessian &H,
        Eigen::VectorXd &tmp_g,
        Eigen::SparseMatrixd &tmp_H) {
    std::chrono::high_resolution_clock::time_point outStart, outEnd;
    std::chrono::duration<double, std::milli> outTime = outEnd - outStart;
    // initialize medial variable
    Eigen::VectorXd x = x0;

    // config parameters
    double alpha = 1.0;
    double p = 0.5;
    double c = 1e-3;
    size_t i = 0;

    // set absolute convergence condition
    double infinityNorm = 0;

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> Time = end - start;
    // iteration
    do {
        outStart = std::chrono::high_resolution_clock::now();
        std::cout << "iteration = " << i << std::endl;

        // prepare for relative convergence condition
        double lastX = x.norm();

        // compute gradient at xi
        start = std::chrono::high_resolution_clock::now();
        g(tmp_g, x);
        end = std::chrono::high_resolution_clock::now();
        Time = end - start;
        // std::cout << "the time of getting gradient = " << Time.count() << std::endl;

        // compute Hessian at xi
        start = std::chrono::high_resolution_clock::now();
        H(tmp_H, x);
        end = std::chrono::high_resolution_clock::now();
        Time = end - start;
        // std::cout << "the time of getting Hessian = " << Time.count() << std::endl;

        // use my solver
        Eigen::VectorXd dir;
        Eigen::VectorXd b = -1 * tmp_g;
        start = std::chrono::high_resolution_clock::now();
        bool solved = PCG(dir, tmp_H, b, "Cholesky");
        // bool solved = gaussSeidel(dir, tmp_H, b, true);
        // bool solved = jacobi(dir, tmp_H, b, true);
        end = std::chrono::high_resolution_clock::now();
        Time = end - start;
        // std::cout << "the time of PCG solver = " << Time.count() << std::endl;

        // backtracking line search
        start = std::chrono::high_resolution_clock::now();
        alpha = 1.0;
        double energyNew = f(x + alpha * dir);
        // double energyUpRange = f(x) + c * alpha * dir.transpose() * tmp_g;
        double energyUpRange = f(x);
        while (energyNew > energyUpRange && alpha > c) {
            alpha *= p;
            energyNew = f(x + alpha * dir);
        }
        end = std::chrono::high_resolution_clock::now();
        Time = end - start;
        // std::cout << "time of line search = " << Time.count() << std::endl;
        x = x + alpha * dir;

        // update iteration condition
        i++;

        // check for convergence
        // by relative convergence conditions
        double relativeNorm = alpha * dir.norm() / lastX;
        if (relativeNorm < c) {
            // std::cout << "the relative norm of solutions = " << relativeNorm << std::endl;
            break;
        }

        // by absolute convergence condition
        infinityNorm = dir.lpNorm<Eigen::Infinity>();
        // std::cout << "the infinity norm of direction = " << infinityNorm << std::endl;

        outEnd = std::chrono::high_resolution_clock::now();
        outTime = outEnd - outStart;
        // std::cout << "time of a single newton = " << outTime.count() << std::endl;

    } while (alpha * infinityNorm > c);

    // write back values
    x0 = x;

    return 0.;
}

// check convergence of newton method for multi-bodies by both relative and absolute conditions
template<typename Objective, typename Jacobian, typename Hessian, typename Collision>
double projectedNewtonMultiBodies(
        Eigen::VectorXd &x0,
        Objective &f,
        Jacobian &g,
        Hessian &H,
        Collision &C,
        Eigen::VectorXd &tmp_g,
        Eigen::SparseMatrixd &tmp_H) {
    // initialize medial variable
    Eigen::VectorXd x = x0;

    // config parameters
    double alpha = 1.0;
    double p = 0.5;
    double tol = 1e-3;
    size_t i = 0;

    // set absolute convergence condition
    double infinityNorm = 0;

    // iteration
    do {
        std::cout << "iteration = " << i << std::endl;

        // prepare for relative convergence condition
        double lastX = x.norm();

        // initialize contact constraints
        C(x);

        // compute gradient at xi
        g(tmp_g, x);

        // compute Hessian at xi
        H(tmp_H, x);

        // use my solver
        Eigen::VectorXd dir;
        Eigen::VectorXd b = -1 * tmp_g;
        bool solved = PCG(dir, tmp_H, b, "Cholesky");
        // bool solved = gaussSeidel(dir, tmp_H, b, true);
        // bool solved = jacobi(dir, tmp_H, b, true);

        // backtracking line search
        alpha = 1.0;
        // double energyUpRange = f(x) + c * alpha * dir.transpose() * tmp_g;
        double energyUpRange = f(x);
        // C(x + alpha * dir);
        double energyNew = f(x + alpha * dir);
        while (energyNew > energyUpRange && alpha > tol) {
            alpha *= p;
            // C(x + alpha * dir);
            energyNew = f(x + alpha * dir);
        }

        x = x + alpha * dir;

        // update iteration condition
        i++;

        // check for convergence
        // by relative convergence conditions
        double relativeNorm = alpha * dir.norm() / lastX;
        if (relativeNorm < tol) {
            // std::cout << "the relative norm of solutions = " << relativeNorm << std::endl;
            break;
        }

        // by absolute convergence condition
        infinityNorm = dir.lpNorm<Eigen::Infinity>();
        // std::cout << "the infinity norm of direction = " << infinityNorm << std::endl;

    } while (alpha * infinityNorm > tol);

    // write back values
    x0 = x;

    return 0.;
}