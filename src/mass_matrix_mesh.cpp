#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
#include <oneapi/tbb.h>
#include <thread>
#include <chrono>

void mass_matrix_mesh(
        Eigen::SparseMatrixd &M,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXi> T,
        double density, Eigen::Ref<const Eigen::VectorXd> v0) {
    // initialize
    M.resize(qdot.size(), qdot.size());
    M.setZero();

    // std::chrono::high_resolution_clock::time_point start, end;
    // start = std::chrono::high_resolution_clock::now();

    // concurrent
    std::mutex mtx;
    auto assemble_mass = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t j = r.begin(); j < r.end(); j++) {
            auto tetIndex = static_cast<Eigen::Index>(j);

            // compute Mj
            Eigen::Matrix1212d Mj;
            Eigen::RowVectorXi element = T.row(tetIndex);
            double volume = v0[tetIndex];
            mass_matrix_linear_tetrahedron(Mj, qdot, element, density, volume);

            // assemble M straightly
            for (size_t row = 0; row < element.size(); row++)
                for (size_t col = 0; col < element.size(); col++) {
                    auto rowGlobal = 3 * element[row];
                    auto colGlobal = 3 * element[col];

                    mtx.lock();
                    M.coeffRef(rowGlobal+0, colGlobal+0) += Mj(3*row+0, 3*col+0);
                    M.coeffRef(rowGlobal+1, colGlobal+1) += Mj(3*row+1, 3*col+1);
                    M.coeffRef(rowGlobal+2, colGlobal+2) += Mj(3*row+2, 3*col+2);
                    mtx.unlock();
                }
        }
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, T.rows()), assemble_mass);

    // end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> frameTime = end - start;
    // double frameRate = 1000.0 / frameTime.count();
    // std::cout << "the FPS of a concurrent mass is " << frameRate << std::endl;

    /*
    bool isSymmetricM = M.isApprox(M.transpose());
    if (!isSymmetricM) {
        // std::cout << "M is not symmetric" << std::endl;
        // return;
    } else {
        // std::cout << "M is symmetric" << std::endl;
    }
     */

    // start = std::chrono::high_resolution_clock::now();
    /*
    // serial
    for (size_t j = 0; j < T.rows(); j++) {
        auto tetIndex = static_cast<Eigen::Index>(j);

        // compute Mj
        Eigen::Matrix1212d Mj;
        Eigen::RowVectorXi element = T.row(tetIndex);
        double volume = v0[tetIndex];
        mass_matrix_linear_tetrahedron(Mj, qdot, element, density, volume);

        // assemble M straightly
        for (size_t row = 0; row < element.size(); row++)
            for (size_t col = 0; col < element.size(); col++) {
                auto rowGlobal = 3 * element[row];
                auto colGlobal = 3 * element[col];

                M.coeffRef(rowGlobal+0, colGlobal+0) += Mj(3*row+0, 3*col+0);
                M.coeffRef(rowGlobal+1, colGlobal+1) += Mj(3*row+1, 3*col+1);
                M.coeffRef(rowGlobal+2, colGlobal+2) += Mj(3*row+2, 3*col+2);
            }
    }
     */
    // end = std::chrono::high_resolution_clock::now();
    // frameTime = end - start;
    // frameRate = 1000.0 / frameTime.count();
    // std::cout << "the FPS of a serial mass is " << frameRate << std::endl;

    /*
    // test concurrent
    bool isSame = mass.isApprox(M.transpose());
    if (!isSame) {
        std::cout << "M is not equal to mass" << std::endl;
        // return;
        // exit(0);
    } else {
        std::cout << "mass can concurrent" << std::endl;
        // exit(0);
    }
     */

}

void massMatrixMultiBodies(
        Eigen::SparseMatrixd &M,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& pointers, double density) {
    M.resize(qdot.size(), qdot.size());
    M.setZero();

    for (int i = 0; i < objects.size(); i++) {
        int initialIndex = pointers[i];
        const SoftBody& obj = objects[i];
        Eigen::Ref<const Eigen::MatrixXi> T = obj.T;
        Eigen::Ref<const Eigen::VectorXd> v = obj.volumes;
        for (int j = 0; j < T.rows(); j++) {
            // compute Mj
            Eigen::Matrix1212d Mj;
            Eigen::RowVectorXi element = T.row(j);
            double volume = v[j];
            mass_matrix_linear_tetrahedron(Mj, qdot, element, density, volume);

            // assemble M straightly
            for (size_t row = 0; row < element.size(); row++)
                for (size_t col = 0; col < element.size(); col++) {
                    auto rowGlobal = initialIndex + 3 * element[row];
                    auto colGlobal = initialIndex + 3 * element[col];

                    M.coeffRef(rowGlobal+0, colGlobal+0) += Mj(3*row+0, 3*col+0);
                    M.coeffRef(rowGlobal+1, colGlobal+1) += Mj(3*row+1, 3*col+1);
                    M.coeffRef(rowGlobal+2, colGlobal+2) += Mj(3*row+2, 3*col+2);
                }
        }
    }

    M.makeCompressed();
}