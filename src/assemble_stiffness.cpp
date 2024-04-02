#include <assemble_stiffness.h>
#include <iostream>
#include <oneapi/tbb.h>
#include <thread>
#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <SparseStorage.h>

void assemble_stiffness(
        Eigen::SparseMatrixd &K,
        Eigen::Ref<const Eigen::MatrixXi> Map,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        double C, double D) {
    size_t rowSize = qdot.size();
    size_t colSize = qdot.size();
    // for speed test
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> localTime = end - start;

    /*
    // concurrent without mutex
    // start = std::chrono::high_resolution_clock::now();
    std::vector<tbb::concurrent_vector<Eigen::Matrix3d>> clusters;
    clusters.resize(V.rows()*V.rows()); // row first
    // end = std::chrono::high_resolution_clock::now();
    // localTime = end - start;
    // std::cout << "time of initializing clusters is " << localTime.count() << std::endl;

    auto collectClusters = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t j = r.begin(); j < r.end(); j++) {
            Eigen::Matrix1212d Hj;
            double volume = v0[j];
            Eigen::RowVectorXi element = T.row(j);
            d2V_linear_tetrahedron_dq2(Hj, q, V, element, volume, C, D);
            // Hj *= -1;

            // divide Hj into subblocks
            Eigen::Matrix3d Hj00 = Hj.block(0, 0, 3, 3);//3*3 matrix
            Eigen::Matrix3d Hj01 = Hj.block(0, 3, 3, 3);
            Eigen::Matrix3d Hj02 = Hj.block(0, 6, 3, 3);
            Eigen::Matrix3d Hj03 = Hj.block(0, 9, 3, 3);

            Eigen::Matrix3d Hj10 = Hj.block(3, 0, 3, 3);
            Eigen::Matrix3d Hj11 = Hj.block(3, 3, 3, 3);
            Eigen::Matrix3d Hj12 = Hj.block(3, 6, 3, 3);
            Eigen::Matrix3d Hj13 = Hj.block(3, 9, 3, 3);

            Eigen::Matrix3d Hj20 = Hj.block(6, 0, 3, 3);
            Eigen::Matrix3d Hj21 = Hj.block(6, 3, 3, 3);
            Eigen::Matrix3d Hj22 = Hj.block(6, 6, 3, 3);
            Eigen::Matrix3d Hj23 = Hj.block(6, 9, 3, 3);

            Eigen::Matrix3d Hj30 = Hj.block(9, 0, 3, 3);
            Eigen::Matrix3d Hj31 = Hj.block(9, 3, 3, 3);
            Eigen::Matrix3d Hj32 = Hj.block(9, 6, 3, 3);
            Eigen::Matrix3d Hj33 = Hj.block(9, 9, 3, 3);

            // collect relevant blocks
            clusters[element[0]*V.rows()+element[0]].push_back(Hj00);
            clusters[element[0]*V.rows()+element[1]].push_back(Hj01);
            clusters[element[0]*V.rows()+element[2]].push_back(Hj02);
            clusters[element[0]*V.rows()+element[3]].push_back(Hj03);

            clusters[element[1]*V.rows()+element[0]].push_back(Hj10);
            clusters[element[1]*V.rows()+element[1]].push_back(Hj11);
            clusters[element[1]*V.rows()+element[2]].push_back(Hj12);
            clusters[element[1]*V.rows()+element[3]].push_back(Hj13);

            clusters[element[2]*V.rows()+element[0]].push_back(Hj20);
            clusters[element[2]*V.rows()+element[1]].push_back(Hj21);
            clusters[element[2]*V.rows()+element[2]].push_back(Hj22);
            clusters[element[2]*V.rows()+element[3]].push_back(Hj23);

            clusters[element[3]*V.rows()+element[0]].push_back(Hj30);
            clusters[element[3]*V.rows()+element[1]].push_back(Hj31);
            clusters[element[3]*V.rows()+element[2]].push_back(Hj32);
            clusters[element[3]*V.rows()+element[3]].push_back(Hj33);
        }
    };
    // start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, T.rows()), collectClusters);
    // end = std::chrono::high_resolution_clock::now();
    // localTime = end - start;
    // std::cout << "time of collecting clusters is " << localTime.count() << std::endl;

    auto fastFillStiffness = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            int rowGlobal = i / V.rows();
            int colGlobal = i % V.rows();

            for (const auto &block : clusters[i]) {
                for (size_t jj = 0; jj < 3; jj++) {
                    for (Eigen::SparseMatrixd::InnerIterator it(K, 3 * colGlobal + jj); it; ++it) {
                        if (it.row() == 3 * rowGlobal) {
                            it.valueRef() = it.value() + block.coeff(0, jj);
                            ++it;
                            it.valueRef() = it.value() + block.coeff(1, jj);
                            ++it;
                            it.valueRef() = it.value() + block.coeff(2, jj);
                            break;
                        }
                    }
                }
            }
        }
    };

    // start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, clusters.size()), fastFillStiffness);
    // end = std::chrono::high_resolution_clock::now();
    // localTime = end - start;
    // std::cout << "time of writing value concurrently is " << localTime.count() << std::endl;

    auto visitClusters = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            int rowGlobal = i / V.rows();
            int colGlobal = i % V.rows();

            for (const auto &block : clusters[i]) {
                for (size_t jj = 0; jj < 3; jj++) {
                    for (Eigen::SparseMatrixd::InnerIterator it(K, 3 * colGlobal + jj); it; ++it) {
                        if (it.row() == 3 * rowGlobal) {
                            it.value() + block.coeff(0, jj);
                            ++it;
                            it.value() + block.coeff(1, jj);
                            ++it;
                            it.value() + block.coeff(2, jj);
                            break;
                        }
                    }
                }
            }
        }
    };
    // start = std::chrono::high_resolution_clock::now();
    // tbb::parallel_for(tbb::blocked_range<size_t>(0, clusters.size()), visitClusters);
    // end = std::chrono::high_resolution_clock::now();
    // localTime = end - start;
    // std::cout << "time of visiting clusters is " << localTime.count() << std::endl;
    */

    // serial with triplet
    Eigen::SparseMatrixd E;
    int r = qdot.rows(); // r = 3n
    E.setZero();
    E.resize(r, r);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(r * r);

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < T.rows(); i++) {

        Eigen::Matrix1212d H_i;

        std::chrono::high_resolution_clock::time_point begin, finish;
        begin = std::chrono::high_resolution_clock::now();
        d2V_linear_tetrahedron_dq2(H_i, q, V, T.row(i), v0(i), C, D);
        // H_i *= -1;
        finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> dt = finish - begin;
        // std::cout << "the time for Hj " << dt.count() << std::endl;

        Eigen::Matrix3d H_i00 = H_i.block(0, 0, 3, 3);//3*3 matrix
        Eigen::Matrix3d H_i01 = H_i.block(0, 3, 3, 3);
        Eigen::Matrix3d H_i02 = H_i.block(0, 6, 3, 3);
        Eigen::Matrix3d H_i03 = H_i.block(0, 9, 3, 3);

        Eigen::Matrix3d H_i10 = H_i.block(3, 0, 3, 3);
        Eigen::Matrix3d H_i11 = H_i.block(3, 3, 3, 3);
        Eigen::Matrix3d H_i12 = H_i.block(3, 6, 3, 3);
        Eigen::Matrix3d H_i13 = H_i.block(3, 9, 3, 3);

        Eigen::Matrix3d H_i20 = H_i.block(6, 0, 3, 3);
        Eigen::Matrix3d H_i21 = H_i.block(6, 3, 3, 3);
        Eigen::Matrix3d H_i22 = H_i.block(6, 6, 3, 3);
        Eigen::Matrix3d H_i23 = H_i.block(6, 9, 3, 3);

        Eigen::Matrix3d H_i30 = H_i.block(9, 0, 3, 3);
        Eigen::Matrix3d H_i31 = H_i.block(9, 3, 3, 3);
        Eigen::Matrix3d H_i32 = H_i.block(9, 6, 3, 3);
        Eigen::Matrix3d H_i33 = H_i.block(9, 9, 3, 3);

        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                tripleList.push_back(Trip(3 * T(i, 0) + ii, 3 * T(i, 0) + jj, H_i00.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 0) + ii, 3 * T(i, 1) + jj, H_i01.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 0) + ii, 3 * T(i, 2) + jj, H_i02.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 0) + ii, 3 * T(i, 3) + jj, H_i03.coeff(ii, jj)));


                tripleList.push_back(Trip(3 * T(i, 1) + ii, 3 * T(i, 0) + jj, H_i10.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 1) + ii, 3 * T(i, 1) + jj, H_i11.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 1) + ii, 3 * T(i, 2) + jj, H_i12.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 1) + ii, 3 * T(i, 3) + jj, H_i13.coeff(ii, jj)));

                tripleList.push_back(Trip(3 * T(i, 2) + ii, 3 * T(i, 0) + jj, H_i20.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 2) + ii, 3 * T(i, 1) + jj, H_i21.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 2) + ii, 3 * T(i, 2) + jj, H_i22.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 2) + ii, 3 * T(i, 3) + jj, H_i23.coeff(ii, jj)));

                tripleList.push_back(Trip(3 * T(i, 3) + ii, 3 * T(i, 0) + jj, H_i30.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 3) + ii, 3 * T(i, 1) + jj, H_i31.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 3) + ii, 3 * T(i, 2) + jj, H_i32.coeff(ii, jj)));
                tripleList.push_back(Trip(3 * T(i, 3) + ii, 3 * T(i, 3) + jj, H_i33.coeff(ii, jj)));

            }
        }

    }
    end = std::chrono::high_resolution_clock::now();
    localTime = end - start;
    std::cout << "the time of obtaining tripletList is  " << localTime.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    // E.setFromTriplets(tripleList.begin(), tripleList.end());
    K.setFromTriplets(tripleList.begin(), tripleList.end());
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> frameTime = end - start;
    std::cout << "the time of a serial K is  " << frameTime.count() << std::endl;
    // double frameRate = 1000.0 / frameTime.count();
    // std::cout << "the FPS of a serial K is " << frameRate << std::endl;

    /*
    // test concurrent
    bool isSame = K.isApprox(E);
    if (!isSame) {
        std::cout << "E is not equal to K" << std::endl;
        // return;
        // exit(0);
    } else {
        std::cout << "stiffness can concurrent" << std::endl;
        // exit(0);
    }
     */
};

void assembleMultiStiffness(
        Eigen::SparseMatrixd& K,
        Eigen::Ref<const Eigen::VectorXd> q,
        const CollidersPtr colliders,
        double C, double D) {
    int rowSize = q.size();
    int colSize = q.size();
    // for serial
    K.resize(rowSize, colSize);
    K.setZero();

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;

    for (int id = 0; id < colliders->objects.size(); id++) {
        const SoftBody& obj = colliders->objects[id];
        int initialIndex = colliders->pointers[id];
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;
        Eigen::Ref<const Eigen::MatrixXi> T = obj.T;
        Eigen::Ref<const Eigen::VectorXd> v = obj.volumes;
        Eigen::Ref<const Eigen::VectorXd> x = q.segment(initialIndex, 3 * V.rows());

        for (int i = 0; i < T.rows(); i++) {
            Eigen::Matrix1212d H_i;

            d2V_linear_tetrahedron_dq2(H_i, x, V, T.row(i), v(i), C, D);

            Eigen::Matrix3d H_i00 = H_i.block(0, 0, 3, 3);//3*3 matrix
            Eigen::Matrix3d H_i01 = H_i.block(0, 3, 3, 3);
            Eigen::Matrix3d H_i02 = H_i.block(0, 6, 3, 3);
            Eigen::Matrix3d H_i03 = H_i.block(0, 9, 3, 3);

            Eigen::Matrix3d H_i10 = H_i.block(3, 0, 3, 3);
            Eigen::Matrix3d H_i11 = H_i.block(3, 3, 3, 3);
            Eigen::Matrix3d H_i12 = H_i.block(3, 6, 3, 3);
            Eigen::Matrix3d H_i13 = H_i.block(3, 9, 3, 3);

            Eigen::Matrix3d H_i20 = H_i.block(6, 0, 3, 3);
            Eigen::Matrix3d H_i21 = H_i.block(6, 3, 3, 3);
            Eigen::Matrix3d H_i22 = H_i.block(6, 6, 3, 3);
            Eigen::Matrix3d H_i23 = H_i.block(6, 9, 3, 3);

            Eigen::Matrix3d H_i30 = H_i.block(9, 0, 3, 3);
            Eigen::Matrix3d H_i31 = H_i.block(9, 3, 3, 3);
            Eigen::Matrix3d H_i32 = H_i.block(9, 6, 3, 3);
            Eigen::Matrix3d H_i33 = H_i.block(9, 9, 3, 3);

            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 0) + ii, initialIndex + 3 * T(i, 0) + jj, H_i00.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 0) + ii, initialIndex + 3 * T(i, 1) + jj, H_i01.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 0) + ii, initialIndex + 3 * T(i, 2) + jj, H_i02.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 0) + ii, initialIndex + 3 * T(i, 3) + jj, H_i03.coeff(ii, jj)));


                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 1) + ii, initialIndex + 3 * T(i, 0) + jj, H_i10.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 1) + ii, initialIndex + 3 * T(i, 1) + jj, H_i11.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 1) + ii, initialIndex + 3 * T(i, 2) + jj, H_i12.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 1) + ii, initialIndex + 3 * T(i, 3) + jj, H_i13.coeff(ii, jj)));

                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 2) + ii, initialIndex + 3 * T(i, 0) + jj, H_i20.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 2) + ii, initialIndex + 3 * T(i, 1) + jj, H_i21.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 2) + ii, initialIndex + 3 * T(i, 2) + jj, H_i22.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 2) + ii, initialIndex + 3 * T(i, 3) + jj, H_i23.coeff(ii, jj)));

                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 3) + ii, initialIndex + 3 * T(i, 0) + jj, H_i30.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 3) + ii, initialIndex + 3 * T(i, 1) + jj, H_i31.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 3) + ii, initialIndex + 3 * T(i, 2) + jj, H_i32.coeff(ii, jj)));
                    tripleList.push_back(Trip(initialIndex + 3 * T(i, 3) + ii, initialIndex + 3 * T(i, 3) + jj, H_i33.coeff(ii, jj)));

                }
            }
        }
    }
    K.setFromTriplets(tripleList.begin(), tripleList.end());
    K.makeCompressed();
}

void assembleHessianFast(
        Eigen::SparseMatrixd &K,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        std::vector<std::map<int, Eigen::Matrix3d>> &blockList,
        double C, double D) {
    int rowSize = V.rows();
    int colSize = V.rows();

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> localTime = end - start;

    start = std::chrono::high_resolution_clock::now();
    // parallel compute a Hessian for a tetrahedron
    tbb::parallel_for(tbb::blocked_range<int>(0, T.rows()), [&](const tbb::blocked_range<int> &r) {
        for (int j = r.begin(); j < r.end(); j++) {
            Eigen::Matrix1212d Hj;
            double volume = v0[j];
            Eigen::RowVectorXi element = T.row(j);
            d2V_linear_tetrahedron_dq2(Hj, q, V, element, volume, C, D);

            // divide Hj into subblocks
            Eigen::Matrix3d Hj00 = Hj.block(0, 0, 3, 3);//3*3 matrix
            Eigen::Matrix3d Hj01 = Hj.block(0, 3, 3, 3);
            Eigen::Matrix3d Hj02 = Hj.block(0, 6, 3, 3);
            Eigen::Matrix3d Hj03 = Hj.block(0, 9, 3, 3);

            Eigen::Matrix3d Hj10 = Hj.block(3, 0, 3, 3);
            Eigen::Matrix3d Hj11 = Hj.block(3, 3, 3, 3);
            Eigen::Matrix3d Hj12 = Hj.block(3, 6, 3, 3);
            Eigen::Matrix3d Hj13 = Hj.block(3, 9, 3, 3);

            Eigen::Matrix3d Hj20 = Hj.block(6, 0, 3, 3);
            Eigen::Matrix3d Hj21 = Hj.block(6, 3, 3, 3);
            Eigen::Matrix3d Hj22 = Hj.block(6, 6, 3, 3);
            Eigen::Matrix3d Hj23 = Hj.block(6, 9, 3, 3);

            Eigen::Matrix3d Hj30 = Hj.block(9, 0, 3, 3);
            Eigen::Matrix3d Hj31 = Hj.block(9, 3, 3, 3);
            Eigen::Matrix3d Hj32 = Hj.block(9, 6, 3, 3);
            Eigen::Matrix3d Hj33 = Hj.block(9, 9, 3, 3);

            // preserve these blocks
            blockList[element[0]*rowSize+element[0]].at(j) = Hj00;
            blockList[element[0]*rowSize+element[1]].at(j) = Hj01;
            blockList[element[0]*rowSize+element[2]].at(j) = Hj02;
            blockList[element[0]*rowSize+element[3]].at(j) = Hj03;

            blockList[element[1]*rowSize+element[0]].at(j) = Hj10;
            blockList[element[1]*rowSize+element[1]].at(j) = Hj11;
            blockList[element[1]*rowSize+element[2]].at(j) = Hj12;
            blockList[element[1]*rowSize+element[3]].at(j) = Hj13;

            blockList[element[2]*rowSize+element[0]].at(j) = Hj20;
            blockList[element[2]*rowSize+element[1]].at(j) = Hj21;
            blockList[element[2]*rowSize+element[2]].at(j) = Hj22;
            blockList[element[2]*rowSize+element[3]].at(j) = Hj23;

            blockList[element[3]*rowSize+element[0]].at(j) = Hj30;
            blockList[element[3]*rowSize+element[1]].at(j) = Hj31;
            blockList[element[3]*rowSize+element[2]].at(j) = Hj32;
            blockList[element[3]*rowSize+element[3]].at(j) = Hj33;
        }
    });
    end = std::chrono::high_resolution_clock::now();
    localTime = end - start;
    // std::cout << "the time of concurrently computing the Hessian of all tetrahedrons = " << localTime.count() << std::endl;

    // assemble complete Hessian
    start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<int>(0, rowSize*colSize), [&](const tbb::blocked_range<int> &r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            int rowGlobal = i / V.rows();
            int colGlobal = i % V.rows();

            for (const auto &block : blockList[i]) {
                for (size_t jj = 0; jj < 3; jj++) {
                    for (Eigen::SparseMatrixd::InnerIterator it(K, 3 * colGlobal + jj); it; ++it) {
                        if (it.row() == 3 * rowGlobal) {
                            it.valueRef() = it.value() + block.second(0, jj);
                            ++it;
                            it.valueRef() = it.value() + block.second(1, jj);
                            ++it;
                            it.valueRef() = it.value() + block.second(2, jj);
                            break;
                        }
                    }
                }
            }
        }
    });
    end = std::chrono::high_resolution_clock::now();
    localTime = end - start;
    // std::cout << "the time of concurrently assembling all Hessians to the system Hessian = " << localTime.count() << std::endl;
}