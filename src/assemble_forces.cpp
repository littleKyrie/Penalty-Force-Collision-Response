#include <assemble_forces.h>
#include <iostream>
#include <oneapi/tbb.h>
#include <thread>

void assemble_forces(
        Eigen::VectorXd &f,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        double C, double D) {
    // initialize inner force f
    f.resize(qdot.size());
    f.setZero();

    // concurrent without mutex
    std::vector<tbb::concurrent_vector<Eigen::Vector3d>> clusters;
    clusters.resize(V.rows());
    auto constructClusters = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t j = r.begin(); j < r.end(); j++) {
            auto tetIndex = static_cast<Eigen::Index>(j);

            Eigen::Vector12d dV;
            Eigen::Vector12d fj;
            Eigen::RowVectorXi element = T.row(tetIndex);
            double volume = v0(tetIndex);
            dV_linear_tetrahedron_dq(dV, q, V, element, volume, C, D);
            fj = -dV;

            for (size_t i = 0; i < element.size(); i++) {
                size_t vertexIndex = element[i];
                Eigen::Vector3d f = fj.segment(3*i, 3);
                clusters[vertexIndex].push_back(f);
            }
        }
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, T.rows()), constructClusters);

    auto assembleF = [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            for (const auto &ele : clusters[i]) {
                f.segment(3*i, 3) += ele;
            }
        }
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, clusters.size()), assembleF);

    // serial
    /*
    for (size_t j = 0; j < T.rows(); j++) {
        auto tetIndex = static_cast<Eigen::Index>(j);

        Eigen::Vector12d dV;
        Eigen::Vector12d fj;
        Eigen::RowVectorXi element = T.row(tetIndex);
        double volume = v0(tetIndex);

        dV_linear_tetrahedron_dq(dV, q, V, element, volume, C, D);
        fj = -dV;

        for (size_t row = 0; row < element.size(); row++) {
            auto rowIndex = element[row];
            f[rowIndex*3+0] += fj[3*row+0];
            f[rowIndex*3+1] += fj[3*row+1];
            f[rowIndex*3+2] += fj[3*row+2];
        }

    }
     */

    /*
    bool isSame = f.isApprox(foo);
    if (!isSame) {
        std::cout << "force can not be concurrent" << std::endl;
        // exit(0);
    } else {
        std::cout << "force can be concurrent" << std::endl;
        // exit(0);
    }
    */
    /*
    std::vector<size_t> nanEle_f;
    for (size_t a = 0; a < 12; a++){
        // check nan
        if (std::isnan(f[a])) {
            nanEle_f.push_back(a);
        }
    }
    std::cout << "the number of nan element in f " << "is " << nanEle_f.size() << std::endl;
    */
};

void assembleMultiForces(
        Eigen::VectorXd &f,
        Eigen::Ref<const Eigen::VectorXd> q,
        const CollidersPtr colliders,
        double C, double D) {
    f.resize(q.size());
    f.setZero();

    for (int id = 0; id < colliders->objects.size(); id++) {
        const SoftBody& obj = colliders->objects[id];
        int initialIndex = colliders->pointers[id];
        Eigen::Ref<const Eigen::MatrixXi> T = obj.T;
        Eigen::Ref<const Eigen::VectorXd> v = obj.volumes;
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;
        Eigen::Ref<const Eigen::VectorXd> x = q.segment(initialIndex, 3*V.rows());

        for (int j = 0; j < T.rows(); j++) {
            Eigen::Vector12d dV;
            Eigen::Vector12d fj;
            Eigen::RowVectorXi element = T.row(j);
            double volume = v(j);

            dV_linear_tetrahedron_dq(dV, x, V, element, volume, C, D);
            fj = -dV;

            for (int row = 0; row < element.size(); row++) {
                auto rowIndex = initialIndex + 3 * element[row];
                f[rowIndex+0] += fj[3*row+0];
                f[rowIndex+1] += fj[3*row+1];
                f[rowIndex+2] += fj[3*row+2];
            }
        }
    }
}

void assembleForcesFast(
        Eigen::VectorXd &f,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        std::vector<std::vector<Eigen::Vector3d>> &forceList,
        const Eigen::MatrixXi &tetPtr,
        double C, double D) {
    // initialize f for assembling force
    f.resize(q.size());
    f.setZero();

    // parallel initialize force list
    /*
   tbb::parallel_for_each(forceList.begin(), forceList.end(), [&](std::vector<Eigen::Vector3d> &elements) {
       for (auto &element : elements) {
           element.setZero();
       }
   });
     */

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> localTime = end - start;
    start = std::chrono::high_resolution_clock::now();
   // parallel compute force for a tetrahedron
   tbb::parallel_for(tbb::blocked_range<int>(0, T.rows()), [&](const tbb::blocked_range<int>& r) {
       for (int j = r.begin(); j < r.end(); j++) {
           // compute local forces
           Eigen::Vector12d dV;
           Eigen::Vector12d fj;
           Eigen::RowVectorXi element = T.row(j);
           double volume = v0(j);
           dV_linear_tetrahedron_dq(dV, q, V, element, volume, C, D);
           fj = -dV;

           // distribute local forces
           for (int i = 0; i < element.size(); i++) {
               int vi = element[i];
               int ids = tetPtr(j, i);
               forceList[vi][ids] = fj.segment(3*i, 3);
           }
       }
   });
   end = std::chrono::high_resolution_clock::now();
   localTime = end - start;
   // std::cout << "the time of concurrently computing forces = " << localTime.count() << std::endl;

   // assemble complete force
   start = std::chrono::high_resolution_clock::now();
   tbb::parallel_for(tbb::blocked_range<int>(0, V.rows()), [&](const tbb::blocked_range<int> &r) {
       for (int i = r.begin(); i < r.end(); i++) {
           for (const auto &force : forceList[i]) {
               f.segment(3*i, 3) += force;
           }
       }
   });
   end = std::chrono::high_resolution_clock::now();
   localTime = end - start;
   // std::cout << "the time of concurrently assembling forces = " << localTime.count() << std::endl;

   // serial code for test speed
   /*
   double comTime = 0;
   double disTime = 0;
   for (int j = 0; j < T.rows(); j++) {
       // compute local forces
       Eigen::Vector12d dV;
       Eigen::Vector12d fj;
       Eigen::RowVectorXi element = T.row(j);
       double volume = v0(j);
       start = std::chrono::high_resolution_clock::now();
       dV_linear_tetrahedron_dq(dV, q, V, element, volume, C, D);
       fj = -dV;
       end = std::chrono::high_resolution_clock::now();
       localTime = end - start;
       comTime += localTime.count();

       // distribute local forces
       start = std::chrono::high_resolution_clock::now();
       for (int i = 0; i < element.size(); i++) {
           int vi = element[i];
           int ids = tetPtr(j, i);
           forceList[vi][ids] = fj.segment(3*i, 3);
       }
       end = std::chrono::high_resolution_clock::now();
       localTime = end - start;
       disTime += localTime.count();
   }
   std::cout << "the average time of computing force of a tet = " << comTime / T.rows() << std::endl;
   std::cout << "the average time of distributing force = " << disTime / T.rows() << std::endl;
   */
}