#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

void T_linear_tetrahedron(
        double &Tj,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::RowVectorXi> element,
        double density, double volume) {

    // construct selective matrix
    Eigen::MatrixXd Ej(12, qdot.size());
    for (size_t i = 0; i < element.size(); i++) {
        auto vertexIndex = static_cast<Eigen::Index>(i);
        Ej(vertexIndex*3+0, element[vertexIndex]*3+0) = 1;
        Ej(vertexIndex*3+1, element[vertexIndex]*3+1) = 1;
        Ej(vertexIndex*3+2, element[vertexIndex]*3+2) = 1;
    }
    Eigen::VectorXd qdot_j = Ej * qdot;


    // get Tj
    Eigen::Matrix1212d Mj;
    Eigen::Vector12d v;
    for (size_t i = 0; i < element.size(); i++) {
        v.segment(3*i, 3) = qdot.segment(3*element[i], 3);
    }
    mass_matrix_linear_tetrahedron(Mj, qdot, element, density, volume);
    Tj = 0.5 * v.dot(Mj * v);
    double T = 0.5 * v.transpose() * Mj * v;

    /*
    bool isSame = v.isApprox(qdot_j);
    if (!isSame) {
        std::cout << "different selective result" << std::endl;
        std::cout << "qj size = " << qdot_j.size() << std::endl;
        // exit(0);
    } else {
        std::cout << "correct selective" << std::endl;
    }
     */

    /*
    double d = abs(abs(Tj) - abs(T));
    if (d > 0.001) {
        std::cout << "error T" << std::endl;
        exit(0);
    } else {
        std::cout << "correct Ti is " << Tj << "and T is " << T << std::endl;
    }
     */
}