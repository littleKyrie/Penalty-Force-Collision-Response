#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(
        Eigen::Vector4d &phi,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::RowVectorXi> element,
        Eigen::Ref<const Eigen::Vector3d> X) {
    // edge matrix
    Eigen::Matrix3d T;

    // vertex X0
    Eigen::VectorXd X0 = V.row(element[0]);

    // compute values and assemble T
    for (size_t i = 0; i < T.rows(); i++)
        for (size_t j = 0; j < T.cols(); j++) {
            // the col index of T will represent the response vertex index of the element
            // the row index of T will represent the response component of vertex coordinate
            auto vertexIndex = static_cast<Eigen::Index>(j);
            auto coordinate = static_cast<Eigen::Index>(i);

            T(coordinate, vertexIndex) =
                    V(element[vertexIndex+1], coordinate) - X0[coordinate];
    }

    // compute the phi values based on T and input X in reference space
    Eigen::Vector3d deltaX;
    deltaX[0] = X[0] - X0[0];
    deltaX[1] = X[1] - X0[1];
    deltaX[2] = X[2] - X0[2];

    Eigen::Vector3d result;
    result = T.reverse() * deltaX;

    // get the linear phi
    phi[0] = 1 - result[0] - result[1] - result[2];
    phi[1] = result[0];
    phi[2] = result[1];
    phi[3] = result[2];
}