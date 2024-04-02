#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(
        Ret &&integrated,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
        Integrand_Func integrand) {
    // compute an integrated value on X(such as the centroid on the tetrahedron)
    // the integrated function can be energy,strain,stress or something else
    // right value reference
    Eigen::Vector3d X;
    X[0] = (q[element[0]*3+0] + q[element[1]*3+0] + q[element[2]*3+0] + q[element[3]*3+0]) / 4.0;
    X[1] = (q[element[0]*3+1] + q[element[1]*3+1] + q[element[2]*3+1] + q[element[3]*3+1]) / 4.0;
    X[2] = (q[element[0]*3+2] + q[element[1]*3+2] + q[element[2]*3+2] + q[element[3]*3+2]) / 4.0;

    integrand(integrated, q, X);
}

