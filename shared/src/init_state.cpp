#include <init_state.h>

//void init_state(Eigen::Ref<Eigen::VectorXd> q, Eigen::Ref<Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> V)
void init_state(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Ref<const Eigen::MatrixXd> V) {
    q.resize(V.rows()*V.cols());
    qdot.resize(V.rows()*V.cols());

    Eigen::MatrixXd Vt = V.transpose();
    q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows()*Vt.cols());
    qdot.setZero();

}

void initMultiState(Eigen::VectorXd &q, Eigen::VectorXd &qdot, const std::vector<Eigen::VectorXd>& localQ, const std::vector<Eigen::VectorXd>& localQdot) {
    int initialIndex = 0;
    size_t numOfObjects = localQ.size();
    int size = 0;
    for (size_t i = 0; i < numOfObjects; i++) {
        size += localQ[i].size();
    }
    q.resize(size);
    qdot.resize(size);
    for (size_t i = 0 ; i < numOfObjects; i++) {
        q.segment(initialIndex, localQ[i].size()) = localQ[i];
        qdot.segment(initialIndex, localQdot[i].size()) = localQdot[i];
        initialIndex += localQ[i].size();
    }
}