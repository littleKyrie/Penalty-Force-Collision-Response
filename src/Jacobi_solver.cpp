//
// Created by Kyrie Zhang on 2023/10/21.
//

/*
#include <Jacobi_solver.h>

template <typename dataType>
bool jacobiSolver (
        Eigen::Matrix<dataType, Eigen::Dynamic, 1>& result,
        const Eigen::SparseMatrix<dataType>& A,
        const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b,
        bool isConcurrent,
        dataType relaxation,
        size_t maxSteps,
        dataType tol) {
    size_t size = b.size();


    // can be solved ?
    tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            if (A.coeff(i, i) == 0) {
                std::cout << "the system matrix is a singular one, can not be solved " << std::endl;
                return;
            }
        }
    });

    // can converge ?
    // s.t.d
    bool diagonallyDominant = true;
    // bool symmetric;
    // bool positiveDefine;

    std::mutex mtx;

    struct SumEle {
        dataType sum = 0;
        const Eigen::SparseMatrix<dataType>& rowVec;
        const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& solution;
        size_t i;

        SumEle(SumEle& x, tbb::split)
                : sum(x.sum), rowVec(x.rowVec), solution(x.solution), i(x.i) {}

        void operator()(const tbb::blocked_range<size_t>& l) {
            dataType localSum = sum;
            for (size_t j = l.begin(); j < l.end(); j++) {
                if (j != i) {
                    localSum += rowVec.coeff(1, j) * solution[j];
                }
            }
            sum = localSum;
        }

        void join(const SumEle& y) {sum += y.sum;}

        SumEle(const Eigen::SparseMatrix<dataType>& sparseVec,
               const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& solutionVec,
               size_t index)
                : rowVec(sparseVec), sum(0), solution(solutionVec), i(index) {}
    };


    tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            tbb::this_task_arena::isolate([&]{
                SumEle s(A.row(i), i);
                tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size), s);
                if (abs(A.coeff(i, i)) <= s.sum) {
                    std::unique_lock<std::mutex> lck(mtx);
                    diagonallyDominant = false;
                }
            });
        }
    });


    // initial guess
    result.setZero();
    result.resize(size);
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x0 = result;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x = x0;

    size_t iterator = 0;
    while(iterator < maxSteps) {

        if (isConcurrent) {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t>& r){
                for (size_t j = r.begin(); j < r.end(); j++) {
                    SumEle s(A.row(j), x0, j);
                    tbb::this_task_arena::isolate([&]{
                        tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size), s);
                    });
                    x[j] = x0[j] + relaxation / A.coeff(j, j) * (b[j] - s.sum - A.coeff(j, j) * x0[j]);

                    if (std::isnan(x[j])) {
                        std::cout << "numerical issue " << std::endl;
                        return false;
                    }
                }
            });
        } else {
            for (size_t i = 0; i < size; i++) {
                dataType sum = 0;
                for (size_t j = 0; j < size; j++) {

                    if (i != j) {
                        sum += A.coeff(i, j) * x0[j];
                    } else {
                        continue;
                    }
                }

                x[i] = x0[i] + relaxation / A.coeff(i, i) * (b[i] - sum - A.coeff(i, i) * x0[i]);
                if (std::isnan(x[i])) {
                    std::cout << "numerical issue " << std::endl;
                    return false;
                }
            }
        }

        Eigen::Matrix<dataType, Eigen::Dynamic, 1> res = x - x0;
        x0 = x;
        if (res.norm() < tol)
            break;

        iterator++;

        if (iterator == maxSteps) {
            std::cout << "the last residual of iteration is " << res.norm() << std::endl;
            return false;
        }
    }

    result = x;

    return true;

}
*/