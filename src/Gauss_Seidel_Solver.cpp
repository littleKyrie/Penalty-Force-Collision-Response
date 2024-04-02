//
// Created by Kyrie Zhang on 2023/10/23.
//

/*
#include <Gauss_Seidel_Solver.h>

template <typename dataType>
void gaussSeidel(Eigen::Matrix<dataType, Eigen::Dynamic, 1>& result,
                 const Eigen::SparseMatrix<dataType>& A,
                 const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b,
                 bool isConcurrent,
                 dataType relaxation,
                 size_t maxSteps,
                 dataType tol) {
    // the matrix A is s.p.d, so the solver will converge
    // the matrix A is not singular, so the system can be solved

    // preparation
    size_t size = b.size();
    result.setZero();
    result.resize(size);

    // initial guess & temple result
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x0 = result;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x = x0;

    if (isConcurrent) {
        // for parallel reduce
        struct SumEle {
            dataType sumEle = 0;
            Eigen::SparseMatrix<dataType>& rowVec;
            Eigen::Matrix<dataType, Eigen::Dynamic, 1>& solution;

            SumEle(const SumEle& x, tbb::split)
                    : sumEle(x.sumEle), rowVec(x.rowVec), solution(x.solution) {}

            void operator()(const tbb::blocked_range<size_t>& r) {
                dataType sumLocalEle = sumEle;
                for (size_t j = r.begin(); j < r.end(); j++) {
                    sumLocalEle += rowVec.coeff(1, j) * solution[j];
                }
                sumEle = sumLocalEle;
            }

            void join(const SumEle& y) {sumEle += y.sumEle;}

            SumEle(const Eigen::SparseMatrix<dataType>& sparseVec,
                   const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& solutionVec)
                    : rowVec(sparseVec), sumEle(0), solution(solutionVec) {}
        };

        // iteration
        size_t iterator = 0;
        while (iterator < maxSteps) {
            for (size_t i = 0; i < size; i++) {
                SumEle newSum(A.row(i), x), oldSum(A.row(i), x0);

                // new sum
                tbb::parallel_reduce(tbb::blocked_range<size_t>(0, i), newSum);

                // old sum
                tbb::parallel_reduce(tbb::blocked_range<size_t>(i+1, size), oldSum);

                // get new xi
                x[i] = x0[i] + relaxation * (b[i] - newSum.sumEle - oldSum.sumEle - A.coeff(i, i) * x0[i]) / A.coeff(i, i);
            }

            // sufficiently converge
            Eigen::Matrix<dataType, Eigen::Dynamic, 1> res = (x - x0);
            if (res.norm() < tol)
                break;

            x0 = x;

            iterator++;
        }
    } else {
        // iteration
        size_t iteration = 0;
        while (iteration < maxSteps) {

            for (size_t i = 0; i < size; i++) {
                dataType newSum = 0;
                dataType oldSum = 0;

                for (size_t j = 0; j < size; j++) {
                    if (j < i) {
                        newSum += A.coeff(i, j) * x[j];
                    } else if (j > i) {
                        oldSum += A.coeff(i, j) * x0[j];
                    } else {
                        continue;
                    }
                }

                x[i] = x0[i] + relaxation * (b[i] - newSum - oldSum - A.coeff(i, i) * x0[i]) / A.coeff(i, i);
            }

            // sufficiently converge
            Eigen::Matrix<dataType, Eigen::Dynamic, 1> res = (x - x0);
            if (res.norm() < tol)
                break;

            x0 = x;

            iteration++;
        }
    }

    // write back
    result = x;
}
*/