#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;


int main() {
    MatrixXd A(4,4);
    A << -0.7265668659069038, 0.22270454583996452, 0.2877257277577091, 0.21613659230923019,
            0.2734331340930963, -0.7772954541600356, 0.2877257277577091, 0.21613659230923019,
            0.2734331340930963, 0.22270454583996452, -0.712274272242291, 0.21613659230923019,
            0.2734331340930963, 0.22270454583996452, 0.2877257277577091, -0.7838634076907699;

    MatrixXd B(3721, 3721);
    B.setZero();
    int index0 = 0, index1 = 100, index2 = 200, index3 = 300;
    B.col(index0)[index0] = A.col(0)[0];
    B.col(index0)[index1] = A.col(0)[1];
    B.col(index0)[index2] = A.col(0)[2];
    B.col(index0)[index3] = A.col(0)[3];

    B.col(index1)[index0] = A.col(1)[0];
    B.col(index1)[index1] = A.col(1)[1];
    B.col(index1)[index2] = A.col(1)[2];
    B.col(index1)[index3] = A.col(1)[3];

    B.col(index2)[index0] = A.col(2)[0];
    B.col(index2)[index1] = A.col(2)[1];
    B.col(index2)[index2] = A.col(2)[2];
    B.col(index2)[index3] = A.col(2)[3];

    B.col(index3)[index0] = A.col(3)[0];
    B.col(index3)[index1] = A.col(3)[1];
    B.col(index3)[index2] = A.col(3)[2];
    B.col(index3)[index3] = A.col(3)[3];
//    std::cout << B << std::endl;
    std::cout << "Here is a 4x4 matrix, A:" << std::endl << A << std::endl << std::endl;

    Eigen::EigenSolver<MatrixXd> es(A);

    std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
    std::cout<< "Test: " << es.eigenvalues()[0] <<  std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
    std::cout << "Test vector 0 :" << std::endl << es.eigenvectors().col(0) << std::endl;
    std::cout << "Test vector 1 :" << std::endl << es.eigenvectors().col(1) << std::endl;

    std::cout << "Pseudo eigenvalues: " << std::endl << es.pseudoEigenvalueMatrix() << std::endl;
    std::cout << "Pseudo eigenvectors: " << std::endl << es.pseudoEigenvectors() << std::endl;
    std::cout << "test (0, 1) of pseudo eigenvectors: " << es.pseudoEigenvectors().col(0)[1] << std::endl;
    Eigen::MatrixXcd sD = es.pseudoEigenvalueMatrix();
    Eigen::MatrixXcd sV = es.pseudoEigenvectors();
    std::cout << "Pseudo V * D * V^(-1) = " << std::endl << sV * sD * sV.inverse() << std::endl;

    Eigen::EigenSolver<MatrixXd> esB(B);
    std::cout << "Pseudo eigenvalues (diagonal): " << std::endl << esB.pseudoEigenvalueMatrix().diagonal() << std::endl;

    std::complex<double> lambda = es.eigenvalues()[0];
    std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
    Eigen::VectorXcd v = es.eigenvectors().col(0);
    std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
    std::cout << "... and A * v = " << std::endl << A.cast<std::complex<double> >() * v << std::endl << std::endl;

    Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
    Eigen::MatrixXcd V = es.eigenvectors();
    std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;
}
