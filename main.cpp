#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;


int main() {
    MatrixXd A = MatrixXd::Random(6,6);
    std::cout << "Here is a random 6x6 matrix, A:" << std::endl << A << std::endl << std::endl;

    Eigen::EigenSolver<MatrixXd> es(A);
    std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
    std::cout<< "Test: " << es.eigenvalues()[0] <<  std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
    std::cout << "Test vector 0 :" << std::endl << es.eigenvectors().col(0) << std::endl;
    std::cout << "Test vector 1 :" << std::endl << es.eigenvectors().col(1) << std::endl;

    std::complex<double> lambda = es.eigenvalues()[0];
    std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
    Eigen::VectorXcd v = es.eigenvectors().col(0);
    std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
    std::cout << "... and A * v = " << std::endl << A.cast<std::complex<double> >() * v << std::endl << std::endl;

    Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
    Eigen::MatrixXcd V = es.eigenvectors();
    std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;
}
