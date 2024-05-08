//
// Created by Xiang Ji on 5/7/24.
//
#include "dr_evomodel_substmodel_eigen_EigenJNIWrapper.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

#include <iostream>

//#define EIGEN_DEBUG_FLOW
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::EigenSolver<MatrixXd> EigenSolver;

class EigenImpl {
private:
    MatrixXd* matrices;
    int kState;
    EigenSolver* eigenSolvers;

    enum EigenReturnCode {
        EIGEN_SUCCESS = 0,
    };

public:
    int createInstance(int matrixCount,
                       int stateCount) {
        matrices = new MatrixXd [matrixCount];
        eigenSolvers = new EigenSolver [matrixCount];
        kState = stateCount;
        for (int i = 0; i < matrixCount; i++) {
            MatrixXd  matrix(stateCount, stateCount);
            matrix.setZero();
            matrices[i] = matrix;
        }
        return EIGEN_SUCCESS;
    }

    int setMatrix(int matrix,
                  int* indices,
                  double* values,
                  int nonZeroCount) {
#ifdef EIGEN_DEBUG_FLOW
        std::cerr<< "Entering setMatrix" << std::endl;
#endif
        matrices[matrix].setZero();

        for (int i = 0; i < nonZeroCount; i++) {
            matrices[matrix](indices[2 * i], indices[2 * i + 1]) = values[i];
        }
#ifdef EIGEN_DEBUG_FLOW
        std::cerr<< "Matrix:" << std::endl << matrices[matrix] << std::endl <<std::endl;
#endif
        return EIGEN_SUCCESS;
    }

    int getEigenDecomposition(int matrix,
                              double* eigenValues,
                              double* eigenVectors,
                              double* inverseEigenVectors) {
#ifdef EIGEN_DEBUG_FLOW
        std::cerr<< "Entering getEigenDecomposition:" << std::endl;
#endif

        EigenSolver es(matrices[matrix]);
        eigenSolvers[matrix] = es;

        MatrixXd V = es.pseudoEigenvectors();
        MatrixXd inverseV = es.pseudoEigenvectors().inverse();
#ifdef EIGEN_DEBUG_FLOW
        std::cerr<< "The eigenvalues of Matrix are:" << std::endl << es.pseudoEigenvalueMatrix() << std::endl;
        std::cerr<< "The matrix of eigenvectors, V, is:" << std::endl << es.pseudoEigenvectors() << std::endl << std::endl;
        std::cerr<< "The matrix of inverse eigenvectors, invV, is:" << std::endl << inverseV << std::endl << std::endl;
#endif
        bool complexPair = false;
        for (int i = 0; i < kState; i++) {
            eigenValues[i] = es.pseudoEigenvalueMatrix().col(i)[i];
            if (complexPair) {
                eigenValues[kState + i] = es.pseudoEigenvalueMatrix().col(i - 1)[i];
                complexPair = false;
            } else {
                double nextValue = i < kState - 1 ? es.pseudoEigenvalueMatrix().col(i + 1)[i] : 0;
                if (nextValue != 0) {
                    complexPair = true;
                    eigenValues[kState + i] = nextValue;
                } else {
                    eigenValues[kState + i] = 0;
                }
            }
            for (int j = 0; j < kState; j++) {
                eigenVectors[kState * i + j] = V.col(j)[i];
                inverseEigenVectors[kState * i + j] = inverseV.col(j)[i];
            }
#ifdef EIGEN_DEBUG_FLOW
            std::cerr<< "The " << i <<"th eigenvalue :" << std::endl << eigenValues[i] << "," << eigenValues[i + kState] << std::endl;
#endif

        }

#ifdef EIGEN_DEBUG_FLOW
        std::cerr<< "The eigenvalues of Matrix are:" << std::endl << es.pseudoEigenvalueMatrix() << std::endl;
        std::cerr<< "The matrix of eigenvectors, V, is:" << std::endl << es.pseudoEigenvectors() << std::endl << std::endl;
        std::cerr<< "The matrix of inverse eigenvectors, invV, is:" << std::endl << inverseV << std::endl << std::endl;
#endif

        return EIGEN_SUCCESS;
    }

    static EigenImpl* factory(int matrixCount,
                                int stateCount,
                                int* errorCode) {
        EigenImpl* impl = new EigenImpl();
        try {
            *errorCode = impl->createInstance(matrixCount, stateCount);
            if (*errorCode == EIGEN_SUCCESS) {
                return impl;
            }
            delete impl;
            return NULL;
        }
        catch(...) {
            delete impl;
            throw;
        }

        delete impl;
        return NULL;
    };


    ~EigenImpl() {
    };

};

typedef std::shared_ptr<EigenImpl> InstancePtr;
std::vector<InstancePtr> instances;

/*
 * Class:     dr_evomodel_substmodel_eigen_EigenJNIWrapper
 * Method:    createInstance
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_eigen_EigenJNIWrapper_createInstance
        (JNIEnv * env, jobject obj, jint matrixCount, jint stateCount) {
    jint errorCode;
    instances.emplace_back(
        EigenImpl::factory(matrixCount, stateCount, &errorCode)
    );
    return errorCode;
}

/*
 * Class:     dr_evomodel_substmodel_eigen_EigenJNIWrapper
 * Method:    setMatrix
 * Signature: (I[I[DI)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_eigen_EigenJNIWrapper_setMatrix
        (JNIEnv * env, jobject obj, jint matrix, jintArray inIndices, jdoubleArray inValues, jint nonZeroCount) {
    jint *indices = env -> GetIntArrayElements(inIndices, NULL);
    jdouble *values = env ->GetDoubleArrayElements(inValues, NULL);

    jint errorCode = (jint) instances[0] -> setMatrix(matrix, indices, values, nonZeroCount);

    env -> ReleaseDoubleArrayElements(inValues, values, JNI_ABORT);
    env ->ReleaseIntArrayElements(inIndices, indices, JNI_ABORT);
    return errorCode;
}

/*
 * Class:     dr_evomodel_substmodel_eigen_EigenJNIWrapper
 * Method:    getVersion
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_dr_evomodel_substmodel_eigen_EigenJNIWrapper_getVersion
        (JNIEnv *env, jobject obj) {
    return env->NewStringUTF("0.0.1");
}

/*
 * Class:     dr_evomodel_substmodel_eigen_EigenJNIWrapper
 * Method:    getEigenDecomposition
 * Signature: (I[D[D[D)I
 */
JNIEXPORT jint JNICALL Java_dr_evomodel_substmodel_eigen_EigenJNIWrapper_getEigenDecomposition
        (JNIEnv * env, jobject obj, jint matrix, jdoubleArray outEigenValues, jdoubleArray outEigenVectors, jdoubleArray outInverseEigenVectors) {
    jdouble * eigenValues = env -> GetDoubleArrayElements(outEigenValues, NULL);
    jdouble * eigenVectors = env -> GetDoubleArrayElements(outEigenVectors, NULL);
    jdouble * inverseEigenVectors = env -> GetDoubleArrayElements(outInverseEigenVectors, NULL);

    jint errorCode = (jint) instances[0] ->getEigenDecomposition(matrix, eigenValues, eigenVectors, inverseEigenVectors);

    env->ReleaseDoubleArrayElements(outEigenValues, eigenValues, 0);
    env->ReleaseDoubleArrayElements(outEigenVectors, eigenVectors, 0);
    env->ReleaseDoubleArrayElements(outInverseEigenVectors, inverseEigenVectors, 0);

    return errorCode;
}



