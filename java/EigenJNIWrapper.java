package dr.evomodel.substmodel.eigen;

/*
 * EigenJNIWrapper
 *
 * @author Xiang Ji
 *
 */
public class EigenJNIWrapper {
    public EigenJNIWrapper () {
    }

    static {
        System.load("/usr/local/lib/libeigen-jni.jnilib");
    }

    public native int createInstance(int matrixCount,
                                     int stateCount);

    public native int setMatrix(int matrix,
                                int[] indices,
                                double[] values,
                                int nonZeroCount,
                                int stateCount);

    public native String getVersion();

    public native int getEigenDecomposition(int matrix,
                                            double[] eigenValues,
                                            double[] eigenVectors,
                                            double[] inverseEigenVectors);
}
