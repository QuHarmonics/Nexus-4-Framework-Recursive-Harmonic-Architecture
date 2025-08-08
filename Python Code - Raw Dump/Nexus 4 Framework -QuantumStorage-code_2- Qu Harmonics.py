using System;
using MathNet.Numerics.LinearAlgebra;

namespace HarmonicCompressionApp
{
    public class HarmonicCompression
    {
        // Constants for harmonic compression
        private const double HarmonicConstant = 0.35;
        private const double Tolerance = 0.01;
        private const int MaxIterations = 100;

        /// <summary>
        /// Compresses input text using harmonic compression and returns the compressed representation.
        /// </summary>
        public double[,] Compress(string inputText, out string base16Data, out string binaryData, out int originalLength)
        {
            // Step 1: Base-16 encoding
            base16Data = TextToBase16(inputText);

            // Step 2: Convert Base-16 to binary
            binaryData = Base16ToBinary(base16Data);
            originalLength = binaryData.Length;

            // Step 3: Prepare binary matrix
            int matrixSize = (int)Math.Ceiling(Math.Sqrt(originalLength));
            var dataMatrix = CreateBinaryMatrix(binaryData, matrixSize);

            // Step 4: Perform harmonic compression
            return HarmonizeData(dataMatrix, HarmonicConstant, true);
        }

        /// <summary>
        /// Expands a compressed matrix back to the original text.
        /// </summary>
        public string Expand(double[,] compressedMatrix, int originalLength)
        {
            // Step 5: Perform harmonic expansion
            var expandedMatrix = HarmonizeData(compressedMatrix, HarmonicConstant, false);

            // Step 6: Flatten and convert to binary string
            string flattenedBinary = FlattenMatrix(expandedMatrix, originalLength);

            // Step 7: Convert binary back to text
            return BinaryToText(flattenedBinary);
        }

        /// <summary>
        /// Converts input text to its Base-16 representation.
        /// </summary>
        private string TextToBase16(string text)
        {
            char[] chars = text.ToCharArray();
            string result = string.Empty;
            foreach (char c in chars)
                result += Convert.ToString(c, 16).PadLeft(2, '0');
            return result;
        }

        /// <summary>
        /// Converts Base-16 string to binary representation.
        /// </summary>
        private string Base16ToBinary(string base16Data)
        {
            string binary = string.Empty;
            foreach (char c in base16Data)
                binary += Convert.ToString(Convert.ToInt32(c.ToString(), 16), 2).PadLeft(4, '0');
            return binary;
        }

        /// <summary>
        /// Creates a binary matrix from a binary string.
        /// </summary>
        private Matrix<double> CreateBinaryMatrix(string binaryData, int matrixSize)
        {
            var matrix = Matrix<double>.Build.Dense(matrixSize, matrixSize, 0);
            for (int i = 0; i < binaryData.Length; i++)
            {
                int row = i / matrixSize;
                int col = i % matrixSize;
                matrix[row, col] = binaryData[i] == '1' ? 1 : 0;
            }
            return matrix;
        }

        /// <summary>
        /// Harmonizes the data matrix for compression or expansion.
        /// </summary>
        private Matrix<double> HarmonizeData(Matrix<double> data, double harmonicConstant, bool compress)
        {
            var result = data.Clone();
            double gain = compress ? 1.0 : -1.0;

            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                double delta = result.RowSums().Average() - harmonicConstant;
                if (Math.Abs(delta) < Tolerance) break;

                result = result.Map(x => Math.Clamp(x - delta * gain, 0.0, 1.0));
            }

            return result;
        }

        /// <summary>
        /// Flattens a matrix back into a binary string.
        /// </summary>
        private string FlattenMatrix(Matrix<double> matrix, int length)
        {
            string result = string.Empty;
            int count = 0;
            foreach (var value in matrix.Enumerate())
            {
                if (count >= length) break;
                result += (int)Math.Round(value);
                count++;
            }
            return result;
        }

        /// <summary>
        /// Converts binary string back to text.
        /// </summary>
        private string BinaryToText(string binary)
        {
            string result = string.Empty;
            for (int i = 0; i < binary.Length; i += 8)
            {
                if (i + 8 > binary.Length) break;
                string byteStr = binary.Substring(i, 8);
                result += (char)Convert.ToInt32(byteStr, 2);
            }
            return result;
        }
    }
}
