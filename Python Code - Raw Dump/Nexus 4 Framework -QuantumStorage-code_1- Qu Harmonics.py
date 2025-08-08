using System;

namespace HarmonicCompression
{
    class Program
    {
        const double HarmonicConstant = 0.35;
        const double Tolerance = 0.01;
        const int MaxIterations = 100;

        static void Main(string[] args)
        {
            // Example input: a large string
            string inputText = string.Join(" ", new string[] { "Quantum Storage Test" }, 50);

            Console.WriteLine("Original Text Length: " + inputText.Length);
            Console.WriteLine("Original Text: " + inputText);

            // Step 1: Base-16 Encoding
            string base16Data = TextToBase16(inputText);
            Console.WriteLine("Base-16 Encoded Data Length: " + base16Data.Length);

            // Step 2: Convert Base-16 Data to Binary
            string binaryData = Base16ToBinary(base16Data);
            int originalBinaryLength = binaryData.Length;
            Console.WriteLine("Binary Data Length: " + originalBinaryLength);

            // Step 3: Prepare Data Matrix
            int matrixSize = (int)Math.Ceiling(Math.Sqrt(originalBinaryLength));
            double[,] dataMatrix = CreateBinaryMatrix(binaryData, matrixSize);

            // Step 4: Harmonic Compression
            double[,] compressedData = HarmonizeData(dataMatrix, HarmonicConstant, true);
            int compressedSize = CountNonZeroElements(compressedData);
            Console.WriteLine("Compressed Data Non-Zero Elements: " + compressedSize);

            // Step 5: Harmonic Expansion
            double[,] expandedData = HarmonizeData(compressedData, HarmonicConstant, false);
            string flattenedBinary = FlattenMatrix(expandedData, originalBinaryLength);

            // Step 6: Convert Binary Back to Text
            string recoveredText = BinaryToText(flattenedBinary);
            Console.WriteLine("Recovered Text: " + recoveredText);

            // Compression Ratio
            double compressionRatio = (double)compressedSize / originalBinaryLength;
            Console.WriteLine("Compression Ratio (Compressed/Original): " + compressionRatio);
        }

        static string TextToBase16(string text)
        {
            char[] chars = text.ToCharArray();
            string result = string.Empty;
            foreach (char c in chars)
                result += Convert.ToString(c, 16).PadLeft(2, '0');
            return result;
        }

        static string Base16ToBinary(string base16Data)
        {
            string binary = string.Empty;
            foreach (char c in base16Data)
                binary += Convert.ToString(Convert.ToInt32(c.ToString(), 16), 2).PadLeft(4, '0');
            return binary;
        }

        static double[,] CreateBinaryMatrix(string binaryData, int matrixSize)
        {
            double[,] matrix = new double[matrixSize, matrixSize];
            for (int i = 0; i < binaryData.Length; i++)
            {
                int row = i / matrixSize;
                int col = i % matrixSize;
                matrix[row, col] = binaryData[i] == '1' ? 1 : 0;
            }
            return matrix;
        }

        static double[,] HarmonizeData(double[,] data, double harmonicConstant, bool compress)
        {
            double[,] result = (double[,])data.Clone();
            int rows = result.GetLength(0);
            int cols = result.GetLength(1);
            double gain = compress ? 1.0 : -1.0;

            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                double delta = CalculateMatrixMean(result) - harmonicConstant;
                if (Math.Abs(delta) < Tolerance) break;

                double adjustment = delta * gain;

                for (int row = 0; row < rows; row++)
                    for (int col = 0; col < cols; col++)
                        result[row, col] = Math.Clamp(result[row, col] - adjustment, 0.0, 1.0);
            }

            return result;
        }

        static double CalculateMatrixMean(double[,] matrix)
        {
            double sum = 0.0;
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            for (int row = 0; row < rows; row++)
                for (int col = 0; col < cols; col++)
                    sum += matrix[row, col];
            return sum / (rows * cols);
        }

        static int CountNonZeroElements(double[,] matrix)
        {
            int count = 0;
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            for (int row = 0; row < rows; row++)
                for (int col = 0; col < cols; col++)
                    if (matrix[row, col] != 0)
                        count++;
            return count;
        }

        static string FlattenMatrix(double[,] matrix, int length)
        {
            string result = string.Empty;
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            for (int i = 0; i < length; i++)
            {
                int row = i / cols;
                int col = i % cols;
                result += (int)Math.Round(matrix[row, col]);
            }
            return result;
        }

        static string BinaryToText(string binary)
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
