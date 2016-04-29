using System.IO;

namespace FiniteElementMethodPE.Helpers
{
    internal class Printer
    {
        public static void WriteLine(string fileName, string text, bool append)
        {
            using (var stream = new StreamWriter(fileName, append))
            {
                stream.WriteLine(text);
            }
        }

        public static void Write(string fileName, double[,] A, int rows, int columns, bool append)
        {
            using (var stream = new StreamWriter(fileName, append))
            {
                for (var i = 0; i < rows; i++)
                {
                    for (var j = 0; j < columns; j++)
                    {
                        stream.Write("{0:E2}  ", A[i, j]);
                        if((j + 1) % 3 == 0) stream.Write("\t");
                    }
                    if ((i + 1) % 3 == 0) stream.WriteLine();
                    stream.WriteLine();
                }
            }
        }

        public static void Write(string fileName, double[] vector, bool append)
        {
            using (var stream = new StreamWriter(fileName, append))
            {
                for (var i = 0; i < vector.Length; i++)
                {
                    stream.Write("{0:E2}  ", vector[i]);
                }
            }
        }
    }
}