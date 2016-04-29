using System;
using FiniteElementMethodPE.Helpers;

namespace FiniteElementMethodPE.FiniteElements
{
    internal class FiniteElementPuassonESolver
    {
        private const int AxisPoints = 3;
        private const int AllPointsNumber = AxisPoints*AxisPoints;
        // межі області
        private readonly double a1, b1, a2, b2;
        private double[,] A; // загальна матриця
        private double[] B; // права частина
        private readonly Func<double, double, double> f;
        private double h;
        private double[,] Ke = new double[AllPointsNumber, AllPointsNumber];
        private double[,] Me = new double[AllPointsNumber, AllPointsNumber];
        private int N;
        private double[] Qe = new double[AllPointsNumber];

        private readonly double[] tempVector = new double[AllPointsNumber];

        //область має бути квадратна
        public FiniteElementPuassonESolver(Func<double, double, double> f, double a1, double b1, double a2, double b2)
        {
            this.a1 = a1;
            this.b1 = b1;
            this.a2 = a2;
            this.b2 = b2;
            this.f = f;
        }

        public double[] Solve(int n)
        {
            N = n;
            h = (b1 - a1)/N;
            if (Math.Abs(h - (b2 - a2)/N) > 0.001) throw new Exception("Must be square");
            int sizeofsystem = (2*n + 1)*(2*n + 1);
            A = new double[sizeofsystem, sizeofsystem];
            B = new double[sizeofsystem];
            for (var i = 0; i < sizeofsystem; i++)
            {
                for (var j = 0; j < sizeofsystem; j++)
                {
                    A[i, j] = 0;
                }
                B[i] = 0;
            }

            for (var p = 0; p < n*n; p++)
            {
                int sigma0 = p/n;
                int sigma1 = p%n;
                double x1a = a1 + sigma1*h;
                double x2a = a2 + sigma0*h;
                int l = sigma0*2*(2*n + 1) + 2*sigma1;
                FillKe(ref Ke, x1a, x2a);
                FillQe(ref Qe, x1a, x2a);

                for (var i = 0; i < AllPointsNumber; i++)
                {
                    int sigma2 = i/3;
                    int sigma3 = i%3;
                    for (var j = 0; j < AllPointsNumber; j++)
                    {
                        int sigma4 = j/3;
                        int sigma5 = j%3;
                        A[l + sigma2*(2*n + 1) + sigma3, l + sigma4*(2*n + 1) + sigma5] += Ke[i, j];
                    }
                    B[l + sigma2*(2*n + 1) + sigma3] += Qe[i];
                }
            }

            // Write marix and right part to file
            Printer.WriteLine("Results.txt", "Matrix : ", false);
            Printer.Write("Results.txt", A, sizeofsystem, sizeofsystem, true);
            Printer.WriteLine("Results.txt", "Right vector : ", true);
            Printer.Write("Results.txt", B, true);

            // Delete all bourder points from matrix A and vector B
            int oldRowSize = 2*n + 1;
            int newRowSize = 2*n - 1;
            int newSize = newRowSize * newRowSize;
            double[,] innerA = new double[newSize, newSize];
            double[] innerB = new double[newSize];

            int indexI = 0, indexJ = 0;
            for (int i = 0; i < sizeofsystem; i++)
            {
                if( i < oldRowSize) continue;
                if( i >= sizeofsystem - oldRowSize) continue;
                if((i % oldRowSize == 0) || (i% oldRowSize == oldRowSize-1)) continue;
                indexJ = 0;
                for (int j = 0; j < sizeofsystem; j++)
                {
                    if (j < oldRowSize) continue;
                    if (j >= sizeofsystem - oldRowSize) continue;
                    if ((j % oldRowSize == 0) || (j % oldRowSize == oldRowSize - 1)) continue;
                    innerA[indexI, indexJ] = A[i, j];
                    indexJ++;
                }
                innerB[indexI] = B[i];
                indexI++;
            }

            //result = SystemOfLinearEquations.SolveWithQRmethod(A, B, B.Length);
            double[] result = SystemOfLinearEquations.SolveWithQRmethod(innerA, innerB, innerB.Length);

            return result;
        }

        // нижче використовуйте клас Integral з Helpers щоб обчислити інтеграли
        // розмірність матриці Ke 9x9 
        // (x1a,x2a) - лівий нижній кут
        private void FillKe(ref double[,] Ke, double x1a, double x2a) // розмірність матриці Ke 9x9 
        {
            // заповнюємо матрицю Ke для кокретної області [x1a,x1b]x[x2a,x2b]
            var E = 0.0005;
            for (var i = 0; i < AllPointsNumber; i++)
            {
                for (var j = 0; j < AllPointsNumber; j++)
                {
                    int sigma0 = i/3;
                    int sigma1 = j/3;
                    int sigma2 = i%3;
                    int sigma3 = j%3;
                    var core = new Core(LagrangeFunc, LagrangeFunc);
                    var core1 = new Core(DLanrangeFunc, DLanrangeFunc);
                    double temp = 0;
                    core1.SetParams(sigma0, x1a, sigma1, x1a);
                    temp = Integral.CalculateWithHauseMethod(core1, x1a, x1a + h, E, 5);
                    core.SetParams(sigma2, x2a, sigma3, x2a);
                    temp *= Integral.CalculateWithHauseMethod(core, x2a, x2a + h, E, 5);
                    Ke[i, j] = temp;
                    core.SetParams(sigma0, x1a, sigma1, x1a);
                    temp = Integral.CalculateWithHauseMethod(core, x1a, x1a + h, E, 5);
                    core1.SetParams(sigma2, x2a, sigma3, x2a);
                    temp *= Integral.CalculateWithHauseMethod(core1, x2a, x2a + h, E, 5);
                    Ke[i, j] += temp;
                }
            }
        }

        private void FillQe(ref double[] Qe, double x1a, double x2a)
        {
            var index = 0;
            for (var i = 0; i < AxisPoints; i++)
            {
                for (var j = 0; j < AxisPoints; j++)
                {
                    tempVector[index] = f(x1a + j*h/2.0, x2a + i*h/2.0);
                    index ++;
                }
            }
            // множення матриці на вектор
            FillMe(ref Me, x1a, x2a);
            double sum = 0;
            for (var i = 0; i < AllPointsNumber; i++)
            {
                sum = 0;
                for (var j = 0; j < AllPointsNumber; j++)
                {
                    sum += Me[i, j]*tempVector[j];
                }
                Qe[i] = sum;
            }
        }

        private void FillMe(ref double[,] Me, double x1a, double x2a)
        {
            // заповнюємо матрицю Me для кокретної області [x1a,x1b]x[x2a,x2b]
            var E = 0.0005;
            for (var i = 0; i < AllPointsNumber; i++)
            {
                for (var j = 0; j < AllPointsNumber; j++)
                {
                    int sigma0 = i/3;
                    int sigma1 = j/3;
                    int sigma2 = i%3;
                    int sigma3 = j%3;
                    var core = new Core(LagrangeFunc, LagrangeFunc);
                    core.SetParams(sigma0, x1a, sigma1, x1a);
                    Me[i, j] = Integral.CalculateWithHauseMethod(core, x1a, x1a + h, E, 5);
                    core.SetParams(sigma2, x2a, sigma3, x2a);
                    Me[i, j] *= Integral.CalculateWithHauseMethod(core, x2a, x2a + h, E, 5);
                }
            }
        }

        // Прописати всі потрібні базисні фукнції l і їх похідні
        private double LagrangeFunc(double x, int i, double xa)
        {
            var chus = 1.0;
            var dob = 1.0;
            if (x < xa || x > xa + h) return 0;
            for (var j = 0; j < 3; j++)
            {
                if (i == j) continue;
                chus *= x - (xa + h*j/2.0);
                dob *= (i - j)*h/2.0;
            }
            return chus/dob;
        }

        private double DLanrangeFunc(double x, int i, double xa) // похідна функції Лагранжа
        {
            double sum = 0;
            var dob = 1.0;
            if (x < xa || x > xa + h) return 0;
            for (var j = 0; j < 3; j++)
            {
                if (i == j) continue;
                sum += x - (xa + h*j/2.0);
                dob *= (i - j)*h/2.0;
            }
            return sum/dob;
        }
    }
}