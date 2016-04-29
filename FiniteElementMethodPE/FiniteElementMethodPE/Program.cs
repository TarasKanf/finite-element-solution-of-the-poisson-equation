using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks; 
using FiniteElementMethodPE.FiniteElements;

namespace FiniteElementMethodPE
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Bottom left corner");
            Console.Write(" x: 0 \n");
            double x = 0;//int.Parse(Console.ReadLine());
            Console.Write(" y: 0 \n");
            double y = 0;//int.Parse(Console.ReadLine());
            Console.Write(" H: 1 \n");
            double h = 1;//int.Parse(Console.ReadLine());
            Console.Write(" N : ");
            int n = int.Parse(Console.ReadLine());
            FiniteElementPuassonESolver solver = new FiniteElementPuassonESolver(F, x, x+h, y, y+ h);
            double[] result = solver.Solve(n);

            double maxDeviation = 0;
            double step = h / (2 * n);
            double a = x + step;
            double b = y + step;
            Console.WriteLine("Result without border:");
            for (int i = 0, index = 0; i < 2 * n - 1; i++)
            {
                for (int j = 0; j < 2 * n - 1; j++)
                {
                    Console.Write("{0:E}  ", result[index]);
                    maxDeviation = Math.Max(maxDeviation,
                        Math.Abs((result[index] - AccurateSolution(a + i*step, b + j*step))/result[index]));
                    index++;
                }
                Console.WriteLine();
            }

            Console.WriteLine("Accurate solution without border: ");
            for (int i = 0; i < 2 * n - 1; i++)
            {
                for (int j = 0; j < 2 * n - 1; j++)
                {
                    Console.Write("{0:E} ({1})({2})  ", AccurateSolution(a + i * step, b + j * step), a + i * step, b + j * step);
                }
                Console.WriteLine();
            }

            Console.WriteLine("Max deviation = {0}",maxDeviation);
            Console.ReadKey();

        }

        // f при u* = xy^2(1-x)^2(1-y)
        public static double F(double x, double y)
        {
            return -2 * Math.Pow(1 - x,2) * x * (1 - y) 
                +4 * Math.Pow(1 - x,2) * x * y 
                +4 * (1 - x) * (1 - y) * y * y 
                - 2 * x * (1 - y) * y *y;
        }

        public static double AccurateSolution(double y, double x)
        {
            return x*y*y*Math.Pow(1 - x,2)*(1 - y);
        }

    }
}
