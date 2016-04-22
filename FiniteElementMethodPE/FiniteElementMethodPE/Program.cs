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
            Console.Write(" x: ");
            double x = int.Parse(Console.ReadLine());
            Console.Write(" y: ");
            double y = int.Parse(Console.ReadLine());
            Console.Write(" H : ");
            double h = int.Parse(Console.ReadLine());
            Console.Write(" N : ");
            int n = int.Parse(Console.ReadLine());
            FiniteElementPuassonESolver solver = new FiniteElementPuassonESolver(F, x, x+h, y, y+ h);
            double[] result = solver.Solve(n);

            for (int i = 0, index = 0; i < 2*n + 1; i++)
            {
                for (int j = 0; j < 2*n + 1; j++)
                {
                    Console.Write("{0:F8}  ", result[index]);
                    index++;
                }
                Console.WriteLine();
            }
            Console.ReadKey();

        }

        // f при u* = xy^2(1-x)(1-y)
        public static double f(double x, double y)
        {
            return 2 * y * y * y - 2 * y * y + 2 * x + 6 * x * x * y - 2 * x * x - 6 * x * y;
        }

        // f при u* = xy^2(1-x)^2(1-y)
        public static double F(double x, double y)
        {
            return -2 * Math.Pow(1 - x,2) * x * (1 - y) 
                +4 * Math.Pow(1 - x,2) * x * y 
                +4 * (1 - x) * (1 - y) * y * y 
                - 2 * x * (1 - y) * y *y;
        }

        public static double AccurateSolution(double x, double y)
        {
            return x*y*y*Math.Pow(1 - x,2)*(1 - y);
        }

    }
}
