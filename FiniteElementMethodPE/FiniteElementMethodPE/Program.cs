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
            //ToDo
            // тестуємо наші класи
        }

        // f при u* = xy^2(1-x)(1-y)
        public static double f(double x, double y)
        {
            return 2 * y * y * y - 2 * y * y + 2 * x + 6 * x * x * y - 2 * x * x - 6 * x * y;
        }

    }
}
