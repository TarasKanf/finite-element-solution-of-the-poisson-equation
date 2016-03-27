using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethodPE.Helpers;

namespace FiniteElementMethodPE.FiniteElements
{    
    class FiniteElementPuassonESolver
    {
        const int SIZE = 9;
        private double[,] A; // загальна матриця
        private double[,] b; // права частина
        private double[,] Ke = new double[SIZE, SIZE];
        private double[,] Me = new double[SIZE, SIZE];
        private double[] Qe = new double[SIZE];
        public FiniteElementPuassonESolver(Func<double,double,double> f, double a1,double b1, double a2,double b2)
        {
            // TODO 
            // приймаємо функ. f для правої частини рівняння
            // межі області [a1,b1]x[a2,b2]                     
        }
        public double[] Solve(int n)
        {
            //TODO
            // будуємо загальну матрицю A на основу малих 9х9 матриць Ke
            // будуємо праву частину на основі Qe
            // розвязуємо СЛАР
            return new double[n * n];
        }
        // нижче використовуйте клас Integral з Helpers щоб обчислити інтеграли
        private void FillKe(ref double[,] Ke, double x1a,double x1b,double x2a,double x2b) // розмірність матриці Ke 9x9
        {
            //TODO 
            // заповнюємо матрицю Ke для кокретної області [x1a,x1b]x[x2a,x2b]
        }
        private void FillQe(ref double[] Qe, double x1a, double x1b, double x2a, double x2b)
        {
            //TODO
            // заповнюємо вектор Qe для кокретної області [x1a,x1b]x[x2a,x2b]
        }
        private void FillMe(ref double[,] Me, double x1a, double x1b, double x2a, double x2b)
        {
            //TODO
            // заповнюємо матрицю Me для кокретної області [x1a,x1b]x[x2a,x2b]
        }

        // TODO
        // Прописати всі потрібні базисні фукнції l і їх похідні
        private double LagrangeFunc(double x,int i)
        {
            //TODO
            return 0;
        }
        private double DLanrangeFunc(double x,int i) // похідна функції Лагранжа
        {
            //TODO
            return 0;
        }

    }
}
