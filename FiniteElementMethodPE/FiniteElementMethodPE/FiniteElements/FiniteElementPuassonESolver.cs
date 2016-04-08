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
        const int Size = 9;
        private double[,] A; // загальна матриця
        private double[,] B; // права частина
        private double[,] Ke = new double[Size, Size];
        private double[,] Me = new double[Size, Size];
        private double[] Qe = new double[Size];
        private double h;
        private int N;
        private Func<double, double, double> f;
        // межі області
        private double a1, b1, a2, b2;

        //область має бути квадратна
        public FiniteElementPuassonESolver(Func<double,double,double> f, double a1,double b1, double a2,double b2)
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
            if(Math.Abs(h - (b2 -a2)/N) < 0.001) throw new Exception("Must be square");
            //TODO
            // будуємо загальну матрицю A на основу малих 9х9 матриць Ke
            // будуємо праву частину на основі Qe
            // розвязуємо СЛАР
            return new double[n * n];
        }
        // нижче використовуйте клас Integral з Helpers щоб обчислити інтеграли
        private void FillKe(ref double[,] Ke, double x1a, double x2a) // розмірність матриці Ke 9x9
        {
            //TODO 
            // заповнюємо матрицю Ke для кокретної області [x1a,x1b]x[x2a,x2b]
        }
        private void FillQe(ref double[] Qe, double x1a, double x2a)
        {
            //TODO
            // заповнюємо вектор Qe для кокретної області [x1a,x1b]x[x2a,x2b]
        }
        private void FillMe(ref double[,] Me, double x1a, double x2a)
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
        private double DLanrangeFunc(double x,int i,double xa) // похідна функції Лагранжа
        {
            double sum = 0;
            double dob = 1.0;
            if (x < xa || x > xa + h) return 0;
            for (var j = 0; j < 3; j++)
            {
                if(i == j) continue;
                sum += x - (xa + h*j/2.0);
                dob *= (i - j)*h/2.0;
            }
            return sum/dob;
        }

    }
}
