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
        private const int AxisPoints = 3;
        const int AllPointsNumber = AxisPoints*AxisPoints;
        private double[,] A; // загальна матриця
        private double[,] B; // права частина
        private double[,] Ke = new double[AllPointsNumber, AllPointsNumber];
        private double[,] Me = new double[AllPointsNumber, AllPointsNumber];
        private double[] Qe = new double[AllPointsNumber];
        private double h;
        private int N;
        private Func<double, double, double> f;
        // межі області
        private readonly double a1, b1, a2, b2;

        private double[] tempVector = new double[AllPointsNumber];

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
            // заповнюємо матрицю Ke для кокретної області [x1a,x1b]x[x2a,x2b]
            double E = 0.0005;
            for (int i = 0; i < AllPointsNumber; i++)
            {
                for (int j = 0; j < AllPointsNumber; j++)
                {
                    int sigma0 = i / 3;
                    int sigma1 = j / 3;
                    int sigma2 = i % 3;
                    int sigma3 = j % 3;
					Core core = new Core(LagrangeFunc, LagrangeFunc);
					Core core1 = new Core(DLanrangeFunc, DLanrangeFunc);
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
                    Ke[i, j] = temp;
                }
            }
        }

        private void FillQe(ref double[] Qe, double x1a, double x2a)
        {
            int index = 0;
            for (int i = 0; i < AxisPoints; i++)
            {
                for (int j = 0; j < AxisPoints; j++)
                {
                    tempVector[index] = f(x1a + j*h/2.0, x2a + i*h/2.0);
                    index ++;
                }
            }
            // множення матриці на вектор
            FillMe(ref Me,x1a,x2a);
            double sum = 0;
            for (int i = 0; i < AllPointsNumber; i++)
            {
                sum = 0;
                for (int j = 0; j < AllPointsNumber; j++)
                {
                    sum += Me[i, j]*tempVector[j];
                }
                Qe[i] = sum;
            }
        }

        private void FillMe(ref double[,] Me, double x1a, double x2a)
        {
            // заповнюємо матрицю Me для кокретної області [x1a,x1b]x[x2a,x2b]
            double E = 0.0005;
            for (int i = 0; i < AllPointsNumber; i++)
            {
                for (int j = 0; j < AllPointsNumber; j++)
                {
                    int sigma0 = i / 3;
                    int sigma1 = j / 3;
                    int sigma2 = i % 3;
                    int sigma3 = j % 3;
                    Core core = new Core(LagrangeFunc, LagrangeFunc);
                    core.SetParams(sigma0, x1a, sigma1, x1a);
                    Me[i, j] = Integral.CalculateWithHauseMethod(core, x1a, x1a + h, E, 5);
                    core.SetParams(sigma2, x2a, sigma3, x2a);
                    Me[i,j] *= Integral.CalculateWithHauseMethod(core, x2a, x2a + h, E, 5);
                }
            }
        }
                
        // Прописати всі потрібні базисні фукнції l і їх похідні
        private double LagrangeFunc(double x, int i, double xa)
        {
			double chus = 1.0;
			double dob = 1.0;
			if(x < xa || x > xa + h) return 0;
			for(var j = 0 ; j < 3 ; j++)
			{
				if(i == j) continue;
				chus *= x - (xa + h * j / 2.0);
				dob *= (i - j) * h / 2.0;
			}
			return chus / dob;
        }

        private double DLanrangeFunc(double x, int i, double xa) // похідна функції Лагранжа
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
