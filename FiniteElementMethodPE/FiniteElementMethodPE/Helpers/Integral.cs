using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethodPE.Helpers
{
    // обчислює інтеграли різними методами від заданих підінтегральних функцій
    // або від занаданого класу ще реалізує інтерфейc ICore
    static class Integral
    {
        private struct HauseKoefitient
        {
            public  double argum;
            public  double koef;
            public HauseKoefitient(double _argum,double _koef)
            {
                argum = _argum;
                koef= _koef;
            }
        }
        private static List<HauseKoefitient> masN3; 
        private static List<HauseKoefitient> masN4;
        private static List<HauseKoefitient> masN5;        
        static Integral()
        {
            masN5 = new List<HauseKoefitient>() { new HauseKoefitient(-0.9061798, 0.2369269), new HauseKoefitient(-0.5384693, 0.4786287), new HauseKoefitient(0, 0.5688889), new HauseKoefitient(0.5384693, 0.4786287), new HauseKoefitient(0.9061798, 0.2369269) };
            masN4 = new List<HauseKoefitient>() { new HauseKoefitient(-0.8611363, 0.3478548), new HauseKoefitient(-0.3399810, 0.6521452), new HauseKoefitient(0.3399810, 0.6521452), new HauseKoefitient(0.8611363, 0.3478548) };    
            masN3 = new List<HauseKoefitient>() { new HauseKoefitient(-0.7745967, 0.5555556),new HauseKoefitient(0, 0.8888889), new HauseKoefitient(0.7745967, 0.5555556) };
        }

        public static double CalculateWithSimpsonMethod(Func<double,double> f, double a, double b, double eps)
        {            
            if (a == b) return 0;
            double I2n = 0, In = 0;
            eps = Math.Abs(eps);
            double rizn = 1;
            double h = 0;
            double koef;
            In = 0;
            h = (b - a);
            In += (h / 6) * (f(a) + 4.0 * f((a + b) / 2) + f(b));
            for (int n = 4; (n < 10000) && (rizn > eps); n *= 2)
            {
                I2n = 0;
                h = (b- a) / n;
                I2n += f(a) + f(b);
                bool ans = true;
                for (int i = 1; i < n; i++)
                {
                    koef = (ans) ? 4.0 : 2.0;
                    I2n += koef* f(a + i * h);
                    ans = !ans;
                }                    
                I2n *= h/3.0;
                rizn = Math.Abs((I2n - In) / I2n);
                In = I2n;
            }            
            return I2n;
        }
        public static double CalculateWithHauseMethod(Func<double,double> f, double a, double b, double eps, int nHause)
        {           
            if (a == b) return 0;
            if ((nHause < 3) || (nHause > 5)) return 0;
            List<HauseKoefitient> hausekoef = null;
            switch (nHause)
            {
                case 3:
                    {
                        hausekoef = masN3;
                        break;
                    }
                case 4:
                    {
                        hausekoef = masN4;
                        break;
                    }
                case 5:
                    {
                        hausekoef = masN5;
                        break;
                    }
            }
            if (hausekoef == null) return 0;
            double h = b - a;
            double I2n = 0, In = 0;
            eps = Math.Abs(eps);
            double rizn = 1;         
            for (int i = 0; i < nHause;i++ )
            {
                In += hausekoef[i].koef * f((a + b) / 2 + (b - a) * hausekoef[i].argum/2.0);
            }
            In *= (b - a) / 2;                       
            for (int n = 4; (n < 10000) && (rizn > eps); n *= 2)
                {
                    I2n = 0;
                    h = (b - a) / n;
                    double a1 = a;
                    for (int j = 0; j < n;j++ )
                    {
                        for (int i = 0; i < nHause; i++)
                        {
                            I2n += hausekoef[i].koef * f((2*a1 + h) / 2.0 + (h) * hausekoef[i].argum / 2);
                        }                        
                        a1 += h;
                    }
                    I2n *= h / 2;
                    rizn = Math.Abs((I2n - In) / I2n);
                    In = I2n;
                }            
            return I2n;
        }
        public static double CalculateWithHauseMethod(ICore f,double a,double b,double eps,int nHause)
        {
            if (a == b) return 0;
            if ((nHause < 3) || (nHause > 5)) return 0;
            List<HauseKoefitient> hausekoef = null;
            switch (nHause)
            {
                case 3:
                    {
                        hausekoef = masN3;
                        break;
                    }
                case 4:
                    {
                        hausekoef = masN4;
                        break;
                    }
                case 5:
                    {
                        hausekoef = masN5;
                        break;
                    }
            }
            if (hausekoef == null) return 0;
            double h = b - a;
            double I2n = 0, In = 0;
            eps = Math.Abs(eps);
            double rizn = 1;
            for (int i = 0; i < nHause; i++)
            {
                In += hausekoef[i].koef * f.GetValue((a + b) / 2 + (b - a) * hausekoef[i].argum / 2.0);
            }
            In *= (b - a) / 2;
            for (int n = 4; (n < 10000) && (rizn > eps); n *= 2)
            {                
                I2n = 0;
                h = (b - a) / n;
                double a1 = a;
                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < nHause; i++)
                    {
                        I2n += hausekoef[i].koef * f.GetValue((2 * a1 + h) / 2.0 + (h) * hausekoef[i].argum / 2);
                    }
                    a1 += h;
                }
                I2n *= h / 2;
                rizn = Math.Abs((I2n - In) / I2n);
                In = I2n;
            }
            return I2n;
        }
    }
}
