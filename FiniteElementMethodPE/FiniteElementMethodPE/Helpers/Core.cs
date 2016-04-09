using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethodPE.Helpers
{
    class Core:ICore
    {
        private Func<double, int,double, double> f1;
        private Func<double, int,double, double> f2;
        private int param1, param3;
        private double param2, param4;
        public Core(Func<double,int, double, double> _f1, Func<double,int, double, double> _f2)
        {
            f1 = _f1;
            f2 = _f2;
        }
        public double GetValue(double x)
        {
            return f1(x,param1,param2)*f1(x,param3, param4);
        }
        public void SetParams(int p1, double p2, int p3, double p4)
        {
            param1 = p1;
            param2 = p2;
            param3 = p3;
            param4 = p4;
        }

    }
}
