using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethodPE.Helpers
{
    class Core:ICore
    {
        private Func<double, int, double> f1;
        private Func<double, int, double> f2;
        private int param1, param2;
        public Core(Func<double,int,double> _f1, Func<double,int,double> _f2)
        {
            f1 = _f1;
            f2 = _f2;
        }
        public double GetValue(double x)
        {
            return f1(x,param1)*f1(x,param2);
        }
        public void SetParams(int p1, int p2)
        {
            param1 = p1;
            param2 = p2;
        }

    }
}
