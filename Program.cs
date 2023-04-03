using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ScottPlot;
using ScottPlot.Drawing.Colormaps;


class Laguerre
{
    private double beta;
    private double sigma;
    private int numberOfPoints { get; set; }
    private double eps;
    private double T;
    private int N;
    private Func<double, double> f;

    public Laguerre(double beta, double sigma, int numberOfPoints, double eps, double T, int n, int N, Func<double, double> f)
    {
        this.beta = beta;
        this.sigma = sigma;
        this.numberOfPoints = numberOfPoints;
        this.eps = eps;
        this.T = T;
        this.N = N;
        this.f = f;
    }

    public static double[] GetIntervals(double T, int numberOfPoints)
    {
        double[] indices = Enumerable.Range(0, numberOfPoints)
            .Select(i => i * T / (numberOfPoints - 1))
            .ToArray();
        return indices;
    }

    public double LaguerreFunc(double t, int n)
    {
        double laguerreResult0 = Math.Sqrt(sigma) * (Math.Exp(-(beta * t) / 2));
        double laguerreResult1 = laguerreResult0 * (1 - (sigma * t));

        if (n == 0)
        {
            return laguerreResult0;
        }
        else if (n == 1)
        {
            return laguerreResult1;
        }
        else
        {
            double laguerreResultN = 0;
            for (int i = 2; i <= n; i++)
            {
                laguerreResultN = ((2 * i - 1 - t * sigma) / i) * laguerreResult1 - ((i - 1) * laguerreResult0) / i;
                laguerreResult0 = laguerreResult1;
                laguerreResult1 = laguerreResultN;
            }
            return laguerreResultN;
        }
    }

    public double[] TabulateLaguerre(double T, int numberOfPoints, int n)
    {
        double[] indices = Laguerre.GetIntervals(T, numberOfPoints);
        double[] results = new double[numberOfPoints];

        for (int i = 0; i < numberOfPoints; i++)
        {
            double t = indices[i];
            results[i] = LaguerreFunc(t, n);
        }

        return results;
    }

    public void PlotTabulated(double[] x, double[] y)
    {
        var plt = new Plot(1000, 600);

        plt.AddScatter(x, y, color: System.Drawing.Color.Blue, lineWidth: 2);

        plt.Title("Tabulated function");

        plt.XLabel("t");
        plt.YLabel("laguerre");

        plt.SaveFig("tabulated_function.png");
    }

    public double[] RectangleIntegral(Func<double, double> f, double[] points, double delta)
    {
        double alpha = sigma - beta;
        double[] result = new double[N + 1];
        //using (StreamWriter writer = new StreamWriter(@"../../../rect.txt"))
        //{
        //    writer.WriteLine("x,sum");
        //    double sum = 0;
        //    for (int i = 0; i < points.Length; i++)
        //    {
        //        sum += f(points[i]) * LaguerreFunc(points[i], 0) * Math.Exp(-alpha * points[i]) * delta;
        //        writer.WriteLine($"{(points[i]).ToString("0.000")},{sum.ToString("0.000000")}");
        //    }
        //}
        for (int k = 0; k <= N; k++)
        {
            double sum = 0;
            for (int i = 0; i < points.Length; i++)
            {

                sum += f(points[i]) * LaguerreFunc(points[i], k) * Math.Exp(-alpha * points[i]) * delta;
            }
            result[k] = sum;
        }
        return result;
    }

    public double[] Transformation(Func<double, double> func)
    {
        int number_of_points = numberOfPoints;
        double delt = T / (number_of_points - 1);
        double half_delta = delt / 2;
        double[] points = Enumerable.Range(0, number_of_points - 1)
        .Select(i => 0 + half_delta + i * delt)
        .ToArray();
        double[] res_0 = new double[N + 1];
        double[] res_1 = RectangleIntegral(func, points, delt);

        while (Enumerable.Max(res_0.Select((x, i) => Math.Abs(x - res_1[i]))) > eps)
        {
            res_0 = res_1;
            number_of_points *= 2;
            double delta_ = T / (number_of_points - 1);
            double half_delta_ = delta_ / 2;
            double[] points_ = Enumerable.Range(0, number_of_points - 1)
                .Select(i => 0 + half_delta + i * delta_)
                .ToArray();
            res_1 = RectangleIntegral(func, points_, delta_);
        }

        return res_1;

    }
    public double InverseTransformation(int t, Func<double, double> func)
    {
        double[] seq = Transformation(func);
        int len = seq.Length;
        double[] k = Enumerable.Range(0, len).Select(x => (double)x).ToArray();
        double[] lag = Enumerable.Range(0, len).Select(x => LaguerreFunc(x, t)).ToArray();
        double h = seq.Zip(lag, (s, l) => s * l).Sum();
        return h;
    }
}

class Program
{
    static double Func(double t)
    {
        if (t >= 0 && t <= 2 * Math.PI)
        {
            return Math.Sin(t - Math.PI / 2) + 1;
        }
        else if (t > 2 * Math.PI)
        {
            return 0;
        }
        return 0;
    }

    static void Main(string[] args)
    {
        double beta = 2;
        double sigma = 4;
        int n = 5;
        int t = 7;
        double T = 7;
        double eps = 0.001;
        int N = 20;
        int numberOfPoints = 100;

        Laguerre laguerre = new Laguerre(beta, sigma, numberOfPoints, eps, T, n, N, Func);

        double laguerreValue = laguerre.LaguerreFunc(T, n);

        Console.WriteLine($"Laguerre({t}, {n}) = {laguerreValue:f2}\n");

        double[] tabulatedValues = laguerre.TabulateLaguerre(T, numberOfPoints, n);

        Console.WriteLine("Tabulated Laguerre");

        using (StreamWriter writer = new StreamWriter(@"../../../laguerre.txt"))
        {
            writer.WriteLine("x,y");
            for (int i = 0; i < numberOfPoints; i++)
            {
                Console.WriteLine($"t = {(i * T / (numberOfPoints - 1)).ToString("0.00")}, Laguerre(t) = {tabulatedValues[i].ToString("0.00")}");
                writer.WriteLine($"{(i * T / (numberOfPoints - 1)).ToString("0.00")},{tabulatedValues[i].ToString("0.00")}");
            }
        }

        Console.WriteLine("\n\n\nData written to file successfully.");

        Console.WriteLine("\n\n");

        laguerre.PlotTabulated(Laguerre.GetIntervals(T, numberOfPoints), tabulatedValues);

        double[] points = Laguerre.GetIntervals(T, numberOfPoints);
        double delta = T / (numberOfPoints - 1);
        double[] result1 = laguerre.RectangleIntegral(Func, points, delta);
        Console.WriteLine($"Integration result: {result1[0]:F2}\n");

        double[] result = laguerre.Transformation(Func);
        Console.WriteLine($"Integration result: {result[0]:F2}\n");

        double inverse_result = laguerre.InverseTransformation(t, Func);
        Console.WriteLine($"Integration result: {inverse_result:F2}");
    }

}