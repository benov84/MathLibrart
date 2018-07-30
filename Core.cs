using System;
using System.Text;

namespace Benov.MathLib
{
    public class Core
    {
        public static double R0 = 200 / Math.PI;

        enum UnitsAngle { Radian, Grads };
        enum UnitsDistance { Meter, Centimeter, Milimeter };

        public static bool IsOdd(int value)
        {
            bool a = value % 2 != 0;
            return a;
        }

        public static string NumberToRoman(int number)
        {
            // Validate
            if (number < 0 || number > 3999)
                throw new ArgumentException("Value must be between 0 - 3,999.");

            if (number == 0) return "N";

            // Set up key numerals and numeral pairs
            int[] values = new int[] { 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1 };
            string[] numerals = new string[] { "M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I" };

            // Initialise the string builder
            StringBuilder result = new StringBuilder();

            // Loop through each of the values to diminish the number
            for (int i = 0; i < 13; i++)
            {
                // If the number being converted is less than the test value, append
                // the corresponding numeral or numeral pair to the resultant string
                while (number >= values[i])
                {
                    number -= values[i];
                    result.Append(numerals[i]);
                }
            }

            // Done
            return result.ToString();
        }

        public static void SortArray(ref double[] Ar, out int n)
        {
            int[] SS = new int[Ar.Length + 1];
            double[] H1 = new double[Ar.Length + 1];

            for (int i = 0; i < Ar.Length; i++)
            {
                n = 0;
                for (int j = 0; j < Ar.Length; j++)
                    if (i != j && Ar[i] >= Ar[j])
                        n++;
                SS[i] = n;
            }

            for (int i = 0; i < Ar.Length; i++)
            {
                n = SS[i];
                H1[n] = Ar[i];
            }

            n = -1;
            for (int i = 0; i < Ar.Length; i++)
            {
                if (H1[i] != 0)
                {
                    n++;
                    Ar[n] = H1[i];
                }
            }
        }

        public static double Area(double[,] PL, int N)
        {
            double P = 0;

            PL[N + 1, 1] = PL[1, 1];
            PL[N + 1, 2] = PL[1, 2];

            for (int j = 1; j <= N; j++)
                P += ((PL[j, 1] + PL[j + 1, 1]) * (PL[j, 2] - PL[j + 1, 2]));
            return Math.Abs(P / 2);
        }

        public static double AreaHor(double[,] PL, int N, ref double[] PLH)
        {
            double[] HC = new double[202];
            for (int i = 1; i <= 201; i++)
                HC[i] = i * 15;

            PLH = new double[202];

            double P = 0;

            double Hmin, Hmax, Xmin, Xmax;
            Hmin = System.Double.MaxValue;
            Hmax = System.Double.MinValue;
            Xmin = System.Double.MaxValue;
            Xmax = System.Double.MinValue;

            for (int i = 1; i <= N; i++)
            {
                if (PL[i, 2] > Hmax) Hmax = PL[i, 2];
                if (PL[i, 2] < Hmin) Hmin = PL[i, 2];
                if (PL[i, 1] > Xmax) Xmax = PL[i, 1];
                if (PL[i, 1] < Xmin) Xmin = PL[i, 1];
            }

            double d = 0.05;

            int Nmin = Convert.ToInt32(Math.Truncate(Hmin / 15)) - 1;
            int Nmax = Convert.ToInt32(Math.Truncate(Hmax / 15)) + 1;
            Xmin = (Math.Round((Xmin * 100) / 5) * 5) / 100;
            Xmax = (Math.Round((Xmax * 100) / 5) * 5) / 100;

            PL[N + 1, 1] = PL[1, 1];
            PL[N + 1, 2] = PL[1, 2];

            for (int i = Nmin; i <= Nmax; i++)
            {
                PLH[i] = 0;
                int m = Convert.ToInt32(15 / d);
                for (int k = 1; k <= m; k++)
                {
                    double HT = i * 15 + (k - 1) * d + d / 2;
                    int q = Convert.ToInt32((Xmax - Xmin) / d);
                    for (int j = 1; j <= q; j++)
                    {
                        double XT = Xmin + (j - 1) * d + d / 2;
                        if (Contour(XT, HT, Xmin - 100, Hmin - 100, N, PL) == 1) PLH[i] += d * d;
                    }
                }
            }

            for (int j = 1; j <= N; j++)
                P += ((PL[j, 1] + PL[j + 1, 1]) * (PL[j, 2] - PL[j + 1, 2]));
            return Math.Abs(P / 2);
        }

        public static int Contour(double XT, double YT, double Xmin, double Ymin, int Nk, double[,] XYHK)
        {
            double XP, YP, XA, YA, XB, YB;
            int Ibr;
            int[] ip;
            ip = new int[4];

            for (int i = 1; i <= Nk; i++) if ((Math.Abs(XYHK[i, 1] - XT) < 0.07) & (Math.Abs(XYHK[i, 2] - YT) < 0.07)) return 1;

            Ibr = 0;
            for (int i = 1; i <= Nk; i++)
            {
                XA = XYHK[i, 1];
                YA = XYHK[i, 2];
                XB = XYHK[i + 1, 1];
                YB = XYHK[i + 1, 2];

                XP = 0;
                YP = 0;
                ip = Prava(XA, YA, XB, YB, XT, YT, Xmin, Ymin, out XP, out YP);

                if ((ip[2] == 1) & (XA == XP) & (YA == YP))
                {
                    Ibr += 1;
                    goto Zvezda;
                }

                if (ip[1] + ip[2] == 2)
                {
                    Ibr += 1;
                    goto Zvezda;
                }

                Zvezda:;
            }
            if ((Ibr & 1) == 1)
                return 1;
            else
                return 2;
        }

        public static double Dist(double XA, double YA, double XB, double YB)
        {
            double DX, DY;
            DX = XB - XA;
            DY = YB - YA;
            return Math.Sqrt(DX * DX + DY * DY);
        }

        public static double Dist(double XA, double YA, double ZA, double XB, double YB, double ZB)
        {
            double DX, DY, DZ;
            DX = XB - XA;
            DY = YB - YA;
            DZ = ZB - ZA;
            return Math.Sqrt(DX * DX + DY * DY + DZ * DZ);
        }

        public static int[] Prava(double XA, double YA, double XB, double YB, double XC, double YC, double XD, double YD, out double XP, out double YP)
        {
            double T1, T2;
            XP = 0;
            YP = 0;
            int[] PravaRes;
            PravaRes = new int[4];

            if ((XA - XB) * (XD - XC) != 0)
            {
                T1 = (YB - YA) / (XB - XA);
                T2 = (YD - YC) / (XD - XC);
                if (T1 - T2 == 0) return PravaRes;
                else PravaRes[3] = 1;

                XP = (YC - YA + T1 * XA - T2 * XC) / (T1 - T2);
                YP = T1 * (XP - XA) + YA;
            }

            if ((XB - XA) == 0 && (XD - XC) != 0)
            {
                T2 = (YD - YC) / (XD - XC);
                XP = XA;
                YP = T2 * (XP - XC) + YC;
                PravaRes[3] = 1;
            }

            if ((XB - XA) != 0 && (XD - XC) == 0)
            {
                T1 = (YB - YA) / (XB - XA);
                XP = XC;
                YP = T1 * (XP - XA) + YA;
                PravaRes[3] = 1;
            }

            if ((XB - XA) == 0 && (XD - XC) == 0)
            {
                PravaRes[3] = 0;
                return PravaRes;
            }

            if ((XA - XP) * (XB - XP) <= 0 && (YA - YP) * (YB - YP) <= 0)
                PravaRes[1] = 1;
            if ((XC - XP) * (XD - XP) <= 0 && (YC - YP) * (YD - YP) <= 0)
                PravaRes[2] = 1;
            return PravaRes;
        }

        public static bool PointOnLine(double XL1, double YL1, double XL2, double YL2, double XP, double YP)
        {
            bool OnLine = Math.Round((YP - YL1) / (YL2 - YL1), 2) == Math.Round((XP - XL1) / (XL2 - XL1), 2);
            if (OnLine)
            {
                if (Math.Round(Dist(XL1, YL1, XP, YP) + Dist(XL2, YL2, XP, YP), 2) > Math.Round(Dist(XL1, YL1, XL2, YL2), 2))
                    return false;
                else
                    return true;
            }
            else
                return false;
        }

        public static double PosAng(double X1, double Y1, double X2, double Y2)
        {
            double ta, a, b, a1, ma;

            a = X2 - X1;
            b = Y2 - Y1;
            a1 = Math.Abs(a);

            if (a1 > 0.001)
            {
                ta = Math.Atan(b / a);
                if (a < 0) { ta += Math.PI; };
                if (ta < 0) { ta += Math.PI * 2; };
                ma = ta;
            }
            else
                if (b == 0)
            {
                ma = 0;
            }
            else
                    if (b > 0)
            {
                ma = Math.PI / 2;
            }
            else
                ma = 3 * Math.PI / 2;

            return ma * R0;
        }

        public static double PosAngRadian(double X1, double Y1, double X2, double Y2)
        {
            double ta, a, b, a1, ma;

            a = X2 - X1;
            b = Y2 - Y1;
            a1 = Math.Abs(a);

            if (a1 > 0.001)
            {
                ta = Math.Atan(b / a);
                if (a < 0) { ta += Math.PI; };
                if (ta < 0) { ta += Math.PI * 2; };
                ma = ta;
            }
            else
                if (b == 0)
            {
                ma = 0;
            }
            else
                    if (b > 0)
            {
                ma = Math.PI / 2;
            }
            else
                ma = 3 * Math.PI / 2;

            return ma;
        }

        public static double OrientUnknown(double XA, double YA, double[] XB, double[] YB, double[] RB, int NB, out double Gr)
        {
            double[] Orn = new double[10];
            double DX, DY, A, Su, V;
            Su = 0;
            double Orient = 0;
            Gr = 0;

            for (int m = 1; m <= NB; m++)
            {
                A = PosAng(XA, YA, XB[m], YB[m]);
                DY = YB[m] - YA;
                DX = XB[m] - XA;
                if (DX * DY != 0)
                    A = Math.Atan(DY / DX) * R0 + (1 - DX / Math.Abs(DX)) * 100 +
                        (1 + DX / Math.Abs(DX)) * (1 - DY / Math.Abs(DY)) * 100;
                else
                {
                    A = 0;
                    if (DX + DY == 0) return 0;
                    if (DX == 0) A = 100 + (1 - DY / Math.Abs(DY)) * 100;
                    if (DY == 0) A = (1 - DX / Math.Abs(DX)) * 100;
                }
                Orn[m] = A - RB[m];
                if (Orn[m] < 0) Orn[m] += 400;
            }

            for (int k = 1; k <= NB; k++)
                Su += Orn[k];

            Orient = Su / NB;

            Su = 0;

            for (int k = 1; k <= NB; k++)
            {
                V = Orient - Orn[k];
                Su += V * V;
            }

            if (NB > 1) Su = Math.Sqrt(Su / (NB - 1));

            Gr = Su;
            return Orient;
        }

        public static double OrientUnknown(double XA, double YA, double[] XB, double[] YB, double[] RB, int NB, bool[] Incl, out double Gr)
        {
            double[] Orn = new double[10];
            double DX, DY, A, Su, V;
            Su = 0;
            double Orient = 0;
            Gr = 0;

            double[] _XB = new double[100];
            double[] _YB = new double[100];
            double[] _RB = new double[100];
            int _NB = 0;
            for (int i = 1; i <= NB; i++)
            {
                if (Incl[i])
                {
                    _NB++;
                    _XB[i] = XB[i];
                    _YB[i] = YB[i];
                    _RB[i] = RB[i];
                }
            }

            for (int m = 1; m <= _NB; m++)
            {
                A = PosAng(XA, YA, _XB[m], _YB[m]);
                DY = _YB[m] - YA;
                DX = _XB[m] - XA;
                if (DX * DY != 0)
                    A = Math.Atan(DY / DX) * R0 + (1 - DX / Math.Abs(DX)) * 100 +
                        (1 + DX / Math.Abs(DX)) * (1 - DY / Math.Abs(DY)) * 100;
                else
                {
                    A = 0;
                    if (DX + DY == 0)
                        return 0;
                    if (DX == 0)
                        A = 100 + (1 - DY / Math.Abs(DY)) * 100;
                    if (DY == 0)
                        A = (1 - DX / Math.Abs(DX)) * 100;
                }
                Orn[m] = A - _RB[m];
                if (Orn[m] < 0) Orn[m] += 400;
            }

            for (int k = 1; k <= _NB; k++)
                Su += Orn[k];

            Orient = Su / _NB;

            Su = 0;

            for (int k = 1; k <= _NB; k++)
            {
                V = Orient - Orn[k];
                Su += V * V;
            }

            if (_NB > 1)
                Su = Math.Sqrt(Su / (_NB - 1));

            Gr = Su;
            return Orient;
        }

        public static double Kota(double XA, double HA, double XB, double HB, double XT)
        {
            if (XB - XA != 0)
                return (XT - XB) * (HB - HA) / (XB - XA) + HB;
            else
                return (HA + HB) / 2;
        }

        public static double GroupUpperLimit(double min, double max, double step, double number)
        {
            if (number < min)
                return min;
            if (number > max)
                return max;

            for (double i = min; i <= max; i += step)
            {
                if (number < i)
                    return i;
            }

            return max;
        }

        public static double DegreeToRadian(double angle)
        {
            return Math.PI * angle / 180.0;
        }

        public static double RadianToDegree(double angle)
        {
            return angle * (180.0 / Math.PI);
        }

        public static int GetQuadrant(double XA, double YA, double XB, double YB)
        {
            if (XA <= XB && YA <= YB) return 1;
            if (XA <= XB && YA > YB) return 2;
            if (XA > XB && YA >= YB) return 3;
            if (XA > XB && YA < YB) return 4;
            return 0;
        }

        public static double FitAngle(double Angle)
        {
            while (Angle > 360)
                Angle -= 360;
            while (Angle < 0)
                Angle += 360;
            return Angle;
        }

        public static double RoundToFraction(double Number, double Fraction)
        {
            return Math.Round(Number / Fraction) * Fraction;
        }
    }
}
