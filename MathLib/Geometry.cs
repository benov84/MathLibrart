using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

namespace Benov.MathLib
{
    public partial class Geometry
    {
        public static double MinDist(double X, double Y, double Z, int Ns, double[] Xs, double[] Ys, double[] Zs)
        {
            double Lpt = double.MaxValue;
            double Lper = double.MaxValue;
            for (int i = 2; i <= Ns; i++)
            {
                //Търсим перпендикуляр
                double Al = Core.PosAng(Xs[i - 1], Ys[i - 1], Xs[i], Ys[i]) + 100;
                double Xt = X + 100 * Math.Cos(Al / Core.R0);
                double Yt = Y + 100 * Math.Sin(Al / Core.R0);
                double Xp, Yp;
                int[] ip = Core.Prava(X, Y, Xt, Yt, Xs[i - 1], Ys[i - 1], Xs[i], Ys[i], out Xp, out Yp);
                if (ip[1] == 1)
                {
                    double Zp = Core.Kota(Xs[i - 1], Zs[i - 1], Xs[i], Zs[i], Xp);
                    double Ltemp = Core.Dist(X, Y, Xp, Yp);
                    if (Ltemp < Lper)
                        Lper = Ltemp;
                }

                //Проверяваме разстоянията до двете точки
                double Ltemp2 = Core.Dist(X, Y, Z, Xs[i - 1], Ys[i - 1], Zs[i - 1]);
                if (Ltemp2 < Lpt)
                    Lpt = Ltemp2;
                Ltemp2 = Core.Dist(X, Y, Z, Xs[i], Ys[i], Zs[i]);
                if (Ltemp2 < Lpt)
                    Lpt = Ltemp2;
            }
            if (Lper != 0)
                return Math.Min(Lpt, Lper);
            else
                return Lpt;
        }

        /*
        public static double FindDistanceToSegment(PointF pt, PointF p1, PointF p2, out PointF closest)
        {
            float dx = p2.X - p1.X;
            float dy = p2.Y - p1.Y;
            if ((dx == 0) && (dy == 0))
            {
                // It's a point not a line segment.
                closest = p1;
                dx = pt.X - p1.X;
                dy = pt.Y - p1.Y;
                return Math.Sqrt(dx * dx + dy * dy);
            }
            // Calculate the t that minimizes the distance.
            float t = ((pt.X - p1.X) * dx + (pt.Y - p1.Y) * dy) / (dx * dx + dy * dy);
            // See if this represents one of the segment'
            // end points or a point in the middle.
            if (t < 0)
            {
                closest = new PointF(p1.X, p1.Y);
                dx = pt.X - p1.X;
                dy = pt.Y - p1.Y;
            }
            else
                if (t > 1)
                {
                    closest = new PointF(p2.X, p2.Y);
                    dx = pt.X - p2.X;
                    dy = pt.Y - p2.Y;
                }
                else
                {
                    closest = new PointF(p1.X + t * dx, p1.Y + t * dy);
                    dx = pt.X - closest.X;
                    dy = pt.Y - closest.Y;
                }
            return Math.Sqrt(dx * dx + dy * dy);
        }
        */

        public static double Radius(double X1, double X2, double X3, double Y1, double Y2, double Y3, out double STR)
        {
            double[] X = new double[] { 0, X1, X2, X3 };
            double[] Y = new double[] { 0, Y1, Y2, Y3 };

            double[] XA = new double[3];
            double[] YA = new double[3];
            double[] XB = new double[3];
            double[] YB = new double[3];

            for (int i = 1; i <= 2; i++)
            {
                XA[i] = (X[i] + X[i + 1]) / 2;
                YA[i] = (Y[i] + Y[i + 1]) / 2;

                double AL = Core.PosAng(X[i], Y[i], X[i + 1], Y[i + 1]);

                XB[i] = (XA[i] + 100) * Math.Cos((AL + 100) / Core.R0);
                YB[i] = (YA[i] + 100) * Math.Sin((AL + 100) / Core.R0);
            }

            double XP, YP;
            int[] ind = Core.Prava(XA[1], YA[1], XB[1], YB[1], XA[2], YA[2], XB[2], YB[2], out XP, out YP);

            STR = Strelka(X1, Y1, X3, Y3, XP, YP);

            if (ind[1] + ind[2] + ind[3] > 0)
                return Core.Dist(XA[1], YA[1], XP, YP);

            return double.MaxValue;
        }

        public static double Strelka(double X1, double Y1, double X2, double Y2, double XC, double YC)
        {
            return Core.Dist(X1, Y1, XC, YC) - Core.Dist((X1 + X2) / 2, (Y1 + Y2) / 2, XC, YC);
        }

        #region Polyline Simpligication
        public static List<Point> DouglasPeuckerReduction(List<Point> Points, double Tolerance)
        {
            if (Points == null || Points.Count < 3)
                return Points;

            Int32 firstPoint = 0;
            Int32 lastPoint = Points.Count - 1;
            List<Int32> pointIndexsToKeep = new List<Int32>();

            //Add the first and last index to the keepers
            pointIndexsToKeep.Add(firstPoint);
            pointIndexsToKeep.Add(lastPoint);


            //The first and the last point can not be the same
            while (Points[firstPoint].Equals(Points[lastPoint]))
            {
                lastPoint--;
            }

            DouglasPeuckerReduction(Points, firstPoint, lastPoint, Tolerance, ref pointIndexsToKeep);

            List<Point> returnPoints = new List<Point>();
            pointIndexsToKeep.Sort();
            foreach (Int32 index in pointIndexsToKeep)
            {
                returnPoints.Add(Points[index]);
            }

            return returnPoints;
        }
        
        private static void DouglasPeuckerReduction(List<Point> points, Int32 firstPoint, Int32 lastPoint, Double tolerance, ref List<Int32> pointIndexsToKeep)
        {
            Double maxDistance = 0;
            Int32 indexFarthest = 0;

            for (Int32 index = firstPoint; index < lastPoint; index++)
            {
                Double distance = PerpendicularDistance(points[firstPoint], points[lastPoint], points[index]);
                if (distance > maxDistance)
                {
                    maxDistance = distance;
                    indexFarthest = index;
                }
            }

            if (maxDistance > tolerance && indexFarthest != 0)
            {
                //Add the largest point that exceeds the tolerance
                pointIndexsToKeep.Add(indexFarthest);

                DouglasPeuckerReduction(points, firstPoint, indexFarthest, tolerance, ref pointIndexsToKeep);
                DouglasPeuckerReduction(points, indexFarthest, lastPoint, tolerance, ref pointIndexsToKeep);
            }
        }
        #endregion

        public static Double PerpendicularDistance(Point Point1, Point Point2, Point Point)
        {
            //Area = |(1/2)(x1y2 + x2y3 + x3y1 - x2y1 - x3y2 - x1y3)|   *Area of triangle
            //Base = √((x1-x2)²+(x1-x2)²)                               *Base of Triangle*
            //Area = .5*Base*H                                          *Solve for height
            //Height = Area/.5/Base

            Double area = Math.Abs(.5 * (Point1.x * Point2.y + Point2.x * Point.y + Point.x * Point1.y - Point2.x * Point1.y - Point.x * Point2.y - Point1.x * Point.y));
            Double bottom = Math.Sqrt(Math.Pow(Point1.x - Point2.x, 2) + Math.Pow(Point1.y - Point2.y, 2));
            Double height = area / bottom * 2;

            return height;

            //Another option
            //Double A = Point.X - Point1.X;
            //Double B = Point.Y - Point1.Y;
            //Double C = Point2.X - Point1.X;
            //Double D = Point2.Y - Point1.Y;

            //Double dot = A * C + B * D;
            //Double len_sq = C * C + D * D;
            //Double param = dot / len_sq;

            //Double xx, yy;

            //if (param < 0)
            //{
            //    xx = Point1.X;
            //    yy = Point1.Y;
            //}
            //else if (param > 1)
            //{
            //    xx = Point2.X;
            //    yy = Point2.Y;
            //}
            //else
            //{
            //    xx = Point1.X + param * C;
            //    yy = Point1.Y + param * D;
            //}

            //Double d = DistanceBetweenOn2DPlane(Point, new Point(xx, yy));

        }

        /// <summary>
        /// Create a perpendicular offset point at a position located along a line segment.
        /// </summary>
        /// <param name="a">Input. PointD(x,y) of p1.</param>
        /// <param name="b">Input. PointD(x,y) of p2.</param>
        /// <param name="position">Distance between p1(0.0) and p2 (1.0) in a percentage.</param>
        /// <param name="offset">Distance from position at 90degrees to p1 and p2- non-percetange based.</param>
        /// <param name="c">Output of the calculated point along p1 and p2. might not be necessary for the ultimate output.</param>
        /// <param name="d">Output of the calculated offset point.</param>
        public static void PerpendicularOffset(Point a, Point b, double position, double offset, out Point c, out Point d, bool left)
        {
            //p3 is located at the x or y delta * position + p1x or p1y original.
            Point p3 = new Point(((b.x - a.x) * position) + a.x, ((b.y - a.y) * position) + a.y);

            //returns an angle in radians between p1 and p2 + 1.5708 (90degress).
            double angle = 90.0.DegreeToRadian();
            if (left)
                angle = -angle;
            double angleRadians = Math.Atan2(a.y - b.y, a.x - b.x) + angle;//1.5708;

            //locates p4 at the given angle and distance from p3.
            Point p4 = new Point(p3.x + Math.Cos(angleRadians) * offset, p3.y + Math.Sin(angleRadians) * offset);

            //send out the calculated points
            c = p3;
            d = p4;
        }

        public static double Azimuth(double x1, double y1, double x2, double y2)
        {
            double degBearing = Core.RadianToDegree(Math.Atan2((x2 - x1), (y2 - y1)));
            degBearing = Core.FitAngle(degBearing);
            return degBearing;
        }

        // Find the point of intersection between
        // the lines p1 --> p2 and p3 --> p4.
        public static void FindIntersection(
            PointF p1, PointF p2, PointF p3, PointF p4,
            out bool lines_intersect, out bool segments_intersect,
            out PointF intersection,
            out PointF close_p1, out PointF close_p2)
        {
            // Get the segments' parameters.
            float dx12 = p2.X - p1.X;
            float dy12 = p2.Y - p1.Y;
            float dx34 = p4.X - p3.X;
            float dy34 = p4.Y - p3.Y;

            // Solve for t1 and t2
            float denominator = (dy12 * dx34 - dx12 * dy34);

            float t1 =
                ((p1.X - p3.X) * dy34 + (p3.Y - p1.Y) * dx34)
                    / denominator;
            if (float.IsInfinity(t1))
            {
                // The lines are parallel (or close enough to it).
                lines_intersect = false;
                segments_intersect = false;
                intersection = new PointF(float.NaN, float.NaN);
                close_p1 = new PointF(float.NaN, float.NaN);
                close_p2 = new PointF(float.NaN, float.NaN);
                return;
            }
            lines_intersect = true;

            float t2 =
                ((p3.X - p1.X) * dy12 + (p1.Y - p3.Y) * dx12)
                    / -denominator;

            // Find the point of intersection.
            intersection = new PointF(p1.X + dx12 * t1, p1.Y + dy12 * t1);

            // The segments intersect if t1 and t2 are between 0 and 1.
            segments_intersect =
                ((t1 >= 0) && (t1 <= 1) &&
                 (t2 >= 0) && (t2 <= 1));

            // Find the closest points on the segments.
            if (t1 < 0)
            {
                t1 = 0;
            }
            else if (t1 > 1)
            {
                t1 = 1;
            }

            if (t2 < 0)
            {
                t2 = 0;
            }
            else if (t2 > 1)
            {
                t2 = 1;
            }

            close_p1 = new PointF(p1.X + dx12 * t1, p1.Y + dy12 * t1);
            close_p2 = new PointF(p3.X + dx34 * t2, p3.Y + dy34 * t2);
        }

        // Find a circle through the three points.
        public static void FindCircle(PointF a, PointF b, PointF c,
            out PointF center, out float radius)
        {
            // Get the perpendicular bisector of (x1, y1) and (x2, y2).
            float x1 = (b.X + a.X) / 2;
            float y1 = (b.Y + a.Y) / 2;
            float dy1 = b.X - a.X;
            float dx1 = -(b.Y - a.Y);

            // Get the perpendicular bisector of (x2, y2) and (x3, y3).
            float x2 = (c.X + b.X) / 2;
            float y2 = (c.Y + b.Y) / 2;
            float dy2 = c.X - b.X;
            float dx2 = -(c.Y - b.Y);

            // See where the lines intersect.
            bool lines_intersect, segments_intersect;
            PointF intersection, close1, close2;
            FindIntersection(
                new PointF(x1, y1), new PointF(x1 + dx1, y1 + dy1),
                new PointF(x2, y2), new PointF(x2 + dx2, y2 + dy2),
                out lines_intersect, out segments_intersect,
                out intersection, out close1, out close2);
            if (!lines_intersect)
            {
                center = new PointF(0, 0);
                radius = 0;
            }
            else
            {
                center = intersection;
                float dx = center.X - a.X;
                float dy = center.Y - a.Y;
                radius = (float)Math.Sqrt(dx * dx + dy * dy);
            }
        }
    }
}
