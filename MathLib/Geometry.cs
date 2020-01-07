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
            Point p1, Point p2, Point p3, Point p4,
            out bool lines_intersect, out bool segments_intersect,
            out Point intersection,
            out Point close_p1, out Point close_p2)
        {
            // Get the segments' parameters.
            double dx12 = p2.x - p1.x;
            double dy12 = p2.y - p1.y;
            double dx34 = p4.x - p3.x;
            double dy34 = p4.y - p3.y;

            // Solve for t1 and t2
            double denominator = (dy12 * dx34 - dx12 * dy34);

            double t1 =
                ((p1.x - p3.x) * dy34 + (p3.y - p1.y) * dx34)
                    / denominator;
            if (double.IsInfinity(t1))
            {
                // The lines are parallel (or close enough to it).
                lines_intersect = false;
                segments_intersect = false;
                intersection = new Point(float.NaN, float.NaN);
                close_p1 = new Point(float.NaN, float.NaN);
                close_p2 = new Point(float.NaN, float.NaN);
                return;
            }
            lines_intersect = true;

            double t2 =
                ((p3.x - p1.x) * dy12 + (p1.y - p3.y) * dx12)
                    / -denominator;

            // Find the point of intersection.
            intersection = new Point(p1.x + dx12 * t1, p1.y + dy12 * t1);

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

            close_p1 = new Point(p1.x + dx12 * t1, p1.y + dy12 * t1);
            close_p2 = new Point(p3.x + dx34 * t2, p3.y + dy34 * t2);
        }

        // Find a circle through the three points.
        public static void FindCircle(Point a, Point b, Point c,
            out Point center, out double radius)
        {
            // Get the perpendicular bisector of (x1, y1) and (x2, y2).
            double x1 = (b.x + a.x) / 2;
            double y1 = (b.y + a.y) / 2;
            double dy1 = b.x - a.x;
            double dx1 = -(b.y - a.y);

            // Get the perpendicular bisector of (x2, y2) and (x3, y3).
            double x2 = (c.x + b.x) / 2;
            double y2 = (c.y + b.y) / 2;
            double dy2 = c.x - b.x;
            double dx2 = -(c.y - b.y);

            // See where the lines intersect.
            bool lines_intersect, segments_intersect;
            Point intersection, close1, close2;
            FindIntersection(
                new Point(x1, y1), new Point(x1 + dx1, y1 + dy1),
                new Point(x2, y2), new Point(x2 + dx2, y2 + dy2),
                out lines_intersect, out segments_intersect,
                out intersection, out close1, out close2);
            if (!lines_intersect)
            {
                center = new Point(0, 0);
                radius = 0;
            }
            else
            {
                center = intersection;
                double dx = center.x - a.x;
                double dy = center.y - a.y;
                radius = Math.Sqrt(dx * dx + dy * dy);
            }
        }

        public static Point PointByDirection(Point A, Point B, double Distance)
        {
            double direction = Core.PosAngRadian(A.x, A.y, B.x, B.y);
            return new Point(A.x + Math.Cos(direction) * Distance, A.y + Math.Sin(direction) * Distance);
        }

        public static Point PointByDirection(Point A, double Direction, double Distance)
        {
            return new Point(A.x + Math.Cos(Direction) * Distance, A.y + Math.Sin(Direction) * Distance);
        }

        /*
        * @return integer code for which side of the line ab c is on.
        * 1 means left turn, -1 means right turn.  Returns
        * 0 if all three are on a line
        */
        public static int findSidePointLine(
            Point lineStart, Point lineEnd, Point testPoint)
        {
            if (lineEnd.x - lineStart.x == 0)
            { // vertical line
                if (testPoint.x < lineEnd.x)
                {
                    return lineEnd.y > lineStart.y ? 1 : -1;
                }
                if (testPoint.x > lineEnd.x)
                {
                    return lineEnd.y > lineStart.y ? -1 : 1;
                }
                return 0;
            }
            if (lineEnd.y - lineStart.y == 0)
            { // horizontal line
                if (testPoint.y < lineEnd.y)
                {
                    return lineEnd.x > lineStart.x ? -1 : 1;
                }
                if (testPoint.y > lineEnd.y)
                {
                    return lineEnd.x > lineStart.x ? 1 : -1;
                }
                return 0;
            }
            double slope = (lineEnd.y - lineStart.y) / (lineEnd.x - lineStart.x);
            double yIntercept = lineStart.y - lineStart.x * slope;
            double cSolution = (slope * testPoint.x) + yIntercept;
            if (slope != 0)
            {
                if (testPoint.y > cSolution)
                {
                    return lineEnd.x > lineStart.x ? 1 : -1;
                }
                if (testPoint.y < cSolution)
                {
                    return lineEnd.x > lineStart.x ? -1 : 1;
                }
                return 0;
            }
            return 0;
        }

        private static Point getPerpendicularPoint(Point lineStart, Point lineEnd, Point point)
        {
            var x1 = lineStart.x;
            var y1 = lineStart.y;
            var x2 = lineEnd.x;
            var y2 = lineEnd.y;
            var x3 = point.x;
            var y3 = point.y;
            var px = x2 - x1;
            var py = y2 - y1;
            var dAB = px * px + py * py;
            var u = ((x3 - x1) * px + (y3 - y1) * py) / dAB;
            var x = x1 + u * px;
            var y = y1 + u * py;
            return new Point(x, y);
        }

        public static bool PointOnLineSegment(Point pt1, Point pt2, Point pt, double epsilon = 0.001)
        {
            if (pt.x - Math.Max(pt1.x, pt2.x) > epsilon ||
                Math.Min(pt1.x, pt2.x) - pt.x > epsilon ||
                pt.y - Math.Max(pt1.y, pt2.y) > epsilon ||
                Math.Min(pt1.y, pt2.y) - pt.y > epsilon)
                return false;

            if (Math.Abs(pt2.x - pt1.x) < epsilon)
                return Math.Abs(pt1.x - pt.x) < epsilon || Math.Abs(pt2.x - pt.x) < epsilon;
            if (Math.Abs(pt2.y - pt1.y) < epsilon)
                return Math.Abs(pt1.y - pt.y) < epsilon || Math.Abs(pt2.y - pt.y) < epsilon;

            double x = pt1.x + (pt.y - pt1.y) * (pt2.x - pt1.x) / (pt2.y - pt1.y);
            double y = pt1.y + (pt.x - pt1.x) * (pt2.y - pt1.y) / (pt2.x - pt1.x);

            return Math.Abs(pt.x - x) < epsilon || Math.Abs(pt.y - y) < epsilon;
        }

        public static double DistanceBetweenPolylines(List<Point> line1, List<Point> line2, double precision = 5.0)
        {
            //Find the direction on first line
            double angle = Benov.MathLib.Core.PosAngRadian(line1.First(), line1.Last()) * (180.0 / Math.PI);
            int position = findSidePointLine(line1.First(), line1.Last(), line2.First());
            double angleLeft = angle + (position * 90);
            double angleRight = angle - (position * 90);
            //double step = Benov.MathLib.Core.Dist(line1.First().x, line1.First().y, line1.Last().x, line1.Last().y) / precision;
            double result = double.MaxValue;
            double maxDistance = Benov.MathLib.Core.Dist(line1.First().x, line1.First().y, line1.Last().x, line1.Last().y);
            double distance = 0;
            while (distance <= maxDistance)
            {
                Point testMiddle = PointByDirection(line1.First(), line1.Last(), distance);
                Point startPoint = PointByDirection(testMiddle, angleLeft * Math.PI / 180.0, maxDistance);
                Point endPoint = PointByDirection(testMiddle, angleRight * Math.PI / 180.0, maxDistance);

                Point IntersectionPoint1 = null;
                Point IntersectionPoint2 = null;

                for (int i = 1; i < line1.Count; i++) {
                    bool lines_intersect;
                    bool segment_intersect;
                    Point pt1, pt2;
                    FindIntersection(startPoint, endPoint, line1[i - 1], line1[i], 
                        out lines_intersect, out segment_intersect, out IntersectionPoint1, out pt1, out pt2);
                    if (lines_intersect)
                        break;
                }

                for (int i = 1; i < line2.Count; i++)
                {
                    bool lines_intersect;
                    bool segment_intersect;
                    Point pt1, pt2;
                    FindIntersection(startPoint, endPoint, line2[i - 1], line2[i],
                        out lines_intersect, out segment_intersect, out IntersectionPoint2, out pt1, out pt2);
                    if (lines_intersect)
                        break;
                }

                if (IntersectionPoint1 != null && IntersectionPoint2 != null &&
                    IntersectionPoint1.DistanceTo(IntersectionPoint2) < result)
                    result = IntersectionPoint1.DistanceTo(IntersectionPoint2);
            }

            return result;
        }
    }
}
