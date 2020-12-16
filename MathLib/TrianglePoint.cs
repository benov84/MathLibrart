using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace Benov.MathLib
{
    public enum RelPos2D
    {
        ll = 1,
        le = 2,
        lg = 3,
        eg = 4,
        gg = 5,
        ge = 6,
        gl = 7,
        el = 8,
        ee = 0
    }

    public static class TrianglePointTools
    {
        public static double Distance(Point3D Point1, Point3D Point2)
        {
            return Math.Sqrt(Math.Pow(Point1.X - Point2.X, 2) + Math.Pow(Point1.Y - Point2.Y, 2));
        }

        public static RelPos2D RelativePosition(Point3D Of, Point3D To)
        {
            int xRel = Of.X < To.X ? -1 : Of.X > To.X ? 1 : 0;
            int yRel = Of.Y < To.Y ? -1 : Of.Y > To.Y ? 1 : 0;

            switch (xRel)
            {
                case -1:
                    switch (yRel)
                    {
                        case -1: return RelPos2D.ll;
                        case 0: return RelPos2D.le;
                        case 1: return RelPos2D.lg;
                    }
                    break;
                case 0:
                    switch (yRel)
                    {
                        case -1: return RelPos2D.el;
                        case 0: return RelPos2D.ee;
                        case 1: return RelPos2D.eg;
                    }
                    break;
                case 1:
                    switch (yRel)
                    {
                        case -1: return RelPos2D.gl;
                        case 0: return RelPos2D.ge;
                        case 1: return RelPos2D.gg;
                    }
                    break;
            }

            return RelPos2D.ee; // never reached
        }

        public static double TriangleArea(Point3D Point1, Point3D Point2, Point3D Point3)
        {
            return 1 / 2d *
                (
                    (Point1.X - Point3.X) * (Point2.Y - Point1.Y) -
                    (Point1.X - Point2.X) * (Point3.Y - Point1.Y)
                );
        }

        public static bool TriangleContainsPoint(Point3D Point1, Point3D Point2, Point3D Point3, Point3D Target)
        {
            var s = Point1.Y * Point3.X - Point1.X * Point3.Y + (Point3.Y - Point1.Y) * Target.X + (Point1.X - Point3.X) * Target.Y;
            var t = Point1.X * Point2.Y - Point1.Y * Point2.X + (Point1.Y - Point2.Y) * Target.X + (Point2.X - Point1.X) * Target.Y;

            if ((s < 0) != (t < 0))
                return false;

            var area = TriangleArea(Point1, Point2, Point3);
            var sign = area < 0 ? -1 : 1;
            s *= sign;
            t *= sign;
            area *= sign;

            return s > 0 && t > 0 && (s + t) < 2 * area;
        }
    }


    public class TrianglePointProblemSolver
    {
        private static RelPos2D[] AllPositions = new RelPos2D[]
        {
            RelPos2D.ee,
            RelPos2D.eg,
            RelPos2D.el,
            RelPos2D.ge,
            RelPos2D.gg,
            RelPos2D.gl,
            RelPos2D.le,
            RelPos2D.lg,
            RelPos2D.ll,
        };

        private static RelPos2D[] NoPositions = new RelPos2D[0];

        private static RelPos2D[] ValidPositions(RelPos2D Pos1, RelPos2D Pos2)
        {
            if (Pos1 == RelPos2D.ee || Pos2 == RelPos2D.ee)
                return AllPositions;

            switch (Pos1)
            {
                case RelPos2D.ll:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg };
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge };
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl };
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.gg:
                            return AllPositions;
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg };
                    }
                    break;
                case RelPos2D.le:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge };
                        case RelPos2D.le:
                            return NoPositions;
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl };
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.ge:
                            return AllPositions.Except(new RelPos2D[] { Pos1, Pos2 }).ToArray();
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg };
                    }
                    break;
                case RelPos2D.lg:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl };
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl };
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl };
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.gl:
                            return AllPositions;
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl };
                    }
                    break;
                case RelPos2D.eg:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el };
                        case RelPos2D.eg:
                            return NoPositions;
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.el:
                            return AllPositions.Except(new RelPos2D[] { Pos1, Pos2 }).ToArray();
                    }
                    break;
                case RelPos2D.gg:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return AllPositions;
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ge, RelPos2D.gl, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll };
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll };
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le, RelPos2D.lg, RelPos2D.eg };
                    }
                    break;
                case RelPos2D.ge:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.le:
                            return AllPositions.Except(new RelPos2D[] { Pos1, Pos2 }).ToArray();
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.gl, RelPos2D.el, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le };
                        case RelPos2D.ge:
                            return NoPositions;
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg, RelPos2D.eg };
                    }
                    break;
                case RelPos2D.gl:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg, RelPos2D.gg, RelPos2D.ge };
                        case RelPos2D.lg:
                            return AllPositions;
                        case RelPos2D.eg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.el, RelPos2D.ll, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg };
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg };
                        case RelPos2D.el:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg };
                    }
                    break;
                case RelPos2D.el:
                    switch (Pos2)
                    {
                        case RelPos2D.ll:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg };
                        case RelPos2D.le:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg, RelPos2D.ge };
                        case RelPos2D.lg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.eg, RelPos2D.gg, RelPos2D.ge, RelPos2D.gl };
                        case RelPos2D.eg:
                            return AllPositions.Except(new RelPos2D[] { Pos1, Pos2 }).ToArray();
                        case RelPos2D.gg:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.ll, RelPos2D.le, RelPos2D.lg, RelPos2D.eg };
                        case RelPos2D.ge:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.le, RelPos2D.lg, RelPos2D.eg };
                        case RelPos2D.gl:
                            return new RelPos2D[] { RelPos2D.ee, RelPos2D.lg, RelPos2D.eg };
                        case RelPos2D.el:
                            return NoPositions;
                    }
                    break;
            }
            return NoPositions;
        }

        private static bool UncertainSolution(RelPos2D Pos1, RelPos2D Pos2, RelPos2D Pos3)
        {
            RelPos2D[] array = new RelPos2D[] { Pos1, Pos2, Pos3 };

            return
                (array.Contains(RelPos2D.ll) && array.Contains(RelPos2D.gg)) ||
                (array.Contains(RelPos2D.lg) && array.Contains(RelPos2D.gl));
        }

        public Tuple<Point3D, Point3D, Point3D> Solve(Point3D Target, params Point3D[] Points)
        {
            Dictionary<Point3D, double> distanceToTarget = new Dictionary<Point3D, double>();
            Dictionary<Point3D, RelPos2D> relativePosition = new Dictionary<Point3D, RelPos2D>();
            List<int> visited = new List<int>();

            Dictionary<RelPos2D, int> countPerPosition = new Dictionary<RelPos2D, int>()
            {
               {RelPos2D.ee,0},
               {RelPos2D.eg,0},
               {RelPos2D.el,0},
               {RelPos2D.ge,0},
               {RelPos2D.gg,0},
               {RelPos2D.gl,0},
               {RelPos2D.le,0},
               {RelPos2D.lg,0},
               {RelPos2D.ll,0}
            };

            foreach (var point in Points)
            {
                distanceToTarget.Add(point, TrianglePointTools.Distance(point, Target));
                RelPos2D position = TrianglePointTools.RelativePosition(point, Target);
                relativePosition.Add(point, position);
                countPerPosition[position]++;
            }

            //check countPerPosition to see if there are solutions
            int pointsCount = Points.Length;
            bool noSolutions = false;
            foreach (var key in countPerPosition.Keys)
            {
                if (countPerPosition[key] == pointsCount)
                {
                    noSolutions = true;
                    break;
                }
            }

            noSolutions = noSolutions ||
                    countPerPosition[RelPos2D.ll] + countPerPosition[RelPos2D.le] + countPerPosition[RelPos2D.lg] == pointsCount ||
                    countPerPosition[RelPos2D.lg] + countPerPosition[RelPos2D.eg] + countPerPosition[RelPos2D.gg] == pointsCount ||
                    countPerPosition[RelPos2D.gg] + countPerPosition[RelPos2D.ge] + countPerPosition[RelPos2D.gl] == pointsCount ||
                    countPerPosition[RelPos2D.ll] + countPerPosition[RelPos2D.el] + countPerPosition[RelPos2D.gl] == pointsCount;

            if (noSolutions)
                throw new Exception("No solutions.");

            var orderedPoints = Points.OrderBy(point => distanceToTarget[point]);
            bool found = false;

            Point3D
                Point1 = null,
                Point2 = null,
                Point3 = null;

            RelPos2D PosPoint1,
                     PosPoint2,
                     PosPoint3;

            foreach (var point1 in orderedPoints)
            {
                Point1 = point1;
                PosPoint1 = relativePosition[Point1];

                var point2Candidates = orderedPoints.Where(p => p != Point1)
                                                    .OrderBy(p => distanceToTarget[p]);

                //this should not happen because we know that we have at least one solution
                if (point2Candidates.Count() == 0)
                    continue;

                foreach (var point2 in point2Candidates)
                {
                    Point2 = point2;
                    PosPoint2 = relativePosition[Point2];

                    var point3ValidPositions = ValidPositions(PosPoint1, PosPoint2);

                    var point3Candidates = orderedPoints.Where(p => p != Point1 && p != Point2 && point3ValidPositions.Contains(relativePosition[p]))
                                                        .OrderBy(p => distanceToTarget[p]);

                    if (point3Candidates.Count() == 0)
                        continue;

                    foreach (var point3 in point3Candidates)
                    {
                        Point3 = point3;
                        PosPoint3 = relativePosition[Point3];

                        //check if already visited
                        //hash subject to conflicts
                        var hash = Point1.GetHashCode() *
                                   Point2.GetHashCode() *
                                   Point3.GetHashCode();

                        if (visited.Contains(hash))
                            continue;

                        if (UncertainSolution(PosPoint1, PosPoint2, PosPoint3))
                        {
                            found = TrianglePointTools.TriangleContainsPoint(Point1, Point2, Point3, Target);
                        }
                        else
                        {
                            found = true;
                        }

                        if (found)
                            break;

                        visited.Add(hash);
                    }

                    if (found)
                        break;
                }

                if (found)
                    break;
            }

            if (found)
                return new Tuple<Point3D, Point3D, Point3D>(Point1, Point2, Point3);

            throw new Exception("No solutions.");
        }
    }
}
