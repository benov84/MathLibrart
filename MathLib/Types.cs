using System;
using System.Collections.Generic;
using System.Text;

namespace Benov.MathLib
{
    public class Line
    {
        public Point StartPoint { get; set; }
        public Point EndPoint { get; set; }

        public Line(Point startPoint, Point endPoint)
        {
            StartPoint = startPoint;
            EndPoint = endPoint;
        }
    }

    public enum Side { Left, Right, Straight }
}
