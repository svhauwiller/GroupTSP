using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;

namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your node data structure and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }


            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        private const int CITY_ICON_SIZE = 5;

        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;
        #endregion

        #region Public members
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed); 
            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        //public void GenerateProblem(int size) // unused
        //{
        //   this.GenerateProblem(size, Modes.Normal);
        //}

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        private class GreedyRoute
        {
            public GreedyRoute(City[] list)
            {
                Cities = list;
                visited = new bool[Cities.Length];
                for (int i = 0; i < visited.Length; i++)
                    visited[i] = false;
                r = new Random();
                startIndex = r.Next() % Cities.Length;
                numVisited = 1;
                route = new ArrayList();
            }

            public void SetRoute()
            {
                City currentCity = Cities[startIndex];
                route.Add(Cities[startIndex]);
                visited[startIndex] = true;
                City nextDest = null;
                while (numVisited < Cities.Length)
                {
                    double minCost = Double.PositiveInfinity;
                    int currentIndex = -1;
                    for (int i = 0; i < Cities.Length; i++)
                    {
                        if (!visited[i])
                        {
                            double tempCost = currentCity.costToGetTo(Cities[i]);
                            if (tempCost < minCost)
                            {
                                minCost = tempCost;
                                nextDest = Cities[i];
                                currentIndex = i;
                            }
                        }
                    }
                    if (minCost == Double.PositiveInfinity)
                        break;
                    else
                    {
                        route.Add(nextDest);
                        visited[currentIndex] = true;
                        currentCity = nextDest;
                    }
                }
                if (currentCity.costToGetTo(Cities[startIndex]) == Double.PositiveInfinity)
                    route = null;
                else if (route.Count < Cities.Length)
                    route = null;
            }

            public ArrayList GetRoute()
            {
                return route;
            }

            public void Reset()
            {
                startIndex = r.Next();
            }

            private Random r;
            private City[] Cities;
            private bool[] visited;
            private int startIndex;
            private int numVisited;
            private ArrayList route;
        }

        /// <summary>
        ///  solve the problem.  This is the entry point for the solver when the run button is clicked
        /// right now it just picks a simple solution. 
        /// </summary>
        public void solveProblem(string algorithm)
        {
            int x;
            Route = new ArrayList(); 
            // this is the trivial solution. 
            if (algorithm == "default")
            {
                for (x = 0; x < Cities.Length; x++)
                {
                    Route.Add(Cities[Cities.Length - x - 1]);
                }
            }
            else if (algorithm == "greedy")
            {
                Route = null;
                GreedyRoute greedy = new GreedyRoute(Cities);
                while (Route == null)
                {
                    greedy.SetRoute();
                    Route = greedy.GetRoute();
                    if (Route == null)
                        greedy.Reset();
                }
            }
            else if (algorithm == "random")
            {
                for (int i = 0; i < Cities.Length; i++)
                {
                    Route.Add(Cities[i]);
                }
            }
            else if (algorithm == "custom")
            {
                Route.Clear();

                //Check for an optimal base case if the problem has 3 or less cities
                if (Cities.Length <= 3)
                {
                    if (Cities.Length > 0)
                    {
                        Route.Add(Cities[0]);

                        if (Cities.Length == 2)
                        {
                            Route.Add(Cities[1]);
                        }
                        else if (Cities.Length == 3)
                        {
                            double forwardDist = 0;
                            double backwardDist = 0;

                            forwardDist += Cities[0].costToGetTo(Cities[1]);
                            forwardDist += Cities[1].costToGetTo(Cities[2]);
                            forwardDist += Cities[2].costToGetTo(Cities[0]);

                            backwardDist += Cities[0].costToGetTo(Cities[2]);
                            backwardDist += Cities[2].costToGetTo(Cities[1]);
                            backwardDist += Cities[1].costToGetTo(Cities[0]);

                            if (forwardDist < backwardDist)
                            {
                                Route.Add(Cities[1]);
                                Route.Add(Cities[2]);
                            }
                            else
                            {
                                Route.Add(Cities[2]);
                                Route.Add(Cities[1]);
                            }
                        }
                    }

                    bssf = new TSPSolution(Route);
                    Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
                    Program.MainForm.Invalidate();
                    return;
                }

                //Add the first city
                Route.Add(Cities[0]);
                List<double> distFromFirstCity = new List<Double>();
                distFromFirstCity.Add(0);

                //Add the second city as the furthest city from the first
                double maxDist = 0;
                City maxDistCity = null;
                int secondCityIndex = 0;
                for (int i = 1; i < Cities.Length; i++)
                {
                    distFromFirstCity.Add(Cities[0].costToGetTo(Cities[i]));
                    if (!Double.IsPositiveInfinity(distFromFirstCity[i]))
                    {
                        if (distFromFirstCity[i] > maxDist)
                        {
                            maxDist = distFromFirstCity[i];
                            maxDistCity = Cities[i];
                            secondCityIndex = i;
                        }
                    }
                }
                Route.Add(maxDistCity);

                //Add the Third City as the furthest city from both of them that will complete the cycle
                maxDist = 0;
                maxDistCity = null;
                int thirdCityIndex = 1;
                for (int i = 1; i < Cities.Length; i++)
                {
                    if (i == secondCityIndex)
                    {
                        continue;
                    }
                    double distFromSecondCity = Cities[secondCityIndex].costToGetTo(Cities[i]);
                    double distToOrigin = Cities[i].costToGetTo(Cities[0]);
                    if (!Double.IsPositiveInfinity(distFromSecondCity) &&
                        !Double.IsPositiveInfinity(distToOrigin))
                    {
                        if (Math.Min(distFromFirstCity[i], distFromSecondCity) > maxDist)
                        {
                            maxDist = Math.Min(distFromFirstCity[i], distFromSecondCity);
                            maxDistCity = Cities[i];
                            thirdCityIndex = i;
                        }
                    }
                }
                Route.Add(maxDistCity);
                double dist = new TSPSolution(Route).costOfRoute();
                List<int> cities = new List<int>();
                cities.Add(0);
                cities.Add(secondCityIndex);
                cities.Add(thirdCityIndex);

                while (Route.Count < Cities.Length) {
                    addNextCity(cities);
                }
            }

            // call this the best solution so far.  bssf is the route that will be drawn by the Draw method. 
            bssf = new TSPSolution(Route);
            // update the cost of the tour. 
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            // do a refresh. 
            Program.MainForm.Invalidate();

        }

        private void addNextCity(List<int> citiesAdded) {
            int bestCity = -1;
            double min = Double.PositiveInfinity;
            int city = -1;
            
            for (int k = 0; k < citiesAdded.Count; k++) {
                for (int i = 1; i < Cities.Length; i++) {
                    int next = 0;
                    if (k + 1 == citiesAdded.Count) {
                        next = citiesAdded[0];
                    } else {
                        next = citiesAdded[k + 1];
                    }
                    if (citiesAdded.Contains(i)) {
                        continue;
                    }
                    double pathLength = Cities[citiesAdded[k]].costToGetTo(Cities[next]);//Length of the path we are considering removing
                    double dist1 = Cities[citiesAdded[k]].costToGetTo(Cities[i]);//Length of path from start city to the city we are considering adding to the path
                    if (next == Cities.Length) {
                        next = 0;
                    }
                    double dist2 = Cities[i].costToGetTo(Cities[next]);//Length of path from city we are considering adding to the next city in the path
                    double increase = dist1 + dist2 - pathLength;
                    if (increase < min) {
                        bestCity = k;
                        min = increase;
                        city = i;
                    }
                }
                k++;//Step through the cities to repeat looking for a valid point to add if there wasn't one with the initially selected point
            }
            if (bestCity == -1){
                throw new Exception("Could not find a possible city to add to the current path");
            } else {
                citiesAdded.Insert(bestCity + 1, city);
                Route.Insert(bestCity + 1, Cities[city]);
            }
        }

        #endregion
    }

}
