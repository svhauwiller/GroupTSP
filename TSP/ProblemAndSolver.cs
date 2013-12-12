using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;

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

        //This is the class that implements the greedy solution for TSP
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
                if (startIndex >= Cities.Length)
                    startIndex = startIndex % Cities.Length;
                
                //start with a random city
                City currentCity = Cities[startIndex];
                route.Add(Cities[startIndex]);
                visited[startIndex] = true;
                City nextDest = null;
                while (numVisited < Cities.Length)
                {
                    double minCost = Double.PositiveInfinity;
                    int currentIndex = -1;

                    //look for the city that has the lowest cost to get to
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

                    //break if a dead-end city is reached
                    if (minCost == Double.PositiveInfinity)
                        break;
                    else
                    {
                        route.Add(nextDest);
                        visited[currentIndex] = true;
                        currentCity = nextDest;
                    }
                }

                //set the route to null if not completed
                if (currentCity.costToGetTo(Cities[startIndex]) == Double.PositiveInfinity)
                    route = null;
                else if (route.Count < Cities.Length)
                    route = null;
            }

            public ArrayList GetRoute()
            {
                return route;
            }

            //get ready to form a new route
            public void Reset()
            {
                startIndex = r.Next();
                route = new ArrayList();
            }

            private Random r;
            private City[] Cities;
            private bool[] visited;
            private int startIndex;
            private int numVisited;
            private ArrayList route;
        }

        //this class implements the branch and bound algorithm using include/exclude states
        private class BranchBound
        {
            public BranchBound(City[] list)
            {
                Cities = list;
                double lBound = InitializeStartMatrix(list);

                bssf = InitialBSSF(list);

                //the state object accepts a matrix object, the number of cities,
                //and its level in the branch and bound tree
                currState = new State(currMatrix, list.Length, 0);
                currState.SetLowerBound(lBound);
                data = new SolutionData();

                //the heap is the priority queue
                heap = new BinaryHeap();
                heap.Insert(currState);

                //the data variable is an object that stores the data needed to fill the table
                data.SetNumMaxStored(1);
                data.AddNumCreatedStates(1);
            }

            //this is the method that finds the solution and returns an object which contains the data for the table
            public SolutionData Solve()
            {
                Stopwatch watch = new Stopwatch();
                watch.Start();

                //the timer is started just before the loop to find the solution begins
                while (!heap.IsEmpty())
                {
                    //the state with the lowest bound is removed from the PQ
                    currState = heap.Delete(0);
                    Matrix currMatrix = currState.GetMatrix();
                    State nextInclude = null;
                    State nextExclude = null;
                    double maxDiff = 0D;

                    int indexI = -1;
                    int indexJ = -1;

                    bool first = true;

                    //this loop cycles through the matrix to compare
                    //the bound differences of edges that have '0' values
                    for (int i = 0; i < Cities.Length; i++)
                    {
                        for (int j = 0; j < Cities.Length; j++)
                        {
                            if (currMatrix.GetValue(i, j) == 0)
                            {

                                //the edge in the matrix is initialized with the first '0'
                                //edge found just in case there turns out to be no edges
                                //which have a bound(exclude - include) > 0
                                if (first)
                                {
                                    indexI = i;
                                    indexJ = j;
                                    first = false;
                                }
                                double tempDiff = GetBoundExclude(i, j) - GetBoundInclude(i, j);
                                if (tempDiff > maxDiff)
                                {
                                    maxDiff = tempDiff;
                                    indexI = i;
                                    indexJ = j;
                                }
                            }
                        }
                    }

                    nextInclude = GetIncludeState(indexI, indexJ);
                    nextExclude = GetExcludeState(indexI, indexJ);

                    //the nextIncludeState should be one level lower 
                    //in the tree than the current state
                    nextInclude.IncrementLevel();

                    bool bottomLevel = false;

                    //check if a full solution has been reached
                    if (nextInclude.GetLevel() == Cities.Length)
                    {
                        bottomLevel = true;

                        //if the solution is less than the current bssf,
                        //change the bssf, prune any edges that can be pruned,
                        //and set the current solution route
                        if (nextInclude.GetLowerBound() < bssf)
                        {
                            bssf = nextInclude.GetLowerBound();
                            route = new ArrayList();
                            int[] exited = nextInclude.GetExited();
                            int index = 0;
                            for (int i = 0; i < Cities.Length; i++)
                            {
                                route.Add(Cities[index]);
                                index = exited[index];
                            }
                            int pruned = heap.CleanQueue(bssf);
                            data.AddNumStatesPruned(pruned);
                            data.IncrementNumBSSFUpdates();
                        }
                    }

                    data.AddNumCreatedStates(2);

                    //if a solution has been reached, don't add the include
                    //and exclude states to the PQ
                    if (bottomLevel)
                        continue;

                    if (nextInclude.GetLowerBound() <= bssf)
                        heap.Insert(nextInclude);
                    else
                        data.AddNumStatesPruned(1);
                    if (nextExclude.GetLowerBound() <= bssf)
                        heap.Insert(nextExclude);
                    else
                        data.AddNumStatesPruned(1);

                    if (heap.GetHeapCount() > data.GetNumMaxStored())
                        data.SetNumMaxStored(heap.GetHeapCount());      
                }

                watch.Stop();
                data.setTime(watch.ElapsedMilliseconds / 1000.0);

                return data;
            }

            //based on the easy way to examine the 0's as given in the specs
            private double GetBoundInclude(int i, int j)
            {
                Matrix currMatrix = currState.GetMatrix();
                double addReduct = 0D;

                for (int n = 0; n < Cities.Length; n++)
                {
                    if (n != j && currMatrix.GetValue(i, n) == 0)
                    {
                        double minCol = Double.PositiveInfinity;
                        for (int m = 0; m < Cities.Length; m++)
                        {
                            if (m != i)
                            {
                                if (currMatrix.GetValue(m, n) < minCol)
                                    minCol = currMatrix.GetValue(m,n);
                            }
                        }
                        addReduct += minCol;
                    }
                    if (n != i && currMatrix.GetValue(n, j) == 0)
                    {
                        double minRow = Double.PositiveInfinity;
                        for (int m = 0; m < Cities.Length; m++)
                        {
                            if (m != j)
                            {
                                if (currMatrix.GetValue(n, m) < minRow)
                                    minRow = currMatrix.GetValue(n,m);
                            }
                        }
                        addReduct += minRow;
                    }
                }

                return currState.GetLowerBound() + addReduct; 
            }

            //based on the easy way to examine the 0's as given in the specs
            private double GetBoundExclude(int i, int j)
            {
                Matrix currMatrix = currState.GetMatrix();

                double minRow = Double.PositiveInfinity;
                double minColumn = Double.PositiveInfinity;

                for (int n = 0; n < Cities.Length; n++)
                {
                    if (n != j)
                    {
                        if (currMatrix.GetValue(i, n) < minRow)
                            minRow = currMatrix.GetValue(i, n);
                    }
                    if (n != i)
                    {
                        if (currMatrix.GetValue(n, j) < minColumn)
                            minColumn = currMatrix.GetValue(n, j);
                    }
                }

                return currState.GetLowerBound() + minRow + minColumn;
            }

            //change the matrix to account for adding the edge (i,j)
            private State GetIncludeState(int i, int j)
            {
                Matrix inc = new Matrix(currState.GetMatrix().GetCities());
                inc.SetRowsColumns(currState.GetMatrix().GetRows(), currState.GetMatrix().GetColumns());
                inc.CopyMatrix(currState.GetMatrix().GetMatrix());
                double includeBound = currState.GetLowerBound();

                includeBound += inc.GetValue(i, j);
                inc.InfiniteRow(i);
                inc.InfiniteColumn(j);
                inc.IgnoreRowColumn(i, j);
                
                State include = new State(inc, Cities.Length, currState.GetLevel());
                include.SetLowerBound(currState.GetLowerBound());
                include.SetEnteredExited(currState.GetEntered(), currState.GetExited());
                if (include.GetLevel() < Cities.Length - 2)
                    include.DeleteEdges(i, j);
                else
                {
                    include.setExited(i, j);
                }
                inc = include.GetMatrix();
                includeBound += inc.ReduceMatrix();

                include.SetMatrix(inc);
                include.SetLowerBound(includeBound);

                return include;
            }

            //change the matrix to account for excluding the edge (i,j)
            private State GetExcludeState(int i, int j)
            {
                Matrix exc = new Matrix(currState.GetMatrix().GetCities());
                exc.SetRowsColumns(currState.GetMatrix().GetRows(), currState.GetMatrix().GetColumns());
                exc.CopyMatrix(currState.GetMatrix().GetMatrix());
                double excludeBound = currState.GetLowerBound();

                exc.ChangeValue(i, j, Double.PositiveInfinity);
                excludeBound += exc.ReduceMatrix();

                State exclude = new State(exc, Cities.Length, currState.GetLevel());
                exclude.SetLowerBound(currState.GetLowerBound());
                exclude.SetEnteredExited(currState.GetEntered(), currState.GetExited());
                exclude.SetMatrix(exc);
                exclude.SetLowerBound(excludeBound);

                return exclude;
            }

            //set the start matrix, doing the initial reduction
            private double InitializeStartMatrix(City[] list)
            {
                List<int> rows = new List<int>();
                List<int> columns = new List<int>();
                for (int i = 0; i < list.Length; i++)
                {
                    rows.Add(i);
                    columns.Add(i);
                }

                currMatrix = new Matrix(list);
                return currMatrix.ReduceMatrix();
            }

            //set the initial bssf, comparing several iterations 
            //of the greedy algorithm for the best choice
            private double InitialBSSF(City[] list)
            {
                double initial = Double.PositiveInfinity;
                for (int i = 0; i < 5; i++)
                {
                    ArrayList Route = null;
                    GreedyRoute greedy = new GreedyRoute(Cities);
                    while (Route == null)
                    {
                        greedy.SetRoute();
                        Route = greedy.GetRoute();
                        if (Route == null)
                            greedy.Reset();
                    }

                    TSPSolution _bssf;
                    _bssf = new TSPSolution(Route);
                    double temp = _bssf.costOfRoute();

                    if (temp < initial)
                    {
                        initial = temp;
                        route = Route;
                    }
                }

                return initial;
            }

            public ArrayList GetRoute()
            {
                return route;
            }

            private City[] Cities;
            ArrayList route;
            private double bssf;
            private Matrix currMatrix;
            private State currState;
            private SolutionData data;
            private BinaryHeap heap;
        }

        //priority queue implemented as a binary heap
        private class BinaryHeap
        {
            private List<State> heapArray;
            private int heapCount;

            public BinaryHeap()
            {
                heapCount = 0;

                //use a dynamic list, cannot know how many
                //states will be stored
                heapArray = new List<State>();
            }

            public int Insert(State state)
            {
                heapArray.Add(state);

                int i = heapCount;

                while (i != 0 && heapArray[i].GetLowerBound() <= heapArray[i / 2].GetLowerBound())
                {

                    //if two states have an equivalent lower bound,
                    //check what level they are in the BB tree, if
                    //one of them is closer to a solution, make
                    //it a priority over the other
                    if (heapArray[i].GetLowerBound() == heapArray[i / 2].GetLowerBound())
                    {
                        if (heapArray[i].GetLevel() > heapArray[i / 2].GetLevel())
                        {
                            heapArray[i] = heapArray[i / 2];
                            heapArray[i / 2] = state;
                            i = i / 2;
                        }
                        else
                            break;
                    }
                    else
                    {
                        heapArray[i] = heapArray[i / 2];
                        heapArray[i / 2] = state;
                        i = i / 2;
                    }
                }

                heapCount++;
                return i;
            }                

            public State Delete(int index)
            {
                State del = heapArray[index];

                heapArray[index] = heapArray[heapCount - 1];
                heapArray.RemoveAt(heapCount - 1);

                heapCount--;

                if (heapCount <= 1)
                    return del;

                int i = index;

                if (i * 2 >= heapCount)
                {
                    return del;
                }
                if (i * 2 + 1 >= heapCount)
                {
                    if (heapArray[i].GetLowerBound() >= heapArray[i * 2].GetLowerBound())
                    {
                        if (heapArray[i].GetLowerBound() == heapArray[i * 2].GetLowerBound())
                        {
                            if (heapArray[i].GetLevel() < heapArray[i * 2].GetLevel())
                            {
                                State temp = heapArray[i * 2];
                                heapArray[i * 2] = heapArray[i];
                                heapArray[i] = temp;
                            }
                        }
                        else
                        {
                            State temp = heapArray[i * 2];
                            heapArray[i * 2] = heapArray[i];
                            heapArray[i] = temp;
                        }
                    }
                    return del;
                }
                bool countBreak = false;

                while (heapArray[i].GetLowerBound() >= heapArray[i * 2].GetLowerBound() ||
                       heapArray[i].GetLowerBound() >= heapArray[i * 2 + 1].GetLowerBound())
                {
                    if (heapArray[i * 2].GetLowerBound() < heapArray[i * 2 + 1].GetLowerBound())
                    {
                        if (heapArray[i].GetLowerBound() == heapArray[i * 2].GetLowerBound())
                        {
                            if (heapArray[i].GetLevel() < heapArray[i * 2].GetLevel())
                            {
                                State temp = heapArray[i * 2];
                                heapArray[i * 2] = heapArray[i];
                                heapArray[i] = temp;
                            }
                            else
                                break;
                        }
                        else
                        {
                            State temp = heapArray[i * 2];
                            heapArray[i * 2] = heapArray[i];
                            heapArray[i] = temp;
                        }

                        i = i * 2;
                    }
                    else
                    {
                        if (heapArray[i].GetLowerBound() == heapArray[i * 2 + 1].GetLowerBound())
                        {
                            if (heapArray[i].GetLevel() < heapArray[i * 2 + 1].GetLevel())
                            {
                                State temp = heapArray[i * 2 + 1];
                                heapArray[i * 2 + 1] = heapArray[i];
                                heapArray[i] = temp;
                            }
                            else
                                break;
                        }
                        else
                        {
                            State temp = heapArray[i * 2 + 1];
                            heapArray[i * 2 + 1] = heapArray[i];
                            heapArray[i] = temp;
                        }

                        i = i * 2 + 1;
                    }

                    if (i * 2 + 1 >= heapCount)
                    {
                        countBreak = true;
                        break;
                    }
                }

                if (i * 2 < heapCount && countBreak)
                {
                    if (heapArray[i].GetLowerBound() < heapArray[i * 2].GetLowerBound())
                        return del;
                    else if (heapArray[i].GetLowerBound() == heapArray[i * 2].GetLowerBound())
                    {
                        if (heapArray[i].GetLevel() < heapArray[i * 2].GetLevel())
                        {
                            State temp = heapArray[i * 2];
                            heapArray[i * 2] = heapArray[i];
                            heapArray[i] = temp;
                        }
                    }
                    else
                    {
                        State temp = heapArray[i * 2];
                        heapArray[i * 2] = heapArray[i];
                        heapArray[i] = temp;
                    }
                } 

                return del;
            }

            //prune any states with a higher bound than the bssf
            public int CleanQueue(double best)
            {
                int numberPruned = 0;
                int initialCount = heapArray.Count;

                for (int i = 0; i < initialCount; i++)
                {
                    if (heapArray[i - numberPruned].GetLowerBound() > best)
                    {
                        Delete(i - numberPruned);
                        numberPruned++;
                    }
                }

                return numberPruned;
            }

            public bool IsEmpty()
            {
                if (heapCount > 0)
                    return false;
                else
                    return true;
            }

            public int GetHeapCount()
            {
                return heapCount;
            }
        }

        //class for storing data needed for table
        private class SolutionData
        {
            public SolutionData()
            {
                numMaxStoredStates = 0;
                numBSSFUpdates = 0;
                numStatesCreated = 0;
                numStatesPruned = 0;
                timeToSolve = 0D;
                optimal = true;
            }

            private int numMaxStoredStates;
            private int numBSSFUpdates;
            private int numStatesCreated;
            private int numStatesPruned;
            private double timeToSolve;
            private bool optimal;

            public void SetNumMaxStored(int stored)
            {
                numMaxStoredStates = stored;
            }

            public int GetNumMaxStored()
            {
                return numMaxStoredStates;
            }

            public void IncrementNumBSSFUpdates()
            {
                numBSSFUpdates = numBSSFUpdates + 1;
            }

            public int GetNumBSSFUpdates()
            {
                return numBSSFUpdates;
            }

            public void AddNumCreatedStates(int num)
            {
                numStatesCreated += num;
            }

            public int GetNumStatesCreated()
            {
                return numStatesCreated;
            }

            public void AddNumStatesPruned(int num)
            {
                numStatesPruned += num;
            }

            public int GetNumStatesPruned()
            {
                return numStatesPruned;
            }

            public void SetNotOptimal()
            {
                optimal = false;
            }

            public bool GetOptimal()
            {
                return optimal;
            }

            public void setTime(double time)
            {
                timeToSolve = time;
            }

            public double GetTime()
            {
                return timeToSolve;
            }
        }

        //represents a state of the branch and bound algorithm
        private class State
        {
            public State(Matrix m, int n, int lvl)
            {
                matrix = m;
                lowerBound = 0D;
                level = lvl;
                numCities = n;
                entered = new int[numCities];
                exited = new int[numCities];
                for (int i = 0; i < numCities; i++)
                {
                    entered[i] = -1;
                    exited[i] = -1;
                }
            }

            public void SetMatrix(Matrix m)
            {
                matrix = m;
            }

            public Matrix GetMatrix()
            {
                return matrix;
            }

            public void SetLowerBound(double low)
            {
                lowerBound = low;
            }

            public double GetLowerBound()
            {
                return lowerBound;
            }

            public void IncrementLevel()
            {
                level += 1;
            }

            public int GetLevel()
            {
                return level;
            }

            public void setExited(int i, int j)
            {
                exited[i] = j;
            }

            public void SetEnteredExited(int[] en, int[] ex)
            {
                for (int i = 0; i < numCities; i++)
                {
                    entered[i] = en[i];
                    exited[i] = ex[i];
                }
            }

            public int[] GetEntered()
            {
                return entered;
            }

            public int[] GetExited()
            {
                return exited;
            }

            //handle premature cycles
            public void DeleteEdges(int i, int j)
            {
                entered[j] = i;
                exited[i] = j;
                int start = i;
                int end = j;

                while (exited[end] != -1)
                    end = exited[end];
                while (entered[start] != -1)
                    start = entered[start];

                if (level < entered.Length - 1)
                {
                    while (start != j)
                    {
                        matrix.ChangeValue(end, start, Double.PositiveInfinity);
                        matrix.ChangeValue(j, start, Double.PositiveInfinity);
                        start = exited[start];
                    }
                }
            }

            private Matrix matrix;
            private double lowerBound;
            private int numCities;
            private int level;
            private int[] entered;
            private int[] exited;
        }

        //stores a matrix of edge cost values
        private class Matrix
        {
            public Matrix(City[] list)
            {
                Cities = list;
                rows = new bool[list.Length];
                columns = new bool[list.Length];
                for (int i = 0; i < list.Length; i++)
                {
                    rows[i] = true;
                    columns[i] = true;
                }
                InitializeMatrix();
            }

            public void CopyMatrix(double[,] m)
            {
                for (int i = 0; i < Cities.Length; i++)
                {
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        matrix[i, j] = m[i, j];
                    }
                }
            }

            private void InitializeMatrix()
            {
                matrix = new double[Cities.Length, Cities.Length];
                for (int i = 0; i < Cities.Length; i++)
                {
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        //a city cannot have an edge to itself
                        if (i == j)
                            matrix[i, j] = Double.PositiveInfinity;
                        else
                            matrix[i, j] = Cities[i].costToGetTo(Cities[j]); 
                    }
                }
            }

            public void InfiniteRow(int row)
            {
                for (int j = 0; j < Cities.Length; j++)
                    matrix[row, j] = Double.PositiveInfinity;
            }

            public void InfiniteColumn(int column)
            {
                for (int i = 0; i < Cities.Length; i++)
                    matrix[i, column] = Double.PositiveInfinity;
            }

            public void IgnoreRowColumn(int i, int j)
            {
                rows[i] = false;
                columns[j] = false;
            }

            //ensure that there are 0's in each applicable
            //column and row
            public double ReduceMatrix()
            {
                double total = 0D;
                double min = Double.PositiveInfinity;
                bool zeroFound = false;
                for (int i = 0; i < rows.Length; i++)
                {
                    if (rows[i])
                    {
                        for (int j = 0; j < columns.Length; j++)
                        {

                            double temp = matrix[i, j];
                            if (temp < min)
                                min = temp;
                            if (temp == 0)
                                zeroFound = true;
                        }
                        if (!zeroFound)
                        {
                            total += min;
                            for (int j = 0; j < columns.Length; j++)
                                matrix[i, j] = matrix[i, j] - min;
                        }
                        else
                            zeroFound = false;

                        min = Double.PositiveInfinity; 
                    }
                }
                for (int j = 0; j < columns.Length; j++)
                {
                    if (columns[j])
                    {
                        for (int i = 0; i < rows.Length; i++)
                        {
                            double temp = matrix[i, j];
                            if (temp < min)
                                min = temp;
                            if (temp == 0)
                                zeroFound = true;
                        }
                        if (!zeroFound)
                        {
                            total += min;
                            for (int i = 0; i < rows.Length; i++)
                                matrix[i, j] = matrix[i, j] - min;
                        }
                        else
                            zeroFound = false;

                        min = Double.PositiveInfinity;
                    }
                }

                return total;
            }

            public void ChangeValue(int i, int j, double value)
            {
                matrix[i, j] = value;
            }

            public double GetValue(int i, int j)
            {
                return matrix[i, j];
            }

            public double[,] GetMatrix()
            {
                return matrix;
            }

            public City[] GetCities()
            {
                return Cities;
            }

            public bool[] GetRows()
            {
                return rows;
            }

            public bool[] GetColumns()
            {
                return columns;
            }

            public void SetRowsColumns(bool[] r, bool[] c)
            {
                for (int i = 0; i < Cities.Length; i++)
                {
                    rows[i] = r[i];
                    columns[i] = c[i];
                }
            }

            private City[] Cities;
            private double[,] matrix;
            private bool[] rows;
            private bool[] columns;
        }

        /// <summary>
        ///  solve the problem.  This is the entry point for the solver when the run button is clicked
        /// right now it just picks a simple solution. 
        /// </summary>
        public void solveProblem(string algorithm)
        {
            int x;
            Route = new ArrayList();
            SolutionData data = null;
            Stopwatch watch = new Stopwatch();
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
                data = new SolutionData(); 
                watch.Start();
                Route = null;
                GreedyRoute greedy = new GreedyRoute(Cities);
                while (Route == null)
                {
                    greedy.SetRoute();
                    Route = greedy.GetRoute();
                    if (Route == null)
                        greedy.Reset();
                }
                watch.Stop();
                data.setTime(watch.ElapsedMilliseconds / 1000.0);
            }
            else if (algorithm == "random")
            {
                for (int i = 0; i < Cities.Length; i++)
                {
                    Route.Add(Cities[i]);
                }
            }
            else if (algorithm == "branchBound")
            {
                watch.Start();
                Route = null;
                BranchBound bb = new BranchBound(Cities);
                data = bb.Solve();
                Route = bb.GetRoute();
                watch.Stop();
                data.setTime(watch.ElapsedMilliseconds / 1000.0);
            }
            else if (algorithm == "custom")
            {
                data = new SolutionData();
                watch.Start();
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

                while (Route.Count < Cities.Length)
                {
                    addNextCity(cities);
                }
                watch.Stop();
                data.setTime(watch.ElapsedMilliseconds / 1000.0);
            }

            // call this the best solution so far.  bssf is the route that will be drawn by the Draw method. 
            bssf = new TSPSolution(Route);
            // update the cost of the tour.
            if (algorithm != "random")
            {
                double timeSpent = Math.Round(data.GetTime(), 2);
                Program.MainForm.tbElapsedTime.Text = " " + timeSpent;
            }

            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            // do a refresh. 
            Program.MainForm.Invalidate();

        }

        private void addNextCity(List<int> citiesAdded)
        {
            int bestCity = -1;
            double min = Double.PositiveInfinity;
            int city = -1;

            for (int k = 0; k < citiesAdded.Count; k++)
            {
                for (int i = 1; i < Cities.Length; i++)
                {
                    int next = 0;
                    if (k + 1 == citiesAdded.Count)
                    {
                        next = citiesAdded[0];
                    }
                    else
                    {
                        next = citiesAdded[k + 1];
                    }
                    if (citiesAdded.Contains(i))
                    {
                        continue;
                    }
                    double pathLength = Cities[citiesAdded[k]].costToGetTo(Cities[next]);//Length of the path we are considering removing
                    double dist1 = Cities[citiesAdded[k]].costToGetTo(Cities[i]);//Length of path from start city to the city we are considering adding to the path
                    if (next == Cities.Length)
                    {
                        next = 0;
                    }
                    double dist2 = Cities[i].costToGetTo(Cities[next]);//Length of path from city we are considering adding to the next city in the path
                    double increase = dist1 + dist2 - pathLength;
                    if (increase < min)
                    {
                        bestCity = k;
                        min = increase;
                        city = i;
                    }
                }
                k++;//Step through the cities to repeat looking for a valid point to add if there wasn't one with the initially selected point
            }
            if (bestCity == -1)
            {
                throw new Exception("Could not find a possible city to add to the current path");
            }
            else
            {
                citiesAdded.Insert(bestCity + 1, city);
                Route.Insert(bestCity + 1, Cities[city]);
            }
        }

        #endregion
    }

}
