// Consider dy/dt = f(y) = f1(y) + f2(y) + ...

// a class to represent y(x,t)
class State
{
  // some stuff here that stores values of y, spatial location information, etc.
  
  // some stuff here that overloads addition and scalar multiplication
};

// an interface to represent g(y) that f or f1 might dependent on
class DiagnosticQuantity
{
public:
  virtual State eval(const State& y) const = 0;
};

// an interface to represent either f(y) or fi(y), i = 1, 2, ...
class Tendency
{
public:
  virtual State eval(const State& y) const = 0;
};

/*************************************************************************************/

// consider solving dy/dt = -y + cos(y)

// a Tendency derived class that represents f1(y) = -y
class Term1 : public Tendency
{
public:
  State eval(const State& y) const override
  {
    return -y;
  }
};

// a Tendency derived class that represents f2(y) = cos(y)
class Term2 : public Tendency
{
public:
  State eval(const State& y) const override
  {
    // code to compute cos(y) from State y and return as a State
  }
};

// a Tendency dervived class that represents f2(y) = cos(y) but only computes cos(y) once
class StaleTerm2 : public Tendency
{
  State result;
  
public:

  StaleTerm2(const State& y)
  {
    // code to compute cos(y) from State y and store as result
   }

  State eval(const State& y) const override
  {
    return result;
  }
};

// a Tendency drived class that represents f(y) = f1(y) + f2(y)
class TwoProcessTendency : public Tendency
{
  const Tendency& f1;
  const Tendency& f2;
  
public:

  RHS(const Tendency& f1, const Tendency& f2) : f1(f1), f2(f2) {}
  
  State eval(const State& y) const override
  {
    return f1.eval(y) + f2.eval(y);
  }
};

void setInitialState(State& y)
{
  // some code here
}

double computeError(State& y)
{
  // some code here
}

int main(int argc, char *argv[])
{
  State y;
  
  Term1 f1;
  Term2 f2;
  
  const double dt = 0.1;
  const double numSteps = 10;
  
  // compute using forward Euler
  setInitialState(y);
  TwoProcessTendency f(f1, f2);
  for (int n=0; n < numSteps; n++)
  {
    y = y + dt * f.eval(y);
  }
  std::cout << computeError(y) << std::endl;
  
  // compute using forward Euler with stale cos(y)
  setInitialState(y);
  StaleTerm2 staleF2(y);
  TwoProcessTendency staleF(f1, staleF2);
  for (int n=0; n < numSteps; n++)
  {
    y = y + dt * staleF.eval(y);
  }
  std::cout << computeError(y) << std::endl;
  
  return 0;
}
  
    

  
  
  







