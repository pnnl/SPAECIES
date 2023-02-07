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

class Cosine : public DiagnosticQuantity
{
public:
  State eval(const State& y) const override
  {
    // code to compute cos(y) from State y and return as a State
  }
};

class CosineStale : public DiagnosticQuantity
{
  State result;
  
public:

  CosineStale(const State& y)
  {
    // code to compute cos(y) from State y and store as result
   }

  State eval(const State& y) const override
  {
    return result;
  }
};

class Term2 : public Tendency
{
  const DiagnosticQuantity& cosine;

public:

  Term2(const DiagnosticQuantity& cosine) : cosine(cosine) {}
  
  State val(const State& y) const override
  {
    return cosine.eval(y);
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
  const double dt = 0.1;
  const double numSteps = 10;

  // compute using forward Euler
  {
    State y;
    setInitialState(y);
    const Term1 f1;
    const Cosine cosine;
    const Term2 f2(cosine);
    TwoProcessTendency f(f1, f2);
    for (int n=0; n < numSteps; n++)
    {
      y = y + dt * f.eval(y);
    }
    std::cout << computeError(y) << std::endl;
  }
  
  // compute using forward Euler with stale cos(y)
  {
    State y;
    setInitialState(y);
    const Term1 f1;
    const CosineStale cosine(y); // <---- note use of stale diagnostic quantity here
    const Term2 f2(cosine);
    TwoProcessTendency f(f1, f2);
    for (int n=0; n < numSteps; n++)
    {
      y = y + dt * f.eval(y);
    }
    std::cout << computeError(y) << std::endl;
  }
  
  return 0;
}
  
    

  
  
  







