struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 0.000001;
  }
};

struct FunctionToApproximate  {
  double operator() (double x)  {
    return x*x - 3*x + 1;  // Replace with your function
  }
};

// ...
using boost::math::tools::bisect;
double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
double to = 1;
std::pair<double, double> result = bisect(FunctionToApproximate(), from, to, TerminationCondition());
double root = (result.first + result.second) / 2;  // = 0.381966...


class root {
public:
  double f(double x) {
    return (x * x - 1);
  }
  void findRoot() {
    double a = 0.0;
    double b = 10.0;
    tolerance tol = 0.00001;

    typedef double (root::*MemFn)(double x);
    std::pair<double, double> found = boost::math::tools::bisect(&root::f,
        a, b, tol);
    std::cout << "==> x = [" << found.first << ',' << found.second << "]\n";
  }
};

