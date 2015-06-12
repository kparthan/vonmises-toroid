#include "Support.h"
#include "UniformRandomNumberGenerator.h"

extern UniformRandomNumberGenerator *uniform_generator;
extern Vector XAXIS,YAXIS,ZAXIS;
extern double MIN_N;

int main(int argc, char **argv)
{
  //Setup();
  XAXIS = Vector(3,0); XAXIS[0] = 1;
  YAXIS = Vector(3,0); YAXIS[1] = 1;
  ZAXIS = Vector(3,0); ZAXIS[2] = 1;

  MIN_N = 10;

  //cout << "ZERO: " << ZERO << endl;
  //cout << "HUGE_VAL: " << HUGE_VAL << endl;
  //cout << "1 - ZERO: " << 1-ZERO << endl;

  UniformReal uniform_distribution(0,1);
  RandomNumberGenerator generator;
  Generator num_gen(generator,uniform_distribution); 
  generator.seed(time(NULL)); // seed with the current time 
  uniform_generator = new UniformRandomNumberGenerator(num_gen);

  srand(time(NULL));

  struct Parameters parameters = parseCommandLineInput(argc,argv);

  if (parameters.test == SET) {
    TestFunctions();
  }


  delete(uniform_generator);

  return 0;
}

