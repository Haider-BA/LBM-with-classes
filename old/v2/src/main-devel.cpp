#include "UnitTest++.h"
#include "UnitTestCustomUtilities.hpp"

int main()
{
  return Unfit::RunAllTheTests();
//  return Unfit::RunOneSuite("FunctionalityAndExceptionTests");
//  return Unfit::RunOneSuite("PrintDemo");
//  return Unfit::RunOneSuite("AnalyticalSolutionTests");
//  return Unfit::RunOneTest("DiffusionEquation");
//  return Unfit::RunOneTest("CoupledNSCDE");
//  return Unfit::RunOneTest("NavierStokesEquation");
//  return Unfit::RunOneTest("PoiseuilleFlow");
//  return Unfit::RunOneTest("DiffusionEquationWithObstacles");
//  return Unfit::RunOneTest("ObstacleNormalReflect");

  // to run only one suite, uncomment the following line
  // and specify the name of the suite. Also, comment out the RunAllTheTests()
  // return Unfit::RunOneSuite("suite_name");

  // to run only one test, uncomment the following line
  // and specify the name of the test.Also, comment out the RunAllTheTests()
  // return Unfit::RunOneTest("test_name");
}
