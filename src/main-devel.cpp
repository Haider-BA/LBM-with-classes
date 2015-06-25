#include "UnitTest++.h"
#include "UnitTestCustomUtilities.hpp"

int main()
{
//  return Unfit::RunAllTheTests();
  return Unfit::RunOneSuite("TestFunctionality");
//  return Unfit::RunOneSuite("TestAnalyticalSolutions");
//  return Unfit::RunOneTest("SimulateConvectionDiffusion");
//  return Unfit::RunOneTest("SimulatePoiseuilleFlow");
//  return Unfit::RunOneTest("SimulateNSCDCoupling");
//  return Unfit::RunOneTest("SimulateNSCDCouplingWithObstacles");
//  return Unfit::RunOneTest("AnalyticalPoiseuille");

  // to run only one suite, uncomment the following line
  // and specify the name of the suite. Also, comment out the RunAllTheTests()
  // return Unfit::RunOneSuite("suite_name");

  // to run only one test, uncomment the following line
  // and specify the name of the test.Also, comment out the RunAllTheTests()
  // return Unfit::RunOneTest("test_name");
}
