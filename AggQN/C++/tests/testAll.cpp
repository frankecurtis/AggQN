// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cstdio>

#include "testOptions.hpp"
#include "testReporter.hpp"
#include "testVector.hpp"

// Main function
int main()
{

  // Initialize result
  int result = 0;

  // Run tests
  printf("testing Options.................... ");
  if (!testOptionsImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testOptions for details)\n");
  }
  printf("testing Reporter................... ");
  if (!testReporterImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testReporter for details)\n");
  }
  printf("testing Vector..................... ");
  if (!testVectorImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testVector for details)\n");
  }

  // Return
  return 0;

} // end main
