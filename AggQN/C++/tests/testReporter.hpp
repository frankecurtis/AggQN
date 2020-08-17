// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTREPORTER_HPP__
#define __TESTREPORTER_HPP__

#include <fstream>
#include <iostream>

#include "FaRSAReporter.hpp"

using namespace FaRSA;

// Implementation of test
int testReporterImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter r;

  // Declare stream report
  std::shared_ptr<StreamReport> rs(new StreamReport("s", R_SOLVER, R_BASIC));

  // Set stream
  rs->setStream(&std::cout);

  // Check option
  if (option == 1) {
    // Add to reporter
    r.addReport(rs);
  }

  // Declare file report
  std::shared_ptr<FileReport> rf(new FileReport("f", R_SOLVER, R_BASIC));

  // Open file
  rf->open("FaRSA_filereport_SOLVER.txt");

  // Add to reporter
  r.addReport(rf);

  // Add file report directly
  r.addFileReport("g", "FaRSA_filereport_SUBSOLVER.txt", R_SUBSOLVER, R_BASIC);

  // Get report that was added directly
  std::shared_ptr<Report> rg = r.report("g");

  // Print something
  r.printf(R_SOLVER, R_BASIC, "This line should appear in all reports.\n");
  r.printf(R_SUBSOLVER, R_BASIC, "This line should appear in all reports.\n");
  r.printf(R_SOLVER, R_BASIC, "This is a string: %s\n", "FaRSA");
  r.printf(R_SUBSOLVER, R_BASIC, "This is a string: %s\n", "FaRSA");
  r.printf(R_SOLVER, R_BASIC, "This is an integer: %d\n", 1);
  r.printf(R_SUBSOLVER, R_BASIC, "This is an integer: %d\n", 1);
  r.printf(R_SOLVER, R_BASIC, "This is a float: %f\n", 2.3);
  r.printf(R_SUBSOLVER, R_BASIC, "This is a float: %f\n", 2.3);
  r.printf(R_SOLVER, R_BASIC, "This is scientific notation: %e\n", 4.56);
  r.printf(R_SUBSOLVER, R_BASIC, "This is scientific notation: %e\n", 4.56);
  r.printf(R_SOLVER, R_BASIC, "Here are all again: %s, %d, %f, %e\n", "FaRSA", 1, 2.3, 4.56);
  r.printf(R_SUBSOLVER, R_BASIC, "Here are all again: %s, %d, %f, %e\n", "FaRSA", 1, 2.3, 4.56);

  // Set types and levels
  rs->setTypeAndLevel(R_SOLVER, R_PER_ITERATION);
  rf->setTypeAndLevel(R_SOLVER, R_PER_INNER_ITERATION);
  rg->setTypeAndLevel(R_SUBSOLVER, R_PER_ITERATION);

  // Print stuff
  r.printf(R_SOLVER, R_BASIC, "SOLVER ALWAYS\n");
  r.printf(R_SOLVER, R_PER_ITERATION, "SOLVER PER ITERATION\n");
  r.printf(R_SOLVER, R_PER_INNER_ITERATION, "SOLVER PER INNER ITERATION\n");
  r.printf(R_SUBSOLVER, R_BASIC, "SUBSOLVER ALWAYS\n");
  r.printf(R_SUBSOLVER, R_PER_ITERATION, "SUBSOLVER PER ITERATION\n");
  r.printf(R_SUBSOLVER, R_PER_INNER_ITERATION, "SUBSOLVER PER INNER ITERATION\n");

  // Delete reports
  r.deleteReports();

  // Read FaRSA_filereport_SOLVER.txt and check values
  std::ifstream infile("FaRSA_filereport_SOLVER.txt");
  std::string line;
  std::getline(infile, line);
  if (line.compare("This line should appear in all reports.") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is a string: FaRSA") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is an integer: 1") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is a float: 2.300000") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is scientific notation: 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("Here are all again: FaRSA, 1, 2.300000, 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("SOLVER ALWAYS") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("SOLVER PER ITERATION") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("SOLVER PER INNER ITERATION") != 0) {
    result = 1;
  }

  // Read FaRSA_filereport_SUBSOLVER.txt and check values
  std::ifstream infile2("FaRSA_filereport_SUBSOLVER.txt");
  std::string line2;
  std::getline(infile2, line2);
  if (line2.compare("This line should appear in all reports.") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is a string: FaRSA") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is an integer: 1") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is a float: 2.300000") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is scientific notation: 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("Here are all again: FaRSA, 1, 2.300000, 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("SUBSOLVER PER ITERATION") != 0) {
    result = 1;
  }

  // Delete files
  remove("FaRSA_filereport_SOLVER.txt");
  remove("FaRSA_filereport_SUBSOLVER.txt");

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      printf("TEST WAS SUCCESSFUL.\n");
    }
    else {
      printf("TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testReporterImplementation

#endif /* __TESTREPORTER_HPP__ */
