// Copyright (C) 2020 Albert Berahas, Frank E. Curtis, Baoyu Zhou
//
// This code is published under the ??? License.
//
// Author(s) : Albert Berahas, Frank E. Curtis, Baoyu Zhou

#ifndef __AGGQNENUMERATIONS_HPP__
#define __AGGQNENUMERATIONS_HPP__

namespace AggQN
{

/** @name Enumerations */
//@{
/**
 * Report type enumerations
 */
enum ReportType
{
  R_SOLVER = 0,
  R_SUBSOLVER
};
/**
 * Report level enumerations
 */
enum ReportLevel
{
  R_BASIC = 0,
  R_PER_ITERATION,
  R_PER_INNER_ITERATION
};
//@}

} // namespace AggQN

#endif /* __AGGQNENUMERATIONS_HPP__ */
