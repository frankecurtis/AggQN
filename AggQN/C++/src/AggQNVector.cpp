// Copyright (C) 2020 Albert Berahas, Frank E. Curtis, Baoyu Zhou
//
// This code is published under the ??? License.
//
// Author(s) : Albert Berahas, Frank E. Curtis, Baoyu Zhou

#include <cassert>
#include <cmath>
#include <cstdio>

#include "AggQNVector.hpp"

namespace AggQN
{

// Constructor with given length; values initialized to zero
Vector::Vector(int length)
  : length_(length)
{

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length; i++) {
    values_[i] = 0.0;
  }

} // end constructor

// Constructor with given length; values initialized to given value
Vector::Vector(int length,
               double value)
  : length_(length)
{

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length; i++) {
    values_[i] = value;
  }

} // end constructor

// Destructor; values array deleted
Vector::~Vector()
{

  // Delete array
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

} // end destructor

// Print array with given name
void Vector::print(const Reporter* reporter,
                   std::string name) const
{

  // Print elements
  for (int i = 0; i < length_; i++) {
    reporter->printf(R_SOLVER, R_BASIC, "%s[%6d]=%+23.16e\n", name.c_str(), i, values_[i]);
  }

} // end print

// Make new Vector as a copy
std::shared_ptr<Vector> Vector::makeNewCopy() const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy elements
  vector->copy(*this);

  // Return
  return vector;

} // end makeNewCopy

// Make new Vector by adding "scalar1" times this Vector to "scalar2" times other_vector
std::shared_ptr<Vector> Vector::makeNewLinearCombination(double scalar1,
                                                         double scalar2,
                                                         const Vector& other_vector) const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy + add elements
  vector->linearCombination(scalar1, *this, scalar2, other_vector);

  // Return
  return vector;

} // end makeNewLinearCombination

// Set length and initialize values to zero
void Vector::setLength(int length)
{

  // Store length
  length_ = length;

  // Delete previous array, if exists
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length_; i++) {
    values_[i] = 0.0;
  }

} // end setLength

// Set element with given index to given value
void Vector::set(int index,
                 double value)
{

  // Asserts
  assert(index >= 0);
  assert(index < length_);

  // Set value
  values_[index] = value;

} // end set

// Copy elements of other_vector
void Vector::copy(const Vector& other_vector)
{

  // Assert
  assert(length_ == other_vector.length());

  // Copy elements
  for (int i = 0; i < length_; i++) {
    values_[i] = other_vector.values()[i];
  }

} // end copy

// Copy elements of double array
void Vector::copyArray(double* array)
{

  // Copy elements
  for (int i = 0; i < length_; i++) {
    values_[i] = array[i];
  }

} // end copyArray

// Scale elements by given scalar
void Vector::scale(double scalar)
{

  // Check for zero scalar
  if (scalar == 0.0) {

    // Scale elements
    for (int i = 0; i < length_; i++) {
      values_[i] = 0.0;
    }

  } // end if
  else if (scalar != 1.0) {

    // Scale elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar*values_[i];
    }

  } // end else

} // end scale

// Add to this Vector "scalar" times other_vector
void Vector::addScaledVector(double scalar,
                             const Vector& other_vector)
{

  // Assert
  assert(length_ == other_vector.length());

  // Add scaled vector
  for (int i = 0; i < length_; i++) {
    values_[i] += scalar*other_vector.values()[i];
  }

} // end addScaledVector

// Set values as linear combination (scalar1*vector1 + scalar2*vector2)
void Vector::linearCombination(double scalar1,
                               const Vector& vector1,
                               double scalar2,
                               const Vector& vector2)
{

  // Asserts
  assert(length_ == vector1.length());
  assert(length_ == vector2.length());

  // Check for nonzero scalars
  if (scalar1 != 0.0 && scalar2 != 0.0) {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar1*vector1.values()[i] + scalar2*vector2.values()[i];
    }

  } // end if
  else if (scalar1 != 0.0) {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar1*vector1.values()[i];
    }

  } // end else if
  else {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar2*vector2.values()[i];
    }

  } // end else

} // end linearCombination

// Inner product with other_vector
double Vector::innerProduct(const Vector& other_vector) const
{

  // Assert
  assert(length_ == other_vector.length());

  // Compute inner product
  double inner_product = 0.0;
  for (int i = 0; i < length_; i++) {
    inner_product += values_[i] * other_vector.values()[i];
  }

  // Return inner product
  return inner_product;

} // end innerProduct

// Maximum element
double Vector::max() const
{

  // Determine maximum
  double maximum = values_[0];
  for (int i = 1; i < length_; i++) {
    maximum = fmax(maximum, values_[i]);
  }

  // Return maximum
  return maximum;

} // end max

// Minimum element
double Vector::min() const
{

  // Determine minimum
  double minimum = values_[0];
  for (int i = 1; i < length_; i++) {
    minimum = fmin(minimum, values_[i]);
  }

  // Return minimum
  return minimum;

} // end max

// 1-norm
double Vector::norm1() const
{

  // Determine 1-norm
  double norm = 0.0;
  for (int i = 0; i < length_; i++) {
    norm += fabs(values_[i]);
  }

  // Return 1-norm
  return norm;

} // end norm1

// 2-norm
double Vector::norm2() const
{

  // Determine 2-norm squared
  double norm = 0.0;
  for (int i = 0; i < length_; i++) {
    norm += pow(values_[i],2.0);
  }

  // Return 1-norm
  return sqrt(norm);

} // end norm2

// inf-norm
double Vector::normInf() const
{

  // Determine inf-norm
  double norm = fabs(values_[0]);
  for (int i = 0; i < length_; i++) {
    norm = fmax(norm, fabs(values_[0]));
  }

  // Return inf-norm
  return norm;

} // end normInf

} // namespace NonOpt
