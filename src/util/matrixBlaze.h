////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_MATRIXBLAZE_H
#define UTIL_MATRIXBLAZE_H

#include <blaze/Blaze.h>

namespace util {

/*! @brief Blaze: Definition of vector */
typedef blaze::DynamicVector<double, blaze::columnVector> VectorXi;
/*! @brief Blaze: Definition of n x m matrix */
typedef blaze::DynamicMatrix<double> Matrixij;
/*! @brief Blaze: Definition of n x n symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> SymMatrixij;
/*! @brief Blaze: Definition of 3 x 3 matrix */
typedef blaze::StaticMatrix<double, 3UL, 3UL> Matrix33;
/*! @brief Blaze: Definition of 3 x 3 symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::StaticMatrix<double, 3UL, 3UL>>
    SymMatrix33;
/*! @brief Blaze: Definition of 3D vector */
typedef blaze::StaticVector<double, 3UL> Vector3;
/*! @brief Blaze: Identity matrix */
typedef blaze::IdentityMatrix<double> IdentityMatrix;

/*! @brief Blaze: Definition of vector */
typedef blaze::DynamicVector<float, blaze::columnVector> VectorFXi;
/*! @brief Blaze: Definition of n x m matrix */
typedef blaze::DynamicMatrix<float> MatrixFij;
/*! @brief Blaze: Definition of n x n symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<float>> SymMatrixFij;
/*! @brief Blaze: Definition of 3 x 3 matrix */
typedef blaze::StaticMatrix<float, 3UL, 3UL> MatrixF33;
/*! @brief Blaze: Definition of 3 x 3 symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::StaticMatrix<float, 3UL, 3UL>>
    SymMatrixF33;
/*! @brief Blaze: Definition of 3D vector */
typedef blaze::StaticVector<float, 3UL> VectorF3;
/*! @brief Blaze: Identity matrix */
typedef blaze::IdentityMatrix<float> IdentityMatrixF;

} // namespace util

#endif // UTIL_MATRIXBLAZE_H
