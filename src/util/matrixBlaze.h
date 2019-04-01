// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef UTIL_MATRIXBLAZE_H
#define UTIL_MATRIXBLAZE_H

#include <blaze/Blaze.h>

namespace util {


/*! @brief Definition of vector */
typedef blaze::DynamicVector<double,blaze::columnVector> VectorXi;
/*! @brief Definition of n x m matrix */
typedef blaze::DynamicMatrix<double> Matrixij;
/*! @brief Definition of n x n symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> SymMatrixij;
/*! @brief Definition of 3 x 3 matrix */
typedef blaze::StaticMatrix<double,3UL,3UL> Matrix33;
/*! @brief Definition of 3 x 3 symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::StaticMatrix<double,3UL,3UL>> SymMatrix33;
/*! @brief Definition of 3D vector */
typedef blaze::StaticVector<double,3UL> Vector3;
/*! @brief Identity matrix */
typedef blaze::IdentityMatrix<double> IdentityMatrix;

/*! @brief Definition of vector */
typedef blaze::DynamicVector<float,blaze::columnVector> VectorFXi;
/*! @brief Definition of n x m matrix */
typedef blaze::DynamicMatrix<float> MatrixFij;
/*! @brief Definition of n x n symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<float>> SymMatrixFij;
/*! @brief Definition of 3 x 3 matrix */
typedef blaze::StaticMatrix<float,3UL,3UL> MatrixF33;
/*! @brief Definition of 3 x 3 symmetric matrix */
typedef blaze::SymmetricMatrix<blaze::StaticMatrix<float,3UL,3UL>>
SymMatrixF33;
/*! @brief Definition of 3D vector */
typedef blaze::StaticVector<float,3UL> VectorF3;
/*! @brief Identity matrix */
typedef blaze::IdentityMatrix<float> IdentityMatrixF;

} // namespace util

#endif