//============================================================================
// Name        : SVD.cpp
// Author      : George Rokos
// Description : SVD solver for 2x2 linear systems - definitions
//============================================================================

#ifndef SVD2X2_HPP_
#define SVD2X2_HPP_

void svd_solve_2x2(int vid, const double A[4], double p[2], const double q[2]);
void svd_init(int nelements);
void svd_teardown(int nelements);
#endif /* SVD_HPP_ */
