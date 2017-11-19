//to compile g++ -I/usr/include/eigen3/ tt.cpp -o tt
#include <istream>
#include <string>
//#include <sstream>
//#include <vector>
#include <stdio.h>
//#include <stdlib.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
//#include <eigen3>
//#include <fstream>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

int main()
{
	int nrows=4;//130;
	int ncols=4;
	int ncols2=9;
VectorXf betai(9);
betai << 0.003847799,0.001000000,0.001000000,0.077725734,0.001000000,6.862597279,4.714029223,2.095145620,4.58666442;
//
//MatrixXf mymatrix(rows,columns); 
MatrixXf X = MatrixXf::Zero(nrows,ncols2);
X << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36;
MatrixXf X1 = MatrixXf::Zero(nrows,ncols2);
MatrixXf X2 = MatrixXf::Zero(nrows,ncols2);
MatrixXf R1 = MatrixXf::Zero(nrows,nrows);
MatrixXf R2 = MatrixXf::Zero(nrows,nrows);
MatrixXf ones = MatrixXf::Ones(1,X2.rows());
MatrixXf ones2 = MatrixXf::Ones(1,X1.rows());
MatrixXf S1, S2, a1, a2; 
X1=X;
X2=X;

MatrixXf pdf = MatrixXf::Zero(ncols2,ncols2);
pdf = betai.asDiagonal();//turn vector to digaonal matrix
R1 = X1*pdf*X1.transpose();
R2 = X2*pdf*X2.transpose();
//gap=R1.diagonal(); //extract diagonal of a matrix
//A = kroneckerProduct(R1.diagonal(),ones).eval();
S1 = R1.diagonal()*ones;
S2 = R2.diagonal()*ones2;
a1 = X1*pdf*X2.transpose(); //t(tx1)%*%pdf%*%tx2
a2 = X2*pdf*X1.transpose();
MatrixXf res = (a1.transpose()+a2-S1.transpose()-S2);
//ofstream a_file ( "example.txt" );
// Outputs to example.txt through a_file
//a_file << X1;
//a_file.close();
//R1.determinant();
//Array3d v(0,-1485.78, -5943.14);
//cout << v.exp() << endl;
ArrayXXf res2=res.array();
//cout << "betai = " << endl << R1.bottomRightCorner(5,5) << endl;
//cout << "betai = " << endl << X1.transpose().bottomRightCorner(9,5) << endl;
//cout << S1.bottomRightCorner(3,2) << endl;
cout << res << endl;
cout << res.array().exp() << endl;

}
