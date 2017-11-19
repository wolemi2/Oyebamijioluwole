//to compile g++ -I/usr/include/eigen3/ emu.cpp -o emu
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
//#include <Eigen/Core>
using namespace std;
using namespace Eigen;

//this function computes distance correlation
MatrixXd mycor(const MatrixXd& X1,const MatrixXd& X2,const VectorXd& betai)
{
MatrixXd ones = MatrixXd::Ones(1,X2.rows());
MatrixXd ones2 = MatrixXd::Ones(1,X1.rows());
MatrixXd pdf = betai.asDiagonal();//turn vector to digaonal matrix
MatrixXd R1 = X1*pdf*X1.transpose();
MatrixXd R2 = X2*pdf*X2.transpose();
MatrixXd S1 = R1.diagonal()*ones;
MatrixXd S2 = R2.diagonal()*ones2;
MatrixXd a1 = X1*pdf*X2.transpose(); //t(tx1)%*%pdf%*%tx2
MatrixXd a2 = X2*pdf*X1.transpose();
MatrixXd res = (a1.transpose()+a2-S1.transpose()-S2);
ArrayXXd res2 = res.array();
return(res2.exp());
}
///////////////////////
//main function
int main()
{
	int nrows=4130;
	int ncols=4;
	int ncols2=9;
	double tau0=0.001000;
	double tau=0.5;
VectorXd betai(9);
betai << 0.003847799,0.001000000,0.001000000,0.077725734,0.001000000,6.862597279,4.714029223,2.095145620,4.58666442;
//read training input dat matrix X
MatrixXd X = MatrixXd::Zero(nrows,ncols2);
ifstream fin1 ("./X.txt");
if (fin1.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin1 >> item;
X(row, col) = item;
}
fin1.close();
}
//read in output matrix Y
MatrixXd Y = MatrixXd::Zero(nrows,ncols);
ifstream fin2 ("./Y.txt");
if (fin2.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < ncols; col++)
{
float item;
fin2 >> item;
Y(row, col) = item;
}
fin2.close();
}

//read in new data point for prediction
MatrixXd newdata = MatrixXd::Zero(1,ncols2);
ifstream fin3 ("./input/input.txt");
if (fin3.is_open())
{
for (int row = 0; row < 1; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin3 >> item;
newdata(row, col) = item;
}
fin3.close();
}

//read in Ainv matrix
MatrixXd Ainv = MatrixXd::Zero(nrows,nrows);
ifstream fin4 ("./Ainv.txt");
if (fin4.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < nrows; col++)
{
float item;
fin4 >> item;
Ainv(row, col) = item;
}
fin4.close();
}

//read in iOmega matrix
MatrixXd iOmega = MatrixXd::Zero(ncols2,ncols2);
ifstream fin5 ("./iOmega.txt");
if (fin5.is_open())
{
for (int row = 0; row < ncols2; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin5 >> item;
iOmega(row, col) = item;
}
fin5.close();
}
//read in betahat matrix
MatrixXd betahat = MatrixXd::Zero(ncols2,ncols);
ifstream fin6 ("./betahat.txt");
if (fin6.is_open())
{
for (int row = 0; row < ncols2; row++)
for (int col = 0; col < ncols; col++)
{
float item;
fin6 >> item;
betahat(row, col) = item;
}
fin6.close();
}
//normalization data
MatrixXd s1(1,4);
MatrixXd s2(1,4);
s2 << 2.166465e-15,2.079999e+01,4.336633e-10,4.304684e-10; 
s1 << 2.467027e-15,2.339083e+01,8.968750e-10,6.864521e-11; 
MatrixXd H = X.array().log();//log-tranformed training inputs
MatrixXd H0 = newdata.array().log();//log of test input
int m=H.cols();
int n=Y.rows();
int n2=H0.rows();
MatrixXd A01 = mycor(H0,H,betai);//###cross correlation
MatrixXd A00 = mycor(H0,H0,betai);//##test point correlation
//MatrixXd A = mycor(H,H,betai);##very big matrix precomputed in R
MatrixXd temp = MatrixXd::Identity(A00.rows(), A00.cols());
MatrixXd A000 = A00+((tau0*tau)*temp);
//MatrixXd iOmega = (H.transpose()*Ainv*H).inverse();//can be precomputed in R
//MatrixXd betahat= iOmega*(H.transpose()*Ainv*Y);// can be precomputed in R
MatrixXd mu_star = H0*betahat+A01.transpose()*Ainv*(Y-H*betahat);
MatrixXd r1 = H0-(A01.transpose()*Ainv*H);
MatrixXd c_star = A000-(A01.transpose()*Ainv*A01)+(((r1)*(iOmega)*r1.transpose()));
MatrixXd Sigma = ((Y-H*betahat).transpose()*Ainv*(Y-H*betahat))/(n-m);
MatrixXd rate = (mu_star.array()*s2.array())+s1.array();
//output results
ofstream a_file ("./output/Rate.txt");
a_file << rate;
a_file.close();
cout << betahat << endl;
}
