#pragma once

#include "glm/glm.hpp"
#include "glm/gtc/constants.hpp"

class Transform
{
public:

	double* position;
	double* rotation;
	double* scale;
	double* matrix;

	Transform()
		: position(new double[3]), rotation(new double[3]), scale(new double[3]), matrix(new double[16])
	{
		for (int i = 0; i < 3; ++i) {
			position[i] = 0.0f;
			rotation[i] = 0.0f;
			scale[i] = 1.0f;
		}

		identity(matrix);
	}

	~Transform()
	{
		delete[] position;
		position = nullptr;

		delete[] rotation;
		rotation = nullptr;

		delete[] scale;
		scale = nullptr;

		delete[] matrix;
		matrix = nullptr;
	}

	// Set position vector
	void setPosition(const double &x, const double &y, const double &z)
	{
		position[0] = x;
		position[1] = y;
		position[2] = z;

		this->composeMatrix();
	}

	// Set rotation vector
	void setRotation(const double &x, const double &y, const double &z)
	{
		rotation[0] = x;
		rotation[1] = y;
		rotation[2] = z;

		this->composeMatrix();
	}

	void setScale(const double &x, const double &y, const double &z)
	{
		scale[0] = x;
		scale[1] = y;
		scale[2] = z;

		this->composeMatrix();
	}

	// Creates an identity matrix
	static void identity(double M[])
	{
		M[0] = 1.0; M[4] = 0.0; M[8] = 0.0;  M[12] = 0.0;
		M[1] = 0.0; M[5] = 1.0; M[9] = 0.0;  M[13] = 0.0;
		M[2] = 0.0; M[6] = 0.0; M[10] = 1.0;  M[14] = 0.0;
		M[3] = 0.0; M[7] = 0.0; M[11] = 0.0;  M[15] = 1.0;
	}

	// Creates a translation matrix
	static void translate(double M[], const double &x, const double &y, const double &z)
	{
		M[0] = 1.0; M[4] = 0.0; M[8] = 0.0;  M[12] = x;
		M[1] = 0.0; M[5] = 1.0; M[9] = 0.0;  M[13] = y;
		M[2] = 0.0; M[6] = 0.0; M[10] = 1.0; M[14] = z;
		M[3] = 0.0; M[7] = 0.0; M[11] = 0.0; M[15] = 1.0;
	}

	// Creates a scale matrix
	static void scaleMatrix(double M[], const double &x, const double &y, const double &z)
	{
		M[0] = x;	M[4] = 0.0; M[8] = 0.0;  M[12] = 0.0;
		M[1] = 0.0; M[5] = y;	M[9] = 0.0;  M[13] = 0.0;
		M[2] = 0.0; M[6] = 0.0; M[10] = z;   M[14] = 0.0;
		M[3] = 0.0; M[7] = 0.0; M[11] = 0.0; M[15] = 1.0;
	}

	// Creates a rotation matrix in X-direction
	static void rotateX(double M[], const double &angle)
	{
		M[0] = 1.0; M[4] = 0.0;        M[8] = 0.0;         M[12] = 0.0;
		M[1] = 0.0; M[5] = cos(angle); M[9] = -sin(angle); M[13] = 0.0;
		M[2] = 0.0; M[6] = sin(angle); M[10] = cos(angle); M[14] = 0.0;
		M[3] = 0.0; M[7] = 0.0;        M[11] = 0.0;        M[15] = 1.0;
	}

	// Creates a rotation matrix in Y-direction
	static void rotateY(double M[], const double &angle)
	{
		M[0] = cos(angle);  M[4] = 0.0; M[8] = sin(angle);  M[12] = 0.0;
		M[1] = 0.0;         M[5] = 1.0; M[9] = 0.0;         M[13] = 0.0;
		M[2] = -sin(angle); M[6] = 0.0; M[10] = cos(angle); M[14] = 0.0;
		M[3] = 0.0;         M[7] = 0.0; M[11] = 0.0;        M[15] = 1.0;
	}

	// Creates a rotation matrix in Z-direction
	static void rotateZ(double M[], const double &angle)
	{
		M[0] = cos(angle); M[4] = -sin(angle); M[8] = 0.0;  M[12] = 0.0;
		M[1] = sin(angle); M[5] = cos(angle);  M[9] = 0.0;  M[13] = 0.0;
		M[2] = 0.0;        M[6] = 0.0;         M[10] = 1.0; M[14] = 0.0;
		M[3] = 0.0;        M[7] = 0.0;         M[11] = 0.0; M[15] = 1.0;
	}

	// Creates a perspective matrix
	static void perspective(double M[], const double &vert_fov, const double &aspect_ratio, const double &z_near, const double &z_far)
	{
		double f = cos(vert_fov / 2) / sin(vert_fov / 2);

		M[0] = f / aspect_ratio; M[4] = 0.0; M[8] = 0.0;                                   M[12] = 0.0;
		M[1] = 0.0;              M[5] = f;   M[9] = 0.0;                                   M[13] = 0.0;
		M[2] = 0.0;              M[6] = 0.0; M[10] = -(z_far + z_near) / (z_far - z_near); M[14] = -(2.0 * z_near*z_far) / (z_far - z_near);
		M[3] = 0.0;              M[7] = 0.0; M[11] = -1.0;                                 M[15] = 0.0;
	}

	// Inverts a given matrix and returns the matrix combination
	static void invertMatrix(double out[], double const a[])
	{
		double a00 = a[0],
			a01 = a[1],
			a02 = a[2],
			a03 = a[3];

		double a10 = a[4],
			a11 = a[5],
			a12 = a[6],
			a13 = a[7];

		double a20 = a[8],
			a21 = a[9],
			a22 = a[10],
			a23 = a[11];

		double a30 = a[12],
			a31 = a[13],
			a32 = a[14],
			a33 = a[15];

		double b00 = a00 * a11 - a01 * a10;
		double b01 = a00 * a12 - a02 * a10;
		double b02 = a00 * a13 - a03 * a10;
		double b03 = a01 * a12 - a02 * a11;
		double b04 = a01 * a13 - a03 * a11;
		double b05 = a02 * a13 - a03 * a12;
		double b06 = a20 * a31 - a21 * a30;
		double b07 = a20 * a32 - a22 * a30;
		double b08 = a20 * a33 - a23 * a30;
		double b09 = a21 * a32 - a22 * a31;
		double b10 = a21 * a33 - a23 * a31;
		double b11 = a22 * a33 - a23 * a32;

		// Calculate the determinant
		double det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

		if (abs(det) > 0.000001)
			det = 1.0 / det;

		out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
		out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
		out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
		out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
		out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
		out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
		out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
		out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
		out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
		out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
		out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
		out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
		out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
		out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
		out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
		out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;
	}

	// Performs a matrix multiplication
	static void multiply(const double M1[], const double M2[], double result[])
	{
		double temp[16];
		int index = 0;
		int i = 0;

		while (i < 16) {
			// Each column
			for (int k = 0; k < 4; ++k) {
				temp[index] = (M1[k] * M2[i]) + (M1[k + 4] * M2[i + 1]) + (M1[k + 8] * M2[i + 2]) + (M1[k + 12] * M2[i + 3]);
				++index;
			}
			i = i + 4; // Jump to next row (go pass 4 indexes)
		}

		// Copy the result to one of the matrices
		for (int n = 0; n < 16; ++n) {
			result[n] = temp[n];
		}
	}

private:
	bool inverted;

	void composeMatrix()
	{
		double dummy[16];

		/*** First, add scale matrix ***/
		scaleMatrix(matrix, scale[0], scale[1], scale[2]);

		/*** Second, add rotation matrix ***/
		// 1. Create the rotation matrix by applying rotation to x, y, z one at the time
		double rotation_matrix[16];
		identity(rotation_matrix);
		rotateZ(dummy, (rotation[2] * glm::pi<double>() / 180.0));
		multiply(dummy, rotation_matrix, rotation_matrix);
		rotateX(dummy, (rotation[0] * glm::pi<double>() / 180.0));
		multiply(dummy, rotation_matrix, rotation_matrix);
		rotateY(dummy, (rotation[1] * glm::pi<double>() / 180.0));
		// 2. Multiply the final rotation matrix to the scale matrix
		multiply(dummy, rotation_matrix, rotation_matrix);
		multiply(rotation_matrix, matrix, matrix);

		/*** Third, translation matrix ***/
		translate(dummy, position[0], position[1], position[2]);
		multiply(dummy, matrix, matrix);
	}
};