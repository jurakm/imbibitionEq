/*
 * tableXY.hh
 *
 *  Created on: 7. lis 2015.
 *      Author: jurak
 */

#ifndef SRC_TABLEXY_HH_
#define SRC_TABLEXY_HH_

#include <iostream>
#include <cassert>

/** \brief A table containing two columns : x and y.
 *
 *
 * @tparams VectorType : a type of dynamic vector. Must support resize(),
 */
template <typename VectorType>
class TableXY{
public:
	using T = typename VectorType::value_type;
	/** Construct table with no data. The size is given explicitely or by default. */
	explicit TableXY(unsigned int dts = 10000) : end_(0), table_size_(dts) {
		resize(table_size_);
	}
	/** Resize the table. */
	void resize(unsigned int N) { assert(end_ < N); x_.resize(N); y_.resize(N); table_size_ = N;}
	/** Remove all elements. Do not alert the size of the table. */
	void clean() { end_ = 0; }
	/** Returns the size of the table. */
	unsigned int size() const { return table_size_; }
	/** Is the table empty? */
	bool empty() const { return end_ == 0; }
	/** Index of the last element in the table or -1 if empty. */
	int last() const { return end_ - 1; }
	/** Index of the place following the last element in the table. */
	unsigned int end() const { return end_; }
    /** push back and resize if necessary. */
	void push_back(T x, T y) {
		if (end_ >= table_size_)
			resize(2 * table_size_);
		x_[ end_ ] = x;
		y_[ end_ ] = y;
		end_++;
	}
    /** Set an element at index. It is assumed index < end_. */
    void set(unsigned int index, T x, T y){
    	assert(index < end_);
    	x_[index] = x;
    	y_[index] = y;
    }
    /** Get x coordinate. */
    T get_x(unsigned int index) const{
    	assert(index < end_);
    	return x_[index];
    }
    /** Get y coordinate. */
    T get_y(unsigned int index) const{
        	assert(index < end_);
        	return y_[index];
    }
    /** Get x-coordinate of the last entry. It is assumed that the table is not empty. */
    T get_last_x() const{
    	assert(end_>0);
    	return x_[end_ -1];
    }
    /** Get y-coordinate of the last entry. It is assumed that the table is not empty. */
    T get_last_y() const{
        	assert(end_> 0);
        	return y_[end_ -1];
        }
    /** Find y-coordinate corresponding to x by linear interpolation. */
    T interpolate_y(T x) const {
    	assert(end_> 0);
    	return interpolate(x_, y_,  x, end_-1);
    }
    /** Find x-coordinate corresponding to y by linear interpolation. Used
     * for calculation of the inverse function. */
    T interpolate_x(T y) const {
    	    assert(end_> 0);
        	return interpolate(y_, x_,  y, end_-1);
        }
    /** Check the ordering of data. Return true if o.k and false otherwise. */
    bool check_order() const{
    	for(unsigned int i = 0; i < end_-1; ++i){
    		if(x_[i] > x_[i+1]) return false;
    		if(y_[i] > y_[i+1]) return false;
    	}
    	return true;
    }
private:
    /** Vector of x-coordinates. */
	VectorType x_;
	 /** Vector of y-coordinates. */
	VectorType y_;
	/** Number of elements stored or index of the position following the last element. */
	unsigned int end_;
	/** Initial size of the table. */
	unsigned int table_size_;

	/** General interpolation routine. */
	double interpolate(VectorType const & X, VectorType const & Y, T x, unsigned int last_index) const {
		const double TOL = 1E-8;
		double x_max = X[last_index];
		if (x > x_max + TOL) {
			std::cerr << "x= " << x << ", x_max = " << x_max << ", idx = " << last_index << std::endl;
			throw std::runtime_error("TableXY::interpolate: value out of bounds, error 541.");
		}
		// Small negative numbers are possible.
	//	if(x < 0.0 and -x < TOL) x = 0.0;
		auto it = std::upper_bound(X.begin(), X.begin() + last_index + 1, x);
		auto index = it - X.begin(); // this is upper bound
		if (index == 0){
	//		std::cout << "x = " << x << ", x_max = " << x_max << ", idx = " << last_index << std::endl;
	//		int bbb; std::cin >> bbb;
	//		return x;
			throw std::logic_error("Internal error 542.");
		}
		if (index > last_index) {
			// There is no bigger element, return largest value.
			// This can happen only on the last element.
			assert(index == last_index + 1);
			return Y[index - 1];
		}
		const double y1 = Y[index - 1];
		const double y2 = Y[index];
		const double x1 = X[index - 1];
		const double x2 = X[index];
		double y = y1;
		if (std::abs(x2 - x1) >= TOL)
			y += (x - x1) * (y2 - y1) / (x2 - x1);
		return y;
	}
};



#endif /* SRC_TABLEXY_HH_ */
