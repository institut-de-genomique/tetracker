/*
 * Point.h
 *
 *  Created on: 6 août 2010
 *      Author: amadoui
 */

#ifndef POINT_H_
#define POINT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>


class Point {
public:
	protected:
	int _x;
	int _y;
	int _id;

	public:
	/* Constructors and Destructors*/
	Point(int x, int y, int id) {
		_x = x;
		_y = y;
		_id = id;
		}
	/* Accessors */
	int x() const { return _x; }
	int y() const { return _y; }
	int id() const { return _id; }
};



#endif /* POINT_H_ */
