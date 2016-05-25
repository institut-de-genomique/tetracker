/*
 * dist.h
 *
 *  Created on: 9 août 2010
 *      Author: amadoui
 */

#ifndef DIST_H_
#define DIST_H_

#include <math.h>

#include "Point.h"

class Dist {
 public:
 
 protected:
  //Point* _i;
  //Point* _j;
  unsigned int _i;
  unsigned int _j;
  //double _dist;
  int _id;
  
 public:
  /* Constructors and Destructors*/
  //Dist(Point* p1, Point* p2, int id) {
  Dist(int id1, int id2, int id) {  
    _i = id1;
    _j = id2;
    //_dist = sqrt( pow( abs(p1->x() - p2->x()), 2) + pow( abs(p1->y() - p2->y()), 2) );//euclidian distance
    _id = id;
  }

  /* Accessors */
  //Point* i() const { return _i; }
  unsigned int i() const { return _i; }
  //Point* j() const { return _j; }
  unsigned int j() const { return _j; }
  //double dist() const {return _dist;}
};

#endif /* DIST_H_ */
