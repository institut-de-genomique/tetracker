/*
 * Cluster.h
 *
 *  Created on: 6 août 2010
 *      Author: amadoui
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <vector>
#include <algorithm> // sort
#include <cstdlib> //abs function
#include <cmath>

#include "Point.h"

using namespace std;

class Cluster {
 public:
 
 protected:
  vector<Point*> _p_vect;
  long _x_min;
  long _x_max;
  long _y_min;
  long _y_max;
  int _id;

 public:
  /* Constructors and Destructors*/
  Cluster(int id) {
    _id = id;
    _x_min = 1000000000000;
    _x_max = 0;
    _y_min = 1000000000000;
    _y_max = 0;
  }

  /* Accessors */
  int id() const { return _id; }
  void id(int i) { _id = i; }
  long x_min() const { return _x_min; }
  long x_max() const { return _x_max; }
  long y_min() const { return _y_min; }
  long y_max() const { return _y_max; }
  long size() const { return _p_vect.size(); }

  void sort();

  bool same_cluster(Point*, int, int, int);
  bool same_cluster_onLastPoint(Point*, int, int, int);
  bool same_cluster(Cluster*, int, int, bool);
  bool compatible_cluster_on_x(Cluster*, int);
  void addPoint(Point*);
  void fuseCluster(Cluster*);
  void getPointList(ostringstream&);
  void getMean_var(int&, double&);
  bool existPoint(int, int);
};



#endif /* CLUSTER_H_ */
