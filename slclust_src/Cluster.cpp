/*
 * Cluster.cpp
 *
 *  Created on: 6 août 2010
 *      Author: amadoui
 */

#include "Cluster.h"

using namespace std;

bool sortPoint(Point* a, Point* b) { return (a->x() < b->x()); }
void Cluster::sort() { std::sort(_p_vect.begin(), _p_vect.end(), sortPoint); }

bool Cluster::same_cluster(Point* p, int x_lim, int y_lim, int debug) {
  bool ended = false;
  ostringstream oss;
  for(vector<Point*>::reverse_iterator it = _p_vect.rbegin() ; it != _p_vect.rend() ; ++it) {
    Point* i = *it;

    int x_diff = abs(p->x() - i->x());
    int y_diff = abs(p->y() - i->y());
    bool same = (x_diff <= x_lim && y_diff <= y_lim + x_diff);
    if(debug) {
      cout << "test all P("<<p->x()<<","<<p->y()<<") and i("<<i->x()<<","<<i->y()<<") => x_diff=" << x_diff << " y_diff=" << y_diff << " same=" << same << endl;
    }
     //if(x_diff > x_lim) return false; // Autre choix (p->x() - i->x()) > x_lim
    if((p->x() - i->x()) > x_lim) return false; 
    // if(x_diff > x_lim)  {
    //   if(!ended) {
    // 	ended=true;
    // 	oss << "ENDED P("<<p->x()<<","<<p->y()<<") and i("<<i->x()<<","<<i->y()<<") => x_diff=" << x_diff << " y_diff=" << y_diff << " same=" << same << " ended=" << ended << endl;
    //   }
    // }
    if(x_diff <= x_lim && y_diff <= y_lim + x_diff) { 
      // if(ended) {
      // 	cout << oss.str() << "FIND AFTER ENDED P("<<p->x()<<","<<p->y()<<") and i("<<i->x()<<","<<i->y()<<") => x_diff=" << x_diff << " y_diff=" << y_diff << " same=" << same << endl << endl;
      // }
      return true;
    }
  }
  return false;
}

bool Cluster::same_cluster_onLastPoint(Point* p, int x_lim, int y_lim, int debug) {
  for(vector<Point*>::reverse_iterator it = _p_vect.rbegin() ; it != _p_vect.rend() ; ++it) {
    Point* i = *it;
    
    int x_diff = abs(p->x() - i->x());
    int y_diff = abs(p->y() - i->y());
    bool same = (x_diff <= x_lim && y_diff <= y_lim + x_diff);
    if(debug) {
      cout << "test lastPoint P("<<p->x()<<","<<p->y()<<") from " << this->id() << "(size=" << this->size() << ") and i("<<i->x()<<","<<i->y()<<") => x_diff=" << x_diff << " y_diff=" << y_diff << " same=" << same << endl;
    }
    return (x_diff <= x_lim && y_diff <= y_lim + x_diff);
  }
}

bool Cluster::same_cluster(Cluster* c, int x_lim, int y_lim, bool lastPoint) {
  Cluster *first, *second;
  if(c->x_min() > _x_min) {
    first = this;
    second = c;
  } else {
    first = c;
    second = this; 
  }
  if( second->x_min() - first->x_max() > x_lim) return false;
  int debug = 0;
  // if(first->id()==17284 || second->id()==17284) { 
  //   debug=1; 
  //   cout << "cluster " << first->id() << " size=" << first->size() << " xmin=" << first->x_min() << " xmax=" << first->x_max() << " ymin=" << first->y_min() << " ymax=" << first->y_max() << endl; 
  //   cout << "cluster " << second->id() << " size=" << second->size() << " xmin=" << second->x_min() << " xmax=" << second->x_max() << " ymin=" << second->y_min() << " ymax=" << second->y_max() << endl; 
  //   cout << "first test ok !" << endl; 
  // }
  for(vector<Point*>::iterator it = (c->_p_vect).begin() ; it != (c->_p_vect).end() ; it++) {
    bool same = lastPoint ? same_cluster_onLastPoint(*it, x_lim, y_lim, debug) : same_cluster(*it, x_lim, y_lim, debug);
    if(same) return true;
  }
  return false;
}

bool Cluster::compatible_cluster_on_x(Cluster* c, int x_lim) {
  Cluster *first, *second;
  if(c->x_min() > _x_min) {
    first = this;
    second = c;
  } else {
    first = c;
    second = this;
  }
  if( second->x_min() - first->x_max() > x_lim ) return false;
  return true;
}

void Cluster::addPoint(Point* p) {
  if(p->x() < _x_min) _x_min = p->x();
  if(p->x() > _x_max) _x_max = p->x();
  if(p->y() < _y_min) _y_min = p->y();
  if(p->y() > _y_max) _y_max = p->y();
  _p_vect.push_back(p);
  //this->sort();
}

void Cluster::fuseCluster(Cluster* c) {
  if(c->x_min() < _x_min) _x_min = c->x_min();
  if(c->x_max() > _x_max) _x_max = c->x_max();
  if(c->y_min() < _y_min) _y_min = c->y_min();
  if(c->y_max() > _y_max) _y_max = c->y_max();
  for(vector<Point*>::iterator it = (c->_p_vect).begin() ; it != (c->_p_vect).end() ; it++) {
    _p_vect.push_back(*it);
    //this->sort();
  }
}

void Cluster::getPointList(ostringstream& oss) {
  for(vector<Point*>::iterator it = _p_vect.begin() ; it != _p_vect.end() ; ++it) {
    Point* p = *it;
    oss << "(" << p->x() << "," << p->y() << ")";
  }
}

void Cluster::getMean_var(int& mean, double& var) {
  int sum = 0;
  for(vector<Point*>::iterator it = _p_vect.begin() ; it != _p_vect.end() ; ++it) {
    Point* p = *it;
    sum += (p->y() - p->x());
  }
  mean = sum / this->size();
  int sumsq = 0;
  for(vector<Point*>::iterator it = _p_vect.begin() ; it != _p_vect.end() ; ++it) {
    Point* p = *it;
    sumsq += ( p->y() - p->x() - mean) * ( p->y() - p->x() - mean);
  }
  var = sqrt(sumsq / this->size());
}

bool Cluster::existPoint(int x, int y) {
  for(vector<Point*>::iterator it = _p_vect.begin() ; it != _p_vect.end() ; ++it) {
    Point* p = *it;
    if(x==p->x() && y==p->y()) return true;
  }
  return false;
}
