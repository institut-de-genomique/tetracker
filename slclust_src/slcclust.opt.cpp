/*
 * pointslc.cpp
 *
 *  Created on: 6 août 2010
 *      Author: MA Madoui & JM Aury
 */
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <list>
#include <algorithm> // atoi
#include "Point.h"
#include "Cluster.h"
#include "Dist.h"
#include <boost/config.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS> GraphU;
bool sortCluster(Cluster* a, Cluster* b) { return (a->x_min() < b->x_min()); }

void usage();

int main (int ac , char** av){
  //Get pointslc options from command line
  int i;
  int min_size_x = 500, min_size_y = 1000, min_clust_size = 10, read_size = 76, mean_cov = 10, type = 1;
  bool h=false , v=false;
  char *input_file, *output_file;
  while ( (i = getopt(ac, av, "i:l:n:m:t:o:r:hv" )) != -1)
    {
      switch (i)
	{
	case 'i':
	  input_file = optarg;
	  break;
	case 'l':
	  min_size_x = atoi(optarg);
	  break;
	case 'n':
	  min_size_y = atoi(optarg);
	  break;
	case 'm':
	  min_clust_size = atoi(optarg);
	  break;
	case 'o':
	  output_file = optarg;
	  break;
	case 't':
	  type = atoi(optarg);
	  break;
	case 'r':
	  read_size = atoi (optarg);
	  break;
	case 'h':
	  h=true;
	  break;
	case 'v':
	  v = true;
	  break;
	default :
	  abort();
	}
    }
  if (v==true){
    cerr << "** verbose mode **" << endl;
  }
  if (h == true){
    usage();
    return (0);
  }

  //Load Points
  ifstream file( input_file, ios::in);
  list<Cluster*> c_list;
  if(file) {
    int n_point = 0;
    while (!file.eof()){
      int x=-1, y=-1;
      string content;
      file >> x >> y >> content;
      if(x==-1 && y==-1) continue;
      Point* new_point = new Point(x, y, n_point++);
      Cluster* clust = new Cluster(n_point);
      clust->addPoint(new_point);
      c_list.push_back(clust);
    }
    file.close();
  }
  else{
    cerr << "Impossible d'ouvrir le fichier !" << endl;
    usage();
    return (0);
  }
  //sort(c_list.begin(), c_list.end(), sortCluster);
  if(v) cerr << "Number of points= "<< c_list.size() << endl;

  // for(list<Cluster*>::iterator itI = c_list.begin() ; itI != c_list.end() ; itI++) {
  //   Cluster* ci = *itI;
  //   if(ci->existPoint(15094157,15099854))
  //     cout << "found in cluster " << ci->id() << " size=" << ci->size() << " xmin=" << ci->x_min() << " xmax=" << ci->x_max() << " ymin=" << ci->y_min() << " ymax=" << ci->y_max() << endl;
  //   if(ci->existPoint(15094237,15101995))
  //     cout << "found in cluster " << ci->id() << " size=" << ci->size() << " xmin=" << ci->x_min() << " xmax=" << ci->x_max() << " ymin=" << ci->y_min() << " ymax=" << ci->y_max() << endl;
  // }
  
  // Fast but not accurate merging of clusters
  for(list<Cluster*>::iterator itI = c_list.begin() ; itI != c_list.end() ; itI++) { 
    Cluster* ci = *itI;
    list<Cluster*>::iterator itJ = itI;
    itJ++;
    for(; itJ != c_list.end() ; itJ++) {
      Cluster* cj = *itJ;
      if(ci->same_cluster(cj, min_size_x, min_size_y, true)) {
	ci->fuseCluster(cj);
  	itJ = c_list.erase(itJ);
      } else {
	if(!ci->compatible_cluster_on_x(cj, min_size_x)) break;
      }
    }
  }
  if(v) cerr << "Number of clusters= "<< c_list.size() << endl;
  
  // Create hash table based on identifiant and find the biggest cluster
  map<int, Cluster*> c_map;
  int max_size=0;
  Cluster* max;
  int new_id = 0;
  for(list<Cluster*>::iterator itI = c_list.begin() ; itI != c_list.end() ; itI++) {
    Cluster* current = *itI;
    current->sort();
    current->id(new_id++);
    c_map[(*itI)->id()] = *itI;
    if((*itI)->size() > max_size) {
      max = *itI;
      max_size = max->size();
    }
  }
  if(v) cerr << "Largest Cluster= "<< max->id() << " size=" << max->size() << " xmin=" << max->x_min() << " xmax=" << max->x_max()<< endl;

  //Get distances vector
  vector<Dist*> d_list;
  int unsigned n_dist = 0;
  int cpt=0;
  for(list<Cluster*>::iterator itI = c_list.begin() ; itI != c_list.end() ; itI++) {
    Cluster* ci = *itI;
    cpt++;
    if(v) if(cpt%500==0) cerr << "\r" << cpt << " of " << c_list.size() << " ; Number of distance added= " << n_dist << flush;
    //if(ci->size() > 5000) continue;
    list<Cluster*>::iterator itJ = itI;
    itJ++;
    for(; itJ != c_list.end() ; itJ++) {
      Cluster* cj = *itJ;
      //if(cj->size() > 5000) continue;
      if(ci->same_cluster(cj, min_size_x, min_size_y, false)) {
	Dist* new_dist = new Dist(ci->id(), cj->id(), ++n_dist);
	d_list.push_back(new_dist);
      }
    }
  }
  if(v) cerr << "\r" << cpt << " of " << c_list.size() << " ; Number of distance added= " << n_dist << endl;

  //Build graph
  GraphU gU(c_list.size());
  if(v) cerr << "Build graph with #vertices= " << num_vertices(gU) << " and #edges= " << d_list.size() << endl;
  for( vector<Dist*>::iterator it = d_list.begin() ; it != d_list.end() ; it++ ) add_edge( (*it)->i(), (*it)->j(), gU );
  std::vector<int> component(num_vertices(gU));
  int _nb_connected_components = connected_components(gU, &component[0]);
  if(v) cerr << "Number of connected components= " << _nb_connected_components << endl;
  vector<list<int> > cc(_nb_connected_components);
  for (int i=0 ; i < component.size() ; i++) cc[component[i]].push_back(i);
  
  //Get clusters
  ofstream outfile ( output_file, ios::out );
  if (!outfile) cerr << "not open " << endl;
  for(vector<list<int> >::iterator it = cc.begin(); it != cc.end() ; it++) {
    list<int> l = *it;
    Cluster *first = NULL;
    for(list<int>::iterator itL = l.begin() ; itL != l.end() ; itL++) {
      int current_id = *itL;
      map<int, Cluster*>::iterator itMap = c_map.find(current_id);
      if(itMap == c_map.end()) { cout << "Empty cluster " << current_id << endl; }
      Cluster *current = itMap->second;
      assert(current);
      if(!first) { first = current; }
      else {
	first->fuseCluster(current);
	// erase current cluster in c_list
	for(list<Cluster*>::iterator it = c_list.begin(); it != c_list.end() ; it++) {
	  if((*it)->id() == current_id) { 
	    it = c_list.erase(it);
	  }
	}
      }
    }
  }
  if(v) cerr << "Number of clusters= "<< c_list.size() << endl;
  
  list<Cluster*> c_new_list;
  for(list<Cluster*>::iterator itI = c_list.begin() ; itI != c_list.end() ; itI++) {
    Cluster* current = *itI;
    if(current->size() > min_clust_size) c_new_list.push_back(current);
  }
  if(v) cerr << "Final number of clusters (size>" << min_clust_size << ")= "<< c_new_list.size() << endl;
  
  int i=0;
  int sum = 0;
  for(list<Cluster*>::iterator it = c_new_list.begin() ; it != c_new_list.end() ; it++) {
    ostringstream oss, oss2;
    int sum = 0;
    Cluster* current = *it;
    current->getPointList(oss2);
    oss << current->x_min() << "\t" << current->x_max() << "\t" << current->y_min() << "\t" << current->y_max();
    //int nb_reads = 2 * current->size();
    int nb_reads = current->size();
    if ( type == 0  ) {
      int mean;
      double var;
      current->getMean_var(mean, var);
      oss << "\t" << mean << "\t" << var;
      outfile << ++i << "\t"<< nb_reads <<"\t"<< oss.str()<<"\t"<< oss2.str() << endl;
    }
    if ( type == 1 ){
      outfile << ++i << "\t"<< nb_reads <<"\t"<< oss.str() <<"\t"<< oss2.str() << endl;
    }
  }

  outfile.close();
}



void usage() {
  cerr << "--------------------------------------------------------------------------------------------" << endl;
  cerr << "slclust - run single linkage clustering (SLC) for points cloud." << endl << endl;
  cerr << "Usage : pointslc -i <infile> -o <outfile> [-options]" << endl;
  cerr << "	-i <file>  	: input tabular file with two column (x and y)" << endl;
  cerr << "	-o <file>  	: output tabular file" << endl;
  cerr << "	-l <int>   	: minimum distance between two x coordinates (default 500)" << endl;
  cerr << "	-n <int> 	: minimum distance between two y coordinates (default 1000)" << endl;
  cerr << "	-m <int>   	: minimum count of points in a cluster (default 10)" << endl;
  cerr << "	-r <int>	: reads size (default 76)" << endl;
  cerr << "	-c <int>	: mean coverage (default 30)" << endl;
  cerr << "	-t <int>	: type of discordant mate pairs,  1 for intrachromosomal, 2 for interchromosomal (default is 1)" << endl;
  cerr << "	-v        	: verbose mode" << endl;
  cerr << "	-h        	: this help" << endl;
  cerr << "--------------------------------------------------------------------------------------------" << endl;
  exit(1);
}
