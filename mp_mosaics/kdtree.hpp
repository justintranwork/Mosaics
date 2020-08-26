/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (curDim >= Dim || curDim < 0) {
      return false;
    }
    if (first[curDim] == second[curDim]) {
      return first < second;
    }
    return first[curDim] < second[curDim];
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
    int potentialDist = 0;
    int bestDist = 0;

    for(int i = 0; i < Dim; i++) {
      potentialDist += ((target[i] - potential[i]) * (target[i] - potential[i]));
      bestDist += ((target[i] - currentBest[i]) * (target[i] - currentBest[i]));
    }
    if (bestDist == potentialDist) {
      return potential < currentBest;
    }
    return potentialDist < bestDist;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{

    /**
     * @todo Implement this function!
     */
    int size = 0;
    vector<Point<Dim>> points_;
    points_.assign(newPoints.begin(), newPoints.end());
    root = makeTree(points_, 0, 0, points_.size() - 1);
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::makeTree(vector<Point<Dim>>& points_, int dimension, unsigned left, unsigned right) {
  	//Check edge cases
    if(points_.empty()||left<0||right>=points_.size()||right<0||left>=points_.size()){
      return NULL; 
    } 
  	if(left>right) { 
      return NULL; 
    }
  	unsigned median_idx = (left+right)/2;
  	KDTreeNode* subroot_ = new KDTreeNode(quickSelect(points_,dimension%Dim,left,right,median_idx));
  	size+=1;
  	dimension++;
  	subroot_->left = makeTree(points_,dimension,left,median_idx-1);  	
    subroot_->right = makeTree(points_,dimension,median_idx+1,right);
    return subroot_;
}

template <int Dim>
Point<Dim>& KDTree<Dim>::quickSelect(vector<Point<Dim>>& list, int dimension, unsigned left, unsigned right, unsigned k) {
  	if(left == right) return list[left];
  	unsigned pivotIndex = k;
  	pivotIndex = quickSelect_position(list,dimension, left, right, pivotIndex);
  	if(k == pivotIndex) return list[k];
  	else if(k < pivotIndex)
    	return quickSelect(list, dimension, left, pivotIndex-1, k);
  	else return quickSelect(list, dimension, pivotIndex+1, right, k);
}

template <int Dim>
unsigned KDTree<Dim>::quickSelect_position(vector<Point<Dim>>& list, int dimension, unsigned left, unsigned right, unsigned pivotIndex) {
  	Point<Dim> pivotValue = list[pivotIndex];
  	Point<Dim> temp = list[pivotIndex];
  	list[pivotIndex] = list[right];
  	list[right] = temp;
  	unsigned storeIndex = left;
  	for(unsigned i = left; i < right; i++){
    	if(smallerDimVal(list[i], pivotValue, dimension)){
      		temp = list[storeIndex];
      		list[storeIndex] = list[i];
      		list[i] = temp;
      		storeIndex++;
    	}
  	}
  	temp = list[storeIndex];
  	list[storeIndex] = list[right];
  	list[right] = temp;
  	return storeIndex;
}

//Helper function implementation
template <int Dim>
void KDTree<Dim>::copy(KDTreeNode * copyOnto, KDTreeNode * original) {
    copyOnto = new KDTreeNode();
    copyOnto->point = original->point;
    copy(copyOnto->left, original->left);
    copy(copyOnto->right, original->right);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  copy(this->root, other->root);
  size = other.size;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  if (this != NULL) {
    obliterate(this->root);
  }
  copy(this->root, rhs->root);
  size = rhs.size;
  return *this;
}


template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
  obliterate(this->root);
}

//Helper function implementation
template <int Dim>
void KDTree<Dim>::obliterate(KDTreeNode * subroot) {
  if (subroot == NULL) {
    return;
  }
  obliterate(subroot->left);
  obliterate(subroot->right);
  delete subroot;
  return;
}



template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    return findNearestNeighbor(root, query, 0);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(KDTreeNode * subroot, const Point<Dim> & query, int dimension) const {
  Point<Dim> currentBest = subroot->point;
	bool flag;
	if (subroot->left == NULL && subroot->right == NULL) {
    return subroot->point;
  }

	if (smallerDimVal(query, subroot->point, dimension)) {
		if (subroot->left == NULL) {
			currentBest = findNearestNeighbor(subroot->right, query, (dimension + 1) % Dim);
    } else {
			currentBest = findNearestNeighbor(subroot->left, query, (dimension + 1) % Dim);
    }
		flag = true;
	}	else {
		if (subroot->right == NULL) {
			currentBest = findNearestNeighbor(subroot->left, query, (dimension + 1) % Dim);
    } else {
			currentBest = findNearestNeighbor(subroot->right, query, (dimension + 1) % Dim);
    }
		flag = false;
	}

	if (shouldReplace(query, currentBest, subroot->point)) {
     currentBest = subroot->point;
  }
	
  double bestDist = 0;

	for (int i = 0; i < Dim; i++) {
		bestDist += (query[i] - currentBest[i]) * (query[i] - currentBest[i]);
	}
	
  double potDist = subroot->point[dimension] - query[dimension];
	potDist = potDist * potDist;

	if (potDist <= bestDist) {
		KDTreeNode * need_to_check = flag ? subroot->right : subroot->left;
		if (need_to_check != NULL) {
			Point<Dim> otherBest = findNearestNeighbor(need_to_check, query, (dimension + 1) % Dim);
			if (shouldReplace(query, currentBest, otherBest)) {
        currentBest = otherBest;
      }
		}
	}
	return currentBest;
}

