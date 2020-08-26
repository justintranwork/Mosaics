/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    //Create a KDTree so that we can compare them
    vector<Point<3>> imageList;
    map<Point<3>, TileImage*> map_;
    for (auto iter = theTiles.begin(); iter != theTiles.end(); iter++) {
    	LUVAPixel pixel_ = iter->getAverageColor();
    	Point<3> point_ = convertToXYZ(pixel_);
    	imageList.push_back(point_);
    	map_[point_] = &*iter;
    }
    KDTree<3> tileKD(imageList);  
    //Create the canvas for the return image


    MosaicCanvas * returnImage = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    // Itreate through the canvas, choosing a best tile based on the average color of each region
    for (int i = 0; i < returnImage->getRows(); i++) {
    	for (int j = 0; j < returnImage->getColumns(); j++) {
    		Point<3> goal_ = convertToXYZ(theSource.getRegionColor(i, j));  // Average color of region
    		Point<3> best_ = tileKD.findNearestNeighbor(goal_);  // Find NN of average color of region
    		TileImage* correct_tile = map_[best_];  // Get the tile that is mapped to by best point
    		returnImage->setTile(i, j, correct_tile);  // Set that tile on the canvas
    	}
    }
    return returnImage;
}
