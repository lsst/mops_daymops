#include "lsst/mops/daymops/linkTracklets/TrackletTree.h"

namespace lsst { namespace mops {



TrackletTree::TrackletTree(const std::vector<MopsDetection> &allDetections,
                           const std::vector<Tracklet> &thisTreeTracklets,
                           double positionalErrorRa, 
                           double positionalErrorDec,
                           unsigned int maxLeafSize,
			   const std::vector<double> &perAxisWidths)
{
    setUpEmptyTree();
    buildFromData(allDetections, thisTreeTracklets, positionalErrorRa,
		  positionalErrorDec, maxLeafSize, perAxisWidths);

}


TrackletTree::TrackletTree(const std::vector<MopsDetection> &allDetections,
                           const std::vector<Tracklet> &thisTreeTracklets,
                           double positionalErrorRa, 
                           double positionalErrorDec,
                           unsigned int maxLeafSize)
{
    setUpEmptyTree();
    std::vector<double> emptyVec;
    buildFromData(allDetections, thisTreeTracklets, positionalErrorRa,
		  positionalErrorDec, maxLeafSize, emptyVec);

}







void TrackletTree::buildFromData(
    const std::vector<MopsDetection> &allDetections,
    const std::vector<Tracklet> &thisTreeTracklets,
    double positionalErrorRa, double positionalErrorDec,
    unsigned int maxLeafSize, 
    const::std::vector<double> &perAxisWidths)
{
    // need to set up fields used by KDTree just the way KDTree would; then 
    // create our set of child TrackletTreesNodes
    myK = 4;


    if (thisTreeTracklets.size() > 0) 
    {

        std::vector<PointAndValue <unsigned int> > parameterizedTracklets;
        std::vector<double> pointsUBounds, pointsLBounds;

        if (maxLeafSize < 1) {
            throw LSST_EXCEPT(BadParameterException, 
       "EE: KDTree: max leaf size must be strictly positive!\n");
        }


        // need to convert tracklets to a parameterized format:

        // (RA_0, Dec_0, RAv, Dec_v) 
        
        // and make a PointAndValue vector.

	// ASSUME all data comes from the same <180 -degree region of sky in both RA and Dec.

        for (uint i = 0; i < thisTreeTracklets.size(); i++) {
            Tracklet myT = thisTreeTracklets.at(i);
            MopsDetection firstDetection = 
                myT.getFirstDetection(allDetections);
            PointAndValue<unsigned int> trackletPav;

            std::vector<double> trackletPoint;
	    const std::vector<double> *raP0Vel = myT.getBestFitFunctionRa();
	    const std::vector<double> *decP0Vel = myT.getBestFitFunctionDec();
            trackletPoint.push_back(raP0Vel->at(0));
            trackletPoint.push_back(decP0Vel->at(0));
            trackletPoint.push_back(raP0Vel->at(1));
            trackletPoint.push_back(decP0Vel->at(1));
            trackletPoint.push_back(myT.getDeltaTime(allDetections));

            trackletPav.setPoint(trackletPoint);
            trackletPav.setValue(myT.getId());

            parameterizedTracklets.push_back(trackletPav);

	    // calculate UBounds, LBounds                                                                                                                                                   
	    if (pointsUBounds.size() == 0) {		 
		 pointsUBounds = trackletPoint;
	    }
	    if (pointsLBounds.size() == 0) {
		 pointsLBounds = trackletPoint;
	    }
	    extendBounds(pointsUBounds, trackletPoint, true);
	    extendBounds(pointsLBounds, trackletPoint, false);
	    
        }

	// C linkTracklets calculates the width of the FULL set of tracklets; 
	// if our parent gave us that set of values use them. otherwise calculate the
	// width of the tracklets in this actual image.
	std::vector<double> widthsToSend = perAxisWidths;
	if (perAxisWidths.size() == 0)
	{
	     widthsToSend.resize(4);
	     for (unsigned int i = 0; i < 4; i++) {		  
		  widthsToSend[i] = (pointsUBounds[i] - pointsLBounds[i]) / 2.0;
	     }
	}
        // build root TrackletTreeNode

        /* create the root of the tree (and the rest of the tree
         * recursively), save it to private var. */
        unsigned int idCounter = 0;
        myRoot = new TrackletTreeNode(parameterizedTracklets, 
                                      positionalErrorRa, 
                                      positionalErrorDec,
                                      maxLeafSize, 
                                      0,
				      widthsToSend,
                                      idCounter,
				      false, 
				      true);
        // don't set hasData until now, when the tree is actually built.
        mySize = idCounter;

        // *read* UBounds, UBounds from TrackletTreeNode (They have
        // *been updated* and extended to account for tracklet
        // *position/velocity error) and set our own
        myUBounds = *(myRoot->getUBounds());
        myLBounds = *(myRoot->getLBounds());
  
        hasData = true;
    }
}

     }}// close lsst::Mops
