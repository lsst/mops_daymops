// jmyers

// this won't compile but it's added here for posterity.  

// this is code Matt wrote to determine most optimal leaf node size on
// your system (allowing leaf node sizes which are powers of two),
// which pretty cool.  For now we just trust him that it's 16 or 32...


/*********************************************************************************
 * 
 * ADAPTIVE LNS ALGORITHM 
 * 
 *********************************************************************************/
int adaptiveAlgorithm(const std::vector<Detection> &allDetections,
		      std::vector<Tracklet> &allTracklets,
		      linkTrackletsConfig searchConfig,
		      std::map<ImageTime, KDTree::KDTree <unsigned int> > &trackletTimeToTreeMap,
		      std::map<unsigned int, std::map<unsigned int, unsigned int> > &parentTable,
		      std::vector<std::vector<int> > *assignment, workItemFile *workItemCache,
		      int LEAF_SIZE, std::string dir)
{
  std::cerr << "Entering adaptive algorithm" << std::endl;

  //this was a #define but i'm trying to minimize globals that aren't needed by both
  //master and workers
  int MAX_ADAPTIVE_LNS = 64;

  //create a map of LNS to work item vectors
  std::map<unsigned int, std::vector<int>*> workItemMap;
  int currentLNS = LEAF_SIZE;
  workItemMap[currentLNS] = &assignment->at(0);

  //keep track of timeUnits at each LNS
  std::map<int, long unsigned int> LNSTimeMap;

  //this multimap stores endpoint pairs that have already been visited
  std::vector<std::pair<treeIdNodeIdPair, treeIdNodeIdPair> > endpointsVisited;
  endpointsVisited.clear();

  //endpoint nodes
  TreeNodeAndTime *firstEndpoint;
  TreeNodeAndTime *secondEndpoint;
  std::vector<TreeNodeAndTime*> newSupportNodes;

  //variables for work item creation
  double RAVelocity, DecVelocity, RAAcceleration, DecAcceleration;
  double RAPosition0, DecPosition0;
  double time0;
  int numCompatible = 0;

  //perform time aggregation on all LNSes
  while(currentLNS <= MAX_ADAPTIVE_LNS){

    std::cerr << "Adaptive algorithm evaluating LNS of " << currentLNS << std::endl;
    if(workItemMap.find(currentLNS) != workItemMap.end()){
      
      unsigned int count = 0;
      //gather work item time for each work item of this LNS
      for(unsigned int i=0; i < workItemMap[currentLNS]->size(); ++i){
	
	//load this work item
	workItemFile adaptiveItem = readWorkItemFileFromDisk(workItemMap[currentLNS]->at(i), workItemCache, dir);

	//nodeId and treeId of the first endpoint
	//unsigned int treeId = workItemMap[currentLNS]->at(i).endpoints.at(0).first;
	//int nodeId = workItemMap[currentLNS]->at(i).endpoints.at(0).second;
	unsigned int treeId = adaptiveItem.endpoints.at(0).first;
	int nodeId = adaptiveItem.endpoints.at(0).second;
	
	//Get the tree and ImageTime for this image (tree) ID
	KDTree::KDTree<unsigned int> *myTree = findTreeById(treeId, trackletTimeToTreeMap);
	ImageTime it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	
	//find this node in its tree
	KDTree::KDTreeNode<unsigned int> *treeNode = NULL;
	if( myTree != NULL ){
	  treeNode = getNodeByIDAndTime(myTree, nodeId);
	}
	
	//ensure proper node was found
	if( treeNode != NULL ){
	  firstEndpoint = new TreeNodeAndTime(treeNode, it);
	}	
	
	//nodeId and treeId of second endpoint
	//treeId = workItemMap[currentLNS]->at(i).endpoints.at(1).first;
	//nodeId = workItemMap[currentLNS]->at(i).endpoints.at(1).second;
	treeId = adaptiveItem.endpoints.at(1).first;
	nodeId = adaptiveItem.endpoints.at(1).second;
	
	//Get the tree and ImageTime for this image (tree) ID
	myTree = findTreeById(treeId, trackletTimeToTreeMap);
	it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	
	//find this node in its tree
	treeNode = NULL;
	if( myTree != NULL ){
	  treeNode = getNodeByIDAndTime(myTree, nodeId);
	}
	
	//ensure proper node was found
	if( treeNode != NULL ){
	  secondEndpoint = new TreeNodeAndTime(treeNode, it);
	}	
	
	//gather up support points
	//for(unsigned int j=2; j < workItemMap[currentLNS]->at(i).endpoints.size(); ++j){
	for(unsigned int j=2; j < adaptiveItem.endpoints.size(); ++j){
	  
	  //treeId = workItemMap[currentLNS]->at(i).endpoints.at(j).first;
	  //nodeId = workItemMap[currentLNS]->at(i).endpoints.at(j).second;
	  treeId = adaptiveItem.endpoints.at(j).first;
	  nodeId = adaptiveItem.endpoints.at(j).second;
	  
	  //Get the tree and ImageTime for this image (tree) ID
	  myTree = findTreeById(treeId, trackletTimeToTreeMap);
	  it = findImageTimeForTree(treeId, trackletTimeToTreeMap);
	  
	  //find this node in its tree
	  treeNode = NULL;
	  if( myTree != NULL ){
	    treeNode = getNodeByIDAndTime(myTree, nodeId);
	  }
	  
	  //ensure proper node was found
	  if( treeNode != NULL ){
	    TreeNodeAndTime *TAT = new TreeNodeAndTime(treeNode, it);
	    newSupportNodes.push_back(TAT);
	  }	
	}
	
	std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator firstEndpointIter;
	std::vector<KDTree::PointAndValue<unsigned int> >::const_iterator secondEndpointIter;
	const std::vector<KDTree::PointAndValue<unsigned int> > *firstEndpointData = firstEndpoint->myTree->getMyData();
	const std::vector<KDTree::PointAndValue<unsigned int> > *secondEndpointData = secondEndpoint->myTree->getMyData();
	
	for (firstEndpointIter = firstEndpointData->begin();
	     firstEndpointIter != firstEndpointData->end();
	     firstEndpointIter++) {
	  
	  for (secondEndpointIter = secondEndpointData->begin();
	       secondEndpointIter != secondEndpointData->end();
	       secondEndpointIter++) {
	    
	    getBestFitVelocityAndAccelerationForTracklets(allDetections,
							  allTracklets, 
							  firstEndpointIter->getValue(),
							  secondEndpointIter->getValue(),
							  RAVelocity,RAAcceleration,RAPosition0,
							  DecVelocity,DecAcceleration,DecPosition0,
							  time0);
	    
	    if (endpointTrackletsAreCompatible(allDetections, 
					       allTracklets,
					       firstEndpointIter->getValue(),
					       secondEndpointIter->getValue(),
					       RAVelocity,RAAcceleration,RAPosition0,
					       DecVelocity,DecAcceleration,DecPosition0,
					       time0,searchConfig)) {
	      ++numCompatible;
	      ++compatibleEndpointsFound;
	    }
	  }
	}
	
	//calculate the amount of time this work item requires and update LNSTimeMap
	//workItemMap[currentLNS]->at(i).timeUnits = ((2 + newSupportNodes.size()) * numCompatible);
	adaptiveItem.timeUnits = ((2 + newSupportNodes.size()) * numCompatible);
	
	/*if(LNSTimeMap.find(currentLNS) == LNSTimeMap.end()){
	  LNSTimeMap[currentLNS] = workItemMap[currentLNS]->at(i).timeUnits;
	}
	else{
	  LNSTimeMap[currentLNS] += workItemMap[currentLNS]->at(i).timeUnits;
	  }*/

	if(LNSTimeMap.find(currentLNS) == LNSTimeMap.end()){
	  LNSTimeMap[currentLNS] = adaptiveItem.timeUnits;
	}
	else{
	  LNSTimeMap[currentLNS] += adaptiveItem.timeUnits;
	}
	
	/******************************************************
	 * create work item for this work item's parent points
	 ******************************************************/
	//these endpoints are those for the work item of this work item's parent
	std::vector<std::pair<unsigned int, unsigned int> > endpoints;
	endpoints.clear();
	
	//first item is first endpoint
	treeId = firstEndpoint->myTime.getImageId();
	nodeId = firstEndpoint->myTree->getId();
	unsigned int parentNodeId = parentTable[treeId][nodeId];
	endpoints.push_back(std::make_pair(treeId, parentNodeId));
	
	//second item is second endpoint
	treeId = secondEndpoint->myTime.getImageId();
	nodeId = secondEndpoint->myTree->getId();
	parentNodeId = parentTable[treeId][nodeId];
	endpoints.push_back(std::make_pair(treeId, parentNodeId));
	
	//add all endpoints
	for(unsigned int m=0; m < newSupportNodes.size(); ++m){
	  treeId = newSupportNodes.at(m)->myTime.getImageId();
	  nodeId = newSupportNodes.at(m)->myTree->getId();
	  parentNodeId = parentTable[treeId][nodeId];
	  endpoints.push_back(std::make_pair(treeId, parentNodeId));
	}
	
	//create workItemFile struct to go on list
	workItemFile myItem;
	myItem.endpoints = endpoints;
	myItem.numNodes = 2 + newSupportNodes.size();
	myItem.numCompatibleEndpoints = 0;                 //numCompatible  -- we don't have this value yet
	myItem.timeUnits = myItem.numNodes * myItem.numCompatibleEndpoints;  //we don't have this value yet
	
	//clear this vector -- it gets big and is workItem-specific
	newSupportNodes.clear();
	
	//add work item to list at corresponding LNS slot
	int parentLNS = currentLNS * 2;
	if(workItemMap.find(parentLNS) == workItemMap.end()){
	  //workItemMap[parentLNS] = new std::vector<workItemFile>();
	  workItemMap[parentLNS] = new std::vector<int>();
	  workItemMap[parentLNS]->clear();
	}
	
	//std::cerr << "checking to see if this endpoint pair exists" << std::endl;
	
	//check to see if this work item already exists, don't add duplicate
	bool exists = false;
	for(unsigned int k=0; k < endpointsVisited.size(); ++k){
	  if((endpointsVisited.at(k).first.treeId == 
	      myItem.endpoints.at(0).first) &&
	     (endpointsVisited.at(k).first.nodeId ==
	      myItem.endpoints.at(0).second)){
	    
	    //only check second endpoint if first was found
	    if((endpointsVisited.at(k).second.treeId ==
		myItem.endpoints.at(1).first) &&
	       (endpointsVisited.at(k).second.nodeId ==
		myItem.endpoints.at(1).second)){
	      exists = true;
	      break;
	    }
	  }
	}
	
	//add this work item to the list if it hasn't already been added
	if(!exists){	  
	  //std::cerr << "doesn't exist, adding to endpointsVisited" << std::endl;
	  
	  //update the list of workitem endpoints visited
	  treeIdNodeIdPair p1;
	  p1.treeId = myItem.endpoints.at(0).first;
	  p1.nodeId = myItem.endpoints.at(0).second;

	  treeIdNodeIdPair p2;
	  p2.treeId = myItem.endpoints.at(1).first;
	  p2.nodeId = myItem.endpoints.at(1).second;

	  endpointsVisited.push_back(std::make_pair(p1, p2));
	  
	  //add this new work item  to thelist associated with the parents
	  //of the  current  work item
	  //workItemMap[parentLNS]->push_back(myItem);
	  workItemMap[parentLNS]->push_back(parentTable[treeId][nodeId]);
	}
	//std::cerr << "work item loop finished iteration " << i << std::endl;
	count = i;
      } /* end work item loop */ 
      std::cerr << "LNS " << currentLNS << " had " << count << " work items!?!" << std::endl;
    } /* LNS exists check */
    
    //this list is unique to each LNS
    endpointsVisited.clear();

    //update the LNS to the next higher power of 2
    currentLNS *= 2;
  } /* end adapative algorithm loop */

  std::cerr << "Finished adaptive loop" << std::endl;

  //find the LNS with lowest total work units and return it
  long unsigned int minWorkUnits = ULONG_MAX;
  int optimalLNS = LEAF_SIZE;
  bool first = true;

  std::map<int, long unsigned int>::const_iterator LNSiter;
  
  for(LNSiter = LNSTimeMap.begin(); LNSiter != LNSTimeMap.begin(); LNSiter++){
    if(((*LNSiter).second < minWorkUnits) || first){
      minWorkUnits = (*LNSiter).second;
      optimalLNS = (*LNSiter).first;
      std::cerr << "updating minWorkUnits to " << minWorkUnits << std::endl;
      std::cerr << "optimalLNS is now " << optimalLNS << std::endl;

      if(first){
	first = false;
      }
    }
  }

  //print adaptive algorithm results
  std::cerr << "LNS Total Work times" << std::endl;
  std::cerr << "----------------------------------------" << std::endl;
  std::map<int, unsigned long int>::iterator it;
  for (it = LNSTimeMap.begin() ; it != LNSTimeMap.end(); it++ )
    std::cerr << (*it).first << " => " << (*it).second << std::endl;

  return optimalLNS;
} 
