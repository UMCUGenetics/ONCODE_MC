import numpy as np

def computeDBIndex(data, centers, labels):

    numClusters = centers.shape[0]
    DB = 0
    for clusterInd in range(0, numClusters):
        maxDij = 0
        for clusterInd2 in range(0, numClusters): #compute the max for each cluster to all other clusters
            #Skip if the cluster is itself
            if clusterInd == clusterInd2:
                continue
             
            #get all points for cluster i
            clusterPoints = data[labels == clusterInd,:]
            center = centers[clusterInd,:]
            
            #Compute the euclidean distance from each point to the center of the cluster
            totalDistanceToCenter = []
            for pointInd in range(0, clusterPoints.shape[0]):
                point = clusterPoints[pointInd, :]
                #Compute the euclidean distance for this point to the center of its cluster
                totalDistanceToCenter.append(np.linalg.norm(point-center))
            #sum these distances and divide it by the total # of points in this cluster
            #This is the average distance from all points to the center. 
            
            di = sum(totalDistanceToCenter) / clusterPoints.shape[0]

            #get all points for cluster j, again compute the average distance of the points to the center. 
            clusterPoints = data[labels == clusterInd2,:]
            center = centers[clusterInd2,:]
            totalDistanceToCenter = []
            for pointInd in range(0, clusterPoints.shape[0]):
                point = clusterPoints[pointInd, :]

                totalDistanceToCenter.append(np.linalg.norm(point-center))

            dj = sum(totalDistanceToCenter) / clusterPoints.shape[0]

            #Compute the euclidean distance between the centers of these two clusters
            dij = np.linalg.norm(centers[clusterInd,:]-centers[clusterInd2,:])
            #compute the ratio within-cluster/between-cluster
            Dij = (di + dj) / dij
            #obtain the maximum ratio. For which two clusters is the within-to-between ratio the worst? 
            if Dij > maxDij:
                maxDij = Dij
        DB += maxDij #sum the maxima obtained for each cluster.
    dbIndex = (1/float(numClusters))*DB
    return dbIndex
