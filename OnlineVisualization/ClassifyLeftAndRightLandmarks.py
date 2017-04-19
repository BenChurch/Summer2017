import numpy as np



# Returns True if k-means is done. K-means terminates either
# because it has run a maximum number of iterations OR the centroids
# stop changing.
def ShouldStopKMeans(oldCentroids, Centroids, iterations):
  MaxIterations = 500
  StoppingDelta = 0.05
  if iterations > MaxIterations: return True
  if iterations == 0:
    return False
  for C in range(len(oldCentroids)):
    if np.linalg.norm(oldCentroids[C] - Centroids[C]) < StoppingDelta:
      return True
  return False

def GetKMeansLabels(DataSet, Centroids, Labels):
  for i, DataPoint in enumerate(DataSet):
    PointVector = np.array(DataPoint)
    minDist = np.linalg.norm(PointVector - Centroids[0])
    Labels[i] = 0
    for j, CentroidVector in enumerate(Centroids):
      PointCentroidDistance = np.linalg.norm(PointVector - CentroidVector)
      if PointCentroidDistance < minDist:
        minDist = PointCentroidDistance
        Labels[i] = j
  return Labels

def GetCentroids(DataSet, Labels, k):
  Centroids = []
  for C in range(k):    # For each centroid
    Centroid = np.random.uniform() * np.ones(len(DataSet[0]))    # Each centroid with as many dimensions as the data
    
    for i, DataPoint in enumerate(DataSet): # Take each data point contributing to the centroid into consideration
      if Labels[i] == C:                    # if it belongs to the current centroid
        for dim in range(len(Centroid)):
          Centroid[dim] += DataPoint[dim]
          
    for dim in range(len(Centroid)):
      Centroid[dim] = Centroid[dim] / np.count_nonzero(Labels == C)
    Centroids.append(Centroid)
  return Centroids

def KMeans(DataSet, k):   # Expects numpy array for DataSet
  DataSetLabels = -1 * np.ones(len(DataSet))
  # Initialize centroids
  numFeatures = len(DataSet[0])
  Centroids = []
  for Cent in range(k):
    Centroids.append(np.random.uniform() * np.ones(numFeatures))
  
# Initialize book keeping variables
  iterations = 0
  oldCentroids = Centroids
  
  # Run the k-means algorithm
  while not ShouldStopKMeans(oldCentroids, Centroids, iterations):
    oldCentroids = Centroids
    iterations += 1
    
    DataSetLabels = GetKMeansLabels(DataSet, Centroids, DataSetLabels)
    print DataSetLabels
    Centroids = GetCentroids(DataSet, DataSetLabels, k)
  return (DataSetLabels, Centroids)

MarkupsNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
MarkupsNode = MarkupsNodes[0]

# Start by finding top and bottom points of spine for verticle normalization
SpineLeft = 10000
SpineRight = -10000
SpineFront = -10000
SpineBack = 10000
SpineTop = -10000
SpineBottom = 10000

for PointIndex in range(MarkupsNode.GetNumberOfFiducials()):
  Point  = MarkupsNode.GetMarkupPointVector(PointIndex, 0)
  if Point[0] < SpineLeft:
    SpineLeft = Point[0]
  if Point[0] > SpineRight:
    SpineRight = Point[0]
  if Point[1] < SpineBack:
    SpineBack = Point[1]
  if Point[1] > SpineFront:
    SpineFront = Point[1]
  if Point[2] > SpineTop:
    SpineTop = Point[2]
  if Point[2] < SpineBottom:
    SpineBottom = Point[2]
  

SpineHeight = SpineTop - SpineBottom
SpineWidth = SpineRight - SpineLeft
SpineDepth = SpineFront - SpineBack

UnclassifiedPoints = []
# Normalize S-I dimension to R-L scale
for PointIndex in range(MarkupsNode.GetNumberOfFiducials()):
  Point = MarkupsNode.GetMarkupPointVector(PointIndex, 0)
  Point[1] = (Point[1]) * (SpineWidth / (SpineDepth * 3.0))
  Point[2] = (Point[2]) * (SpineWidth / (SpineHeight * 2.0))
  UnclassifiedPoints.append(np.array(Point))

(KmLabels, KmCentroids) = KMeans(UnclassifiedPoints, 2)

LeftMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
LeftMarkupsNode.SetName(MarkupsNode.GetName() + 'Left')

RightMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
RightMarkupsNode.SetName(MarkupsNode.GetName() + 'Right')

# If KmLabel == 0 indicates a left-side point
if KmCentroids[0][0] < KmCentroids[1][0]:
  for i, UnclassifiedPoint in enumerate(UnclassifiedPoints):
    if KmLabels[i] == 0:
      LeftMarkupsNode.AddFiducialFromArray(MarkupsNode.GetMarkupPointVector(i,0))
      LeftMarkupsNode.SetNthFiducialLabel(LeftMarkupsNode.GetNumberOfFiducials()-1, MarkupsNode.GetNthFiducialLabel(i) + '_Left')
    else:
      RightMarkupsNode.AddFiducialFromArray(MarkupsNode.GetMarkupPointVector(i,0))
      RightMarkupsNode.SetNthFiducialLabel(RightMarkupsNode.GetNumberOfFiducials()-1, MarkupsNode.GetNthFiducialLabel(i) + '_Right')
else: # If KmLabel == 0 indicates a right-side point
  for i, UnclassifiedPoint in enumerate(UnclassifiedPoints):
    if KmLabels[i] == 0:
      RightMarkupsNode.AddFiducialFromArray(MarkupsNode.GetMarkupPointVector(i,0))
      RightMarkupsNode.SetNthFiducialLabel(RightMarkupsNode.GetNumberOfFiducials()-1, MarkupsNode.GetNthFiducialLabel(i) + '_Right')
    else:
      LeftMarkupsNode.AddFiducialFromArray(MarkupsNode.GetMarkupPointVector(i,0))
      LeftMarkupsNode.SetNthFiducialLabel(LeftMarkupsNode.GetNumberOfFiducials()-1, MarkupsNode.GetNthFiducialLabel(i) + '_Left')

slicer.mrmlScene.AddNode(LeftMarkupsNode)
slicer.mrmlScene.AddNode(RightMarkupsNode)
