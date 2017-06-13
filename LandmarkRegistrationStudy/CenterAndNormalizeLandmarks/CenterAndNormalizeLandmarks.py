import os
import unittest
import vtk, qt, ctk, slicer, numpy
from slicer.ScriptedLoadableModule import *
import logging

#
# CenterAndNormalizeLandmarks
#

class CenterAndNormalizeLandmarks(ScriptedLoadableModule):

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Center and normalize markups points"
    self.parent.categories = ["Scoliosis"]
    self.parent.dependencies = []
    self.parent.contributors = ["Ben Church"]
    self.parent.helpText = """

    """
    self.parent.acknowledgementText = """

""" # replace with organization, grant and thanks.

#
# CenterAndNormalizeLandmarksWidget
#

class CenterAndNormalizeLandmarksWidget(ScriptedLoadableModuleWidget):

  def setup(self):
    
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Interface Area
    #
    InterfaceCollapsibleButton = ctk.ctkCollapsibleButton()
    InterfaceCollapsibleButton.text = "Interface"
    self.layout.addWidget(InterfaceCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(InterfaceCollapsibleButton)

    #
    # Center Button
    #
    self.CenterButton = qt.QPushButton("Center all vtkMRMLMarkupsFiducialNodes")
    self.CenterButton.toolTip = "Subtracts average coordinate from each dimension, in each landmark set, centering the set."
    self.CenterButton.enabled = True
    parametersFormLayout.addRow(self.CenterButton)
    
    #
    # Normalize Button
    #
    self.NormalizeButton = qt.QPushButton("Normalize all vtkMRMLMarkupsFiducialNodes")
    self.NormalizeButton.toolTip = "Divides each point's coordinate by the maximum value in that dimension of that set, normalizing over [0,1]"
    self.NormalizeButton.enabled = True
    parametersFormLayout.addRow(self.NormalizeButton)
    
    #
    # Save Button
    #
    self.saveNewPointSets = qt.QPushButton("Save centered/normalized point sets")
    self.saveNewPointSets.toolTip = "Saves transverse process locations to .csv file containing their names and positions."
    self.saveNewPointSets.enabled = True
    parametersFormLayout.addRow(self.saveNewPointSets)
    
    # reload button
    self.reloadButton = qt.QPushButton('Restart module')
    self.layout.addWidget(self.reloadButton)
    
    # connections
    self.CenterButton.connect('clicked(bool)', self.onCenterButton)
    self.NormalizeButton.connect('clicked(bool)', self.onNormalizeButton)
    self.reloadButton.connect('clicked(bool)',self.onReloadButton)
    self.saveNewPointSets.connect('clicked(bool)', self.onSaveButton)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    #self.onSelect()

  def cleanup(self):
    pass

  #def onSelect(self):
   # self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onCenterButton(self):
    logic = CenterPointsLogic()
    logic.CenterPointSets()
    
  def onNormalizeButton(self):
    logic = NormalizeDataLogic()
    logic.NormalizeData()
    
  def onReloadButton(self):
    self.MarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for PointSet in range(self.MarkupNodes.__len__()):
      if (self.MarkupNodes.__getitem__(PointSet).GetNthFiducialLabel(0)[-1] == "-Z"):
        slicer.mrmlScene.RemoveNode(self.MarkupNodes.__getitem__(PointSet))
      if (self.MarkupNodes.__getitem__(PointSet).GetNthFiducialLabel(0)[-1] == "-N"):
        slicer.mrmlScene.RemoveNode(self.MarkupNodes.__getitem__(PointSet))
    slicer.util.reloadScriptedModule('CenterAndNormalizeLandmarks')
    
  # Assumes CalculateAngles of DegradeTransverseProcesses has been done
  def onSaveButton(self):
    import csv
    CenteredDataOutput = qt.QFileDialog.getSaveFileName(0, "Zeroed point sets", "", "CSV File (*.csv)")
    NormalizedDataOutput = qt.QFileDialog.getSaveFileName(0, "Normalized point sets", "", "CSV File (*.csv)")
    self.MarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    self.MarkupNodes = sorted(self.MarkupNodes, key = lambda DataSet: DataSet.GetName())
    AnglesTable = slicer.util.getModuleGui('DegradeTransverseProcesses').children()[4].children()[1]  
    #print (AnglesTable.item(0,1)[0])
    MaxAngles = []
    SupCritVert = []
    InfCritVert = []
    for InputSet in range(AnglesTable.rowCount):
      MaxAngles.append(float(AnglesTable.item(InputSet,1).data(0)))
      InfCritVert.append(AnglesTable.item(InputSet,2).data(0))
      SupCritVert.append(AnglesTable.item(InputSet,3).data(0))
      
            
    with open(CenteredDataOutput, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      CompleteSetNum = 0
      for MarkupsNode in range(self.MarkupNodes.__len__()):
        CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        if(CurrentLandmarkSet.GetName()[-1] == "Z"):
          # If the landmark set is a zeroed one
          writer.writerow(['Max angle:', str(MaxAngles[CompleteSetNum]), CurrentLandmarkSet.GetName(), '']) 
          CompleteSetNum += 1
          for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      writer.writerow(['EOF', '', '', ''])
      
    with open(NormalizedDataOutput, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      CompleteSetNum = 0
      for MarkupsNode in range(self.MarkupNodes.__len__()):
        CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        if(CurrentLandmarkSet.GetName()[-1] == "N"):
          # If the landmark set is a normalized one
          writer.writerow(['Max angle:', str(MaxAngles[CompleteSetNum]), InfCritVert[CompleteSetNum], SupCritVert[CompleteSetNum], CurrentLandmarkSet.GetName()]) 
          CompleteSetNum += 1
          for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      writer.writerow(['EOF', '', '', ''])
 
#
# CenterPointsLogic
#

class CenterPointsLogic(ScriptedLoadableModuleLogic):

  def __init__(self): 
    self.LandmarkPointSets = []      # Will contain tuples: (TrXFiducial point labels, TrXFiducial### point sets)
    self.AverageCoordinates = []      # Will contain one 3-tuple for each landmark set
    
  def CenterPointSets(self): 
    self.InputData = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    
    for InputSet in range(self.InputData.__len__()):
      CurrentLandmarkSet = self.InputData.__getitem__(InputSet)
      self.LandmarkPointSets.append([])
      print " "   #empty line
      print "Landmark set #" + str(InputSet)
      for InputPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
        self.LandmarkPointSets[InputSet].append((CurrentLandmarkSet.GetNthFiducialLabel(InputPoint),CurrentLandmarkSet.GetMarkupPointVector(InputPoint,0)))
        print self.LandmarkPointSets[InputSet][InputPoint] 
    
    self.FindAverageCoordinates()
    
    for i, CurrentLandmarkSet in enumerate(self.LandmarkPointSets):
      for LandmarkPoint in CurrentLandmarkSet:
        for dim in range(3):
          LandmarkPoint[1][dim] -= self.AverageCoordinates[i][dim]
    
    # Create new FiducialMarkupNodes for Slicer scene
    for InputSet in range(self.LandmarkPointSets.__len__()):
      NewMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewMarkupsNode.SetName(self.InputData.__getitem__(InputSet).GetName() + "-Z")
      CurrentLandmarkSet = self.LandmarkPointSets[InputSet]
      for InputPoint in range(CurrentLandmarkSet.__len__()):
        CurrentLandmarkPoint = CurrentLandmarkSet[InputPoint][1]
        NewMarkupsNode.AddFiducial(CurrentLandmarkPoint[0], CurrentLandmarkPoint[1], CurrentLandmarkPoint[2])
        NewPointLabel = self.LandmarkPointSets[InputSet][InputPoint][0]
        NewMarkupsNode.SetNthFiducialLabel(InputPoint, NewPointLabel)
      slicer.mrmlScene.AddNode(NewMarkupsNode)
           
    return True
    
  def FindAverageCoordinates(self):
    for CurrentLandmarkSet in self.LandmarkPointSets:
      CurrentRLAverage = 0
      CurrentAPAverage = 0
      CurrentSIAverage = 0
      for LandmarkPoint in CurrentLandmarkSet:
        CurrentRLAverage += LandmarkPoint[1][0]
        CurrentAPAverage += LandmarkPoint[1][1]
        CurrentSIAverage += LandmarkPoint[1][2]
      CurrentRLAverage = CurrentRLAverage / CurrentLandmarkSet.__len__()
      CurrentAPAverage = CurrentAPAverage / CurrentLandmarkSet.__len__()
      CurrentSIAverage = CurrentSIAverage / CurrentLandmarkSet.__len__()
      self.AverageCoordinates.append((CurrentRLAverage, CurrentAPAverage, CurrentSIAverage))

class NormalizeDataLogic(ScriptedLoadableModuleLogic):
  
  def __init__(self):
    self.LandmarkPointSets = []      # Will contain tuples: (TrXFiducial point labels, TrXFiducial### point sets)
    self.MaxAbsValCoordinates = []      # Will contain one 3-tuple for each maximum absolute value in each dimension
    
  def NormalizeData(self):
    self.InputData = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    
    for InputSet in range(self.InputData.__len__()):
      CurrentLandmarkSet = self.InputData.__getitem__(InputSet)
      self.LandmarkPointSets.append([])
      print " "   #empty line
      print "Landmark set #" + str(InputSet)
      for InputPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
        self.LandmarkPointSets[InputSet].append((CurrentLandmarkSet.GetNthFiducialLabel(InputPoint),CurrentLandmarkSet.GetMarkupPointVector(InputPoint,0)))
        print self.LandmarkPointSets[InputSet][InputPoint]
  
    self.FindMaxAbsVals()
    
    for i, CurrentLandmarkSet in enumerate(self.LandmarkPointSets):
      for LandmarkPoint in CurrentLandmarkSet:
        for dim in range(3):
          LandmarkPoint[1][dim] = LandmarkPoint[1][dim] / max(self.MaxAbsValCoordinates[i])
          
    # Create new FiducialMarkupNodes for Slicer scene
    for InputSet in range(self.LandmarkPointSets.__len__()):
      NewMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewMarkupsNode.SetName(self.InputData.__getitem__(InputSet).GetName() + "-N")
      CurrentLandmarkSet = self.LandmarkPointSets[InputSet]
      for InputPoint in range(CurrentLandmarkSet.__len__()):
        CurrentLandmarkPoint = CurrentLandmarkSet[InputPoint][1]
        NewMarkupsNode.AddFiducial(CurrentLandmarkPoint[0], CurrentLandmarkPoint[1], CurrentLandmarkPoint[2])
        NewPointLabel = self.LandmarkPointSets[InputSet][InputPoint][0]
        NewMarkupsNode.SetNthFiducialLabel(InputPoint, NewPointLabel)
      slicer.mrmlScene.AddNode(NewMarkupsNode)
          
  def FindMaxAbsVals(self):
    for CurrentLandmarkSet in self.LandmarkPointSets:
      CurrentSetRLMax = 0
      CurrentSetAPMax = 0
      CurrentSetSIMax = 0
      for CurrentLandmarkPoint in CurrentLandmarkSet:
        if abs(CurrentLandmarkPoint[1][0]) > CurrentSetRLMax:
          CurrentSetRLMax = abs(CurrentLandmarkPoint[1][0])
        if abs(CurrentLandmarkPoint[1][1]) > CurrentSetAPMax:
          CurrentSetAPMax = abs(CurrentLandmarkPoint[1][1])
        if abs(CurrentLandmarkPoint[1][2]) > CurrentSetSIMax:
          CurrentSetSIMax = abs(CurrentLandmarkPoint[1][2])
      self.MaxAbsValCoordinates.append((CurrentSetRLMax, CurrentSetAPMax, CurrentSetSIMax))
    
class CenterAndNormalizeLandmarksTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_CenterAndNormalizeLandmarks1()

  def test_CenterAndNormalizeLandmarks1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = CenterAndNormalizeLandmarksLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
