import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np

#
# DegradeTransverseProcesses
#

class DegradeTransverseProcesses(ScriptedLoadableModule):

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Degrade Transverse Processes"
    self.parent.categories = ["Scoliosis"]
    self.parent.dependencies = []
    self.parent.contributors = ["Ben Church"]
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# DegradeTransverseProcessesWidget
#

class DegradeTransverseProcessesWidget(ScriptedLoadableModuleWidget):

  def setup(self):
    
    self.DegradedNodeTag = "~"
    
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...
    
    #
    # Config Area
    #
    ConfigArea = ctk.ctkCollapsibleButton()
    ConfigArea.text = "Config"
    self.layout.addWidget(ConfigArea)

    # Layout within the dummy collapsible button
    ConfigLayout = qt.QGridLayout(ConfigArea)
    ConfigLayout.setHorizontalSpacing(12)
    #parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    
    #
    # input volume selector
    #
    # Don't use this - iterate through entire batch of TrXFiducials loaded
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    #parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # noise standard deviation
    #
    self.StdDevSliderWidget = ctk.ctkSliderWidget()
    self.StdDevSliderWidget.singleStep = 0.01
    self.StdDevSliderWidget.minimum = 0
    self.StdDevSliderWidget.maximum = 10
    self.StdDevSliderWidget.value = 1
    self.StdDevSliderWidget.setToolTip("Set standard deviation (mm) of noise to introduce to all points.")
    ConfigLayout.addWidget(qt.QLabel("Noise StdDev (mm^2)"), 2, 0, 1, 1)
    ConfigLayout.addWidget(self.StdDevSliderWidget, 2, 1, 1, 3)
    #parametersFormLayout.addRow("Noise standard deviation", self.StdDevSliderWidget)

    #
    # deletion fraction
    #
    self.DeletionSliderWidget = ctk.ctkSliderWidget()
    self.DeletionSliderWidget.singleStep = 0.01
    self.DeletionSliderWidget.minimum = 0
    self.DeletionSliderWidget.maximum = 1
    self.DeletionSliderWidget.value = 0.5
    self.DeletionSliderWidget.setToolTip("Set fraction of points from original sets to delete.")
    ConfigLayout.addWidget(qt.QLabel("Deletion fraction"), 3, 0, 1, 1)
    ConfigLayout.addWidget(self.DeletionSliderWidget, 3, 1, 1, 3)
    #parametersFormLayout.addRow("Deletion fraction", self.DeletionSliderWidget)

    #
    # displacement fraction
    #
    self.DisplacementSliderWidget = ctk.ctkSliderWidget()
    self.DisplacementSliderWidget.singleStep = 0.01
    self.DisplacementSliderWidget.minimum = 0
    self.DisplacementSliderWidget.maximum = 1
    self.DisplacementSliderWidget.value = 0.5
    self.DisplacementSliderWidget.setToolTip("Set fraction of points from original data set to misplace.")
    ConfigLayout.addWidget(qt.QLabel("Displacement fraction"), 4, 0, 1, 1)
    ConfigLayout.addWidget(self.DisplacementSliderWidget, 4, 1, 1, 3)
    #parametersFormLayout.addRow("Displacement fraction", self.DisplacementSliderWidget)
    
    #
    # Checkbox to enable protection of boundaries from deletion
    #
    self.DontDeleteBoundaryPointsBox = qt.QCheckBox()
    self.DontDeleteBoundaryPointsBox.toolTip = "If checked, four outer-most points of spine will not be deleted"
    self.DontDeleteBoundaryPointsBox.enabled = True
    ConfigLayout.addWidget(qt.QLabel("Keep boundaries"), 5, 0, 1, 1)
    ConfigLayout.addWidget(self.DontDeleteBoundaryPointsBox, 5, 1, 1, 1)
    #parametersFormLayout.addRow(self.DontDeleteBoundaryPointsBox)
    
    #
    # Degrade Button
    #
    self.DegradeButton = qt.QPushButton("Degrade point sets")
    self.DegradeButton.toolTip = "Apply noise, delete fraction of random points, and misplace fraction of random points."
    self.DegradeButton.enabled = True
    ConfigLayout.addWidget(self.DegradeButton, 6, 0, 1, 4)
    #parametersFormLayout.addRow(self.DegradeButton)
    
    # Output area
    anglePanel = ctk.ctkCollapsibleButton()
    anglePanel.text = "Angle measurement"
    self.layout.addWidget(anglePanel)
    
    outputVerticalLayout = qt.QGridLayout(anglePanel)

    TableHeaders = ["Landmark set", "Max angle", "InfCritVert", "SupCritVert"]
    self.angleTable = qt.QTableWidget((slicer.mrmlScene.GetNodesByClass('vtkMRMLMarkupsFiducialNode').GetNumberOfItems()), 4)
    self.angleTable.sortingEnabled = False
    self.angleTable.setEditTriggers(0)
    self.angleTable.setMinimumHeight(self.angleTable.verticalHeader().length() + 25)
    self.angleTable.horizontalHeader().setResizeMode(qt.QHeaderView.Stretch)
    outputVerticalLayout.addWidget(self.angleTable, 2, 0, 1, 5)
    self.angleTable.setHorizontalHeaderLabels(TableHeaders)
    
    self.CalculateAnglesButton = qt.QPushButton("Calculate angles")
    self.savePointsWithAnglesButton = qt.QPushButton("Save altered data + angles")
    outputVerticalLayout.addWidget(self.CalculateAnglesButton, 0, 0, 1, 3)
    outputVerticalLayout.addWidget(self.savePointsWithAnglesButton, 0, 3, 1, 3)
    
    self.ReloadModuleButton = qt.QPushButton("Reload module")
    self.layout.addWidget(self.ReloadModuleButton)
    
    # connections
    self.DegradeButton.connect('clicked(bool)', self.onDegradeButton)
    self.CalculateAnglesButton.connect('clicked(bool)', self.onCalculateAnglesButton)
    self.ReloadModuleButton.connect('clicked(bool)',self.onReloadModuleButton)
    self.savePointsWithAnglesButton.connect('clicked(bool)', self.onSaveButton)
    #self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    #self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    #self.onSelect()

  def cleanup(self):
    pass

  #def onSelect(self):
   # self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onDegradeButton(self):
    # Get parameters
    NoiseStdDev = self.StdDevSliderWidget.value
    DeletionFraction = self.DeletionSliderWidget.value
    DisplacementFraction = self.DisplacementSliderWidget.value
    PreserveBoundaries = self.DontDeleteBoundaryPointsBox.isChecked()
    
    # Get nodes to degrade (all of them)
    OriginalNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    
    # Degrade each original node
    logic = DegradeTransverseProcessesLogic()
    for Node in OriginalNodes:
      NewestDegradedNode = logic.DegradeInputData(Node, NoiseStdDev, DeletionFraction, DisplacementFraction, PreserveBoundaries)
      NewestDegradedNode.SetName(Node.GetName() + self.DegradedNodeTag)
      slicer.mrmlScene.AddNode(NewestDegradedNode)
    
    #imageThreshold = self.imageThresholdSliderWidget.value
    
  def onCalculateAnglesButton(self):
    AllMarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    logic = CalculateAnglesLogic()
    self.Angles = logic.CalculateAngles()
    (self.MaxAngles, self.CriticalVertebrae) = logic.FindMaxCoronalAngles()
    self.PopulateAnglesTable(self.MaxAngles, self.CriticalVertebrae, AllMarkupNodes)
    
  def onReloadButton(self):
    AllMarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for PointSet in range(AllMarkupNodes.__len__()):
      if (AllMarkupNodes.__getitem__(PointSet).GetNthFiducialLabel(0)[-1] == "~"):
        slicer.mrmlScene.RemoveNode(AllMarkupNodes.__getitem__(PointSet))
    slicer.util.reloadScriptedModule('DegradeTransverseProcesses')
    
  # Assumes CalculateAngles has been done
  def onSaveButton(self):
    import csv
    OriginalDataOutput = qt.QFileDialog.getSaveFileName(0, "Unmodified point data", "", "CSV File (*.csv)")
    ModifiedDataOutput = qt.QFileDialog.getSaveFileName(0, "Modified point data", "", "CSV File (*.csv)")
    self.MarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    
    with open(OriginalDataOutput, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      for MarkupsNode in range(self.MarkupNodes.__len__()):
        CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        if(CurrentLandmarkSet.GetName()[-1] != "~"):
          # If the landmark set is not a modified one
          #writer.writerow(['', '', '', ''])
          writer.writerow(['MaxAngle:', str(self.MaxAngles[MarkupsNode]), CurrentLandmarkSet.GetName(), ''])
          #writer.writerow(['', '', '', ''])
          for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      writer.writerow(['EOF', '', '', ''])
            
    with open(ModifiedDataOutput, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      for MarkupsNode in range(self.MarkupNodes.__len__()):
        CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        if(CurrentLandmarkSet.GetName()[-1] == "~"):
          # If the landmark set is a modified one
          #writer.writerow(['', '', '', ''])
          writer.writerow(['MaxAngle:', str(self.MaxAngles[MarkupsNode]), CurrentLandmarkSet.GetName(), '']) 
          #writer.writerow(['', '', '', ''])
          for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      writer.writerow(['EOF', '', '', ''])
    
    
  def PopulateAnglesTable(self, MaxAngles, CriticalVertebrae, MarkupPoints):
    SortedAngles  = sorted(zip(MarkupPoints, MaxAngles), key = lambda DataSet: DataSet[0].GetName())
    SortedMinVertebraeZip = sorted(zip(MarkupPoints, CriticalVertebrae[0]), key = lambda DataSet: DataSet[0].GetName())
    SortedMinVertebrae = zip(*SortedMinVertebraeZip)[1]
    #print SortedMinVertebrae
    SortedMaxVertebraeZip = sorted(zip(MarkupPoints, CriticalVertebrae[1]), key = lambda DataSet: DataSet[0].GetName())
    SortedMaxVertebrae = zip(*SortedMaxVertebraeZip)[1]
    #print SortedMaxVertebrae
    SortedMarkupPoints = zip(*SortedAngles)[0]
    SortedMaxAngles = zip(*SortedAngles)[1]
    for i, Angle in enumerate(SortedMaxAngles):
      CurrentLandmarkSet = SortedMarkupPoints.__getitem__(i)
      self.angleTable.setItem(i, 0, qt.QTableWidgetItem())
      self.angleTable.setItem(i, 1, qt.QTableWidgetItem())
      self.angleTable.setItem(i, 2, qt.QTableWidgetItem())
      self.angleTable.setItem(i, 3, qt.QTableWidgetItem())
      self.angleTable.item(i, 0).setText(CurrentLandmarkSet.GetName())
      self.angleTable.item(i, 1).setText(str(Angle))
      self.angleTable.item(i, 2).setText(SortedMinVertebrae[i])
      self.angleTable.item(i, 3).setText(SortedMaxVertebrae[i])

  def onReloadModuleButton(self):
    # Delete degraded markups nodes
    AllNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    AllNodeNames = [AllNodes[i].GetName() for i in range(len(AllNodes))]
    DegradedNodes = [Node for (i,Node) in enumerate(AllNodes) if (AllNodeNames[i]).__contains__(self.DegradedNodeTag)]
    for Node in DegradedNodes:
      slicer.mrmlScene.RemoveNode(Node)
      #if Node.__contains__(self.DegradedNodeTag): 
      #  slicer.mrmlScene.RemoveNode(AllNodes[i])
      
    # Reload the module
    slicer.util.reloadScriptedModule(slicer.moduleNames.DegradeTransverseProcesses)
    print slicer.moduleNames.DegradeTransverseProcesses, " reloaded"
    
#
# DegradeTransverseProcessesLogic
#

class DegradeTransverseProcessesLogic(ScriptedLoadableModuleLogic):

  def __init__(self):
    
    self.NodeLabelCoords = []      # Will contain tuples: (TrXFiducial point labels, TrXFiducial### point sets) of degraded points (has omission, noise, etc)
    self.DeletionLabelIndex = [[],[]]   # Will ensure we don't try to misplace deleted points or vice-versa
    self.SortedDeletionLandmarkLabels = [[],[]]
    self.MisplacedLandmarkLabels = [[],[]]
    self.RightLeftScale = 1000          # Length of vector in coronal plane used to misplace point to behind a rib
    self.AntPostScale = 20   # Length of vector in parasagital plane used to misplace point from behind rib onto rib
    self.MaxNoiseMM = 20    # 
    
  def DegradeInputData(self, Node, NoiseStdDev, DeletionFraction, DisplacementFraction, PreserveBoundaries=False):
    
    if(DeletionFraction + DisplacementFraction > 1):
      print "Error - DeletionFraction + DisplacementFraction > 1"
      print " Cannot displace points once deleted, or reduce DisplacementFraction by deleteing displaced points"
      print " Returning without degrading points"
      return
      
    self.InputData = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    self.NoiseStdDev = NoiseStdDev
    self.DeletionFraction = DeletionFraction
    self.DisplacementFraction = DisplacementFraction
    
    
    #for InputSet in range(self.InputData.__len__()):
      #CurrentLandmarkSet = self.InputData.__getitem__(InputSet)
      
    # If PreserveBoundaries, must make sure there are enough non-boundary points to delete to achieve DeletionFraction
    if PreserveBoundaries and ((Node.GetNumberOfFiducials()*DeletionFraction) - 4) <= 0:
      print "Error - cannot preserve boundaries for ", Node.GetName()
      print " Not enough non-boundary points to satisyfy DeletionFraction=", self.DeletionFraction
      print " Markups nodes left un-modified"
      return
      
    #self.NodeLabelCoords.append([])
    #print " "   #empty line
    #print "Landmark set #" + str(InputSet)
    
    #LabelsCoords = [(Node.GetNthFiducialLabel(i), Node.GetMarkupPointVector(i,0)) for i in range(Node.GetNumberOfFiducials())]
    
    for InputPoint in range(Node.GetNumberOfFiducials()):
      self.NodeLabelCoords.append((Node.GetNthFiducialLabel(InputPoint),Node.GetMarkupPointVector(InputPoint,0)))
      #print self.NodeLabelCoords[InputSet][InputPoint]
    
    # Initialize these lists here so their functions can be called in either order
    #for InputSet in range(self.NodeLabelCoords.__len__()):
    #  self.SortedDeletionLandmarkLabels.append(([],[]))
    #  self.DeletionLandmarkLabels.append(([],[]))
    #  self.MisplacedLandmarkLabels.append(([],[]))
    
    if self.DisplacementFraction > 0:
      self.SelectDisplacementPoints()
      self.DisplacePoints()
    
    if self.DeletionFraction > 0:
      self.SelectDeletionPoints(PreserveBoundaries)
      self.DeletePoints()
    
    if self.NoiseStdDev > 0:  # Check prevents division by zero if no noise is wanted
      self.AddNoise()
    
    # Create new FiducialMarkupNodes
    NewMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
    
    # Populate it with degraded landmark set
    for i, InputPoint in enumerate(self.NodeLabelCoords): #range(CurrentLandmarkSet.__len__()):
      NewMarkupsNode.AddFiducialFromArray(InputPoint[1])
      NewMarkupsNode.SetNthFiducialLabel(i, InputPoint[0])
           
    return NewMarkupsNode
   
  def SelectDisplacementPoints(self):
    #for InputSet in range(self.NodeLabelCoords.__len__()):
    DisplacementAmount = (int)(self.NodeLabelCoords.__len__() * self.DisplacementFraction)
    # could go into infinite loop if there are very few points?
    while (DisplacementAmount > 0):
      DisplacementIndex = (int)(np.random.uniform(0,self.NodeLabelCoords.__len__()))
      DisplacementLabel = self.NodeLabelCoords[DisplacementIndex][0]
      while((DisplacementIndex in self.MisplacedLandmarkLabels[1]) or (DisplacementIndex in self.DeletionLabelIndex[1])):
        DisplacementIndex = (int)(np.random.uniform(0,self.NodeLabelCoords.__len__()))
        DisplacementLabel = self.NodeLabelCoords[DisplacementIndex][0]

      self.MisplacedLandmarkLabels[0].append(DisplacementLabel)
      self.MisplacedLandmarkLabels[1].append(DisplacementIndex)
      DisplacementAmount -= 1
      
  
  def SelectDeletionPoints(self, PreserveBoundaries=False):
    #for InputSet in range(self.NodeLabelCoords.__len__()):
    DeletionAmount = (int)(self.NodeLabelCoords.__len__() * self.DeletionFraction)
    # could go into infinite loop if there are very few points?
    while (DeletionAmount > 0):
    
      if PreserveBoundaries:  
        # Preserving the boundaries means not selecting the first or last two points for deletion
        DeletionIndex = (int)(np.random.uniform(2, self.NodeLabelCoords.__len__()-2))
        
        # If point (index) was already selected for deletion or displacement, find one that wasn't
        while((DeletionIndex in self.DeletionLabelIndex[1]) or (DeletionIndex in self.MisplacedLandmarkLabels[1])):
          DeletionIndex = (int)(np.random.uniform(2,self.NodeLabelCoords.__len__()-2))
          #DeletionLabel = self.NodeLabelCoords[DeletionIndex][0]
      
      else:   # IF we can delete from the corner points as well
        DeletionIndex = (int)(np.random.uniform(0,self.NodeLabelCoords.__len__()))
        
        # If point (index) was already selected for deletion or displacement, find one that wasn't
        while((DeletionIndex in self.DeletionLabelIndex[1]) or (DeletionIndex in self.MisplacedLandmarkLabels[1])):
          DeletionIndex = (int)(np.random.uniform(0,self.NodeLabelCoords.__len__()))
          #DeletionLabel = self.NodeLabelCoords[DeletionIndex][0]
      
      # Get point label corresponding to deletion index, store both
      DeletionLabel = self.NodeLabelCoords[DeletionIndex][0]
      print "Will delete: ", DeletionLabel
      # This search is blind - pretty slow
      self.DeletionLabelIndex[0].append(DeletionLabel)
      self.DeletionLabelIndex[1].append(DeletionIndex)
      
      # One fewer point to delete
      DeletionAmount -= 1
        
  def DisplacePoints(self):       # Systematically displaces points laterally outwards, then forwards (anteriorally), attempting to emulate accidental placement on a rib
    #for InputSet in range(self.LandmarkPointSets.__len__()):
    CurrentLandmarkSet = self.NodeLabelCoords
    for i, DisplacementLabel in enumerate(self.MisplacedLandmarkLabels[0]):       # FOR each point in the landmark set having points displaced
      DisplacementPoint = CurrentLandmarkSet[self.MisplacedLandmarkLabels[1][(self.MisplacedLandmarkLabels[0]).index(DisplacementLabel)]]
      DisplacementIndex = CurrentLandmarkSet.index(DisplacementPoint)
      
      # Vector pointing from neighbor to displacement point
      RightLeftVector = self.EstimateRightLeftVector(DisplacementPoint, CurrentLandmarkSet, self.RightLeftScale)
      
      # Make sure vector points outwards
      if(DisplacementPoint[0][-1] == "L" and RightLeftVector[0] > 0):
        RightLeftVector = [-1.0*RightLeftVector[dim] for dim in range(3)]

      if(DisplacementPoint[0][-1] == "R" and RightLeftVector[0] < 0):
        RightLeftVector = [-1.0*RightLeftVector[dim] for dim in range(3)]

      # Also need a vector pointing in anterior direction
      AntPostVector = self.EstimateAntPostVector(DisplacementPoint, CurrentLandmarkSet, self.AntPostScale)
      
      # Displace the point
      self.NodeLabelCoords[DisplacementIndex][1] = [(self.NodeLabelCoords[DisplacementIndex][1][dim] + RightLeftVector[dim] + AntPostVector[dim]) for dim in range(3)]
    
  def DeletePoints(self):
    #for InputSet in range(self.LandmarkPointSets.__len__()):
    #CurrentLandmarkSet = self.NodeLabelCoords
    for i, (Label,Index) in enumerate(zip(self.DeletionLabelIndex[0], self.DeletionLabelIndex[1])): 
      Labels = [self.NodeLabelCoords[j][0] for j in range(len(self.NodeLabelCoords))]
      DeletionIndex = Labels.index(Label)
      self.NodeLabelCoords.__delitem__(DeletionIndex)
    
  # Add noise to point locations 
  #
  # MAKE SURE YOU CORRECTLY DESCRIBE THIS NOISE'S STATISTICS
  #
  def AddNoise(self): 
    #for InputSet in range(self.LandmarkPointSets.__len__()):
    CurrentLandmarkSet = self.NodeLabelCoords
    for InputPoint in range(CurrentLandmarkSet.__len__()):      # FOR each landmark point in the landmark set
      CurrentLandmarkPoint = CurrentLandmarkSet[InputPoint][1]
      for dim in range(3):                                     # FOR each of the point's spatial dimensions, noise must be added seperately
        CandidateNoise = np.random.uniform(self.MaxNoiseMM)
        while (self.CumulativeDistribution(self.NoiseStdDev, CandidateNoise) > np.random.uniform()):
          # We accept CandidateNoise probabilistically
          CandidateNoise = np.random.uniform(self.MaxNoiseMM)
        if (np.random.uniform() > 0.5):             # Random 50% chance to have noise reflected
          CandidateNoise = (-1.0) * CandidateNoise  # Noise must operate in both directions
        CurrentLandmarkPoint[dim] += CandidateNoise    
          
  def NormalDistribution(self, StdDev, Arg):
    return (1.0/(np.sqrt(2*(StdDev*StdDev)*np.pi))) * (np.exp((-(Arg*Arg))/(2*StdDev*StdDev)))
   
  def CumulativeDistribution(self, StdDev, Arg):      # Returns probablity that uniformly random measurement taken from Arg domain is less than Arg, assuming mean=0
    import math
    return (1.0/2) * (1 + math.erf(Arg / (StdDev * math.sqrt(2))))
   
  # returns vector pointing from symmetric neighbor to PointToAnchor, estimates it from nearby points if symmetric neighbor is missing
  def EstimateRightLeftVector(self, PointToAnchor, LabelsCoords, RightLeftScale):
    for LabelPoint in LabelsCoords:
      if((PointToAnchor[0][:-1] == LabelPoint[0][:-1]) and (PointToAnchor[0] != LabelPoint[0])):
        # IF the PointToAnchor has a symmetric partner
        
        RightLeftVector = [(PointToAnchor[1][dim] - LabelPoint[1][dim]) for dim in range(3)]
        RightLeftVectorLength = np.linalg.norm(RightLeftVector)

        # Normalize into unit vector
        RightLeftUnitVector = [(RightLeftVector[dim] / RightLeftVectorLength) for dim in range(3)]

        # Apply desired length scale
        RightLeftDisplacementVector = [(RightLeftUnitVector[dim] * RightLeftScale)]
        
        #print RightLeftDisplacementVector
        return RightLeftDisplacementVector
       
      else:
        # IF that PointToAnchor has no symmetric neighbor
        print "Error - spine sides not symmetric?"
        print " Could not find neighbor for point ", LabelPoint[0]
        print " Returning zero vector - RESULTS INVALID"
        return [0,0,0]
  
  def EstimateAntPostVector(self, PointToAnchor, LabelsCoords, AntPostScale):    # Anterior-Posterior vector determined by cross-product of right-left and superior-inferior vectors
    # Need  RightLeftVector for cross-product
    RightLeftVector = self.EstimateRightLeftVector(PointToAnchor, LabelsCoords, 1)

    # Find index of 
    DisplacementPointIndex = LabelsCoords[:].index(PointToAnchor)
    
    # Inititialize as-of-yet-unknown vectors
    AntPostVector = [0, 0, 0]
    InfSupVector = [0, 0, 0]
    
    if (not self.IsTopPoint(PointToAnchor, LabelsCoords)):      # IF the point to have a corresponding anchor point added is NOT at the top of the spine
      AbovePoint = self.ReturnNeighborAbove(PointToAnchor, LabelsCoords)
      SupVector = [PointToAnchor[1][dim] - AbovePoint[1][dim] for dim in range(3)]

      if(not self.IsBottomPoint(PointToAnchor, LabelsCoords)):  # IF the point is not at the bottom of the spine (we already know it's not at the top)
        BelowPoint = self.ReturnNeighborBelow(PointToAnchor, LabelsCoords)
        InfVector = [BelowPoint[1][dim] - PointToAnchor[1][dim] for dim in range(3)]
        
        # Compute Anterior-Posterior vector as cross product
        AntPostVector = np.cross(RightLeftVector, InfSupVector)
        
        # Normalize, and reflect into anterior direction, if necessary
        AntPostVectorLength = np.linalg.norm(AntPostVector)
        if(AntPostVector[1] < 0):                                                                 # IF the offset is into the posterior direction (RAS = 012)
          AntPostUnitVector = [((-1*AntPostVector[dim])/AntPostVectorLength) for dim in range(3)] # THEN reflect and normalize the vector
        else:           
          AntPostUnitVector = [(AntPostVector[dim]/AntPostVectorLength) for dim in range(3)]      # Simply normalize the vector

        AnchorPointOffsetVector = AntPostUnitVector * AntPostScale
        return AnchorPointOffsetVector
        
      else:                                                     # IF displacement point is at bottom, but not top
        AntPostVector = np.cross(RightLeftVector, SupVector)
        
        # Normalize
        AntPostVectorLength = np.linalg.norm(AntPostVector)
        if(AntPostVector[1] < 0):                               # IF the offset is into the posterior direction (RAS = 012)
          AntPostUnitVector = [((-1.0*AntPostVector[dim])/AntPostVectorLength) for dim in range(3)]   # THEN normalize and reflect
        else:                                                   # IF offset vector already points in anterior direction, only need to normalize
          AntPostUnitVector = [(AntPostVector[dim]/AntPostVectorLength) for dim in range(3)]   # THEN normalize

        AnchorPointOffsetVector = AntPostUnitVector * AntPostScale
        return AnchorPointOffsetVector
    
    else:                                                       # IF displacement point is at top, not necessarily bottom
      if(self.IsBottomPoint(PointToAnchor, LabelsCoords)):  # PointToAnchor is only one on its side
        print "ERROR - only one point on one side of spine, cannot compute anterior direction"
        print "   returning zero-AntPostVector"
        return [0,0,0]
        
      else: # PointToAnchor is top point and not bottom point
        BelowPoint = self.ReturnNeighborBelow(PointToAnchor, LabelsCoords)
        InfVector = [BelowPoint[1][0] - PointToAnchor[1][0], BelowPoint[1][1] - PointToAnchor[1][1], BelowPoint[1][2] - PointToAnchor[1][2]]
        AntPostVector = np.cross(RightLeftVector, InfVector)
        
        # Normalize
        AntPostVectorLength = np.linalg.norm(AntPostVector)
        if(AntPostVector[1] < 0):   # if the offset is into the posterior direction (RAS = 012)
          AntPostUnitVector = [((-1.0*AntPostVector[dim])/AntPostVectorLength) for dim in range(3)]   # THEN reflect and normalize the vector
        else:                                                 # IF offset vector already points in anterior direction, only need to normalize
          AntPostUnitVector = [(AntPostVector[dim]/AntPostVectorLength) for dim in range(3)]    # THEN just normalize
        
        AnchorPointOffsetVector = AntPostUnitVector * AntPostScale
        return AnchorPointOffsetVector
    
  # Is___Point functions used to check for boundary conditions, to avoid referencing non-existent points in creating vectors  
  def IsTopPoint(self, Point, PointSet):
    LandmarkLabels = []
    for LandmarkPoint in range(PointSet.__len__()):
      LandmarkLabels.append(PointSet[LandmarkPoint][0]) 
    PointIndex = LandmarkLabels.index(Point[0])
    for PotentialNeighbor in range(0,PointIndex):
      if(PointSet[PotentialNeighbor][0][-1] == Point[0][-1]):   # If 'L' == 'L' or 'R' == 'R'
        return False
    return True
        
  def IsBottomPoint(self, Point, PointSet):
    LandmarkLabels = []
    for LandmarkPoint in range(PointSet.__len__()):
      LandmarkLabels.append(PointSet[LandmarkPoint][0]) 
    PointIndex = LandmarkLabels.index(Point[0])
    for PotentialNeighbor in range(PointIndex+1,PointSet.__len__()):
      if(PointSet[PotentialNeighbor][0][-1] == Point[0][-1]):   # If 'L' == 'L' or 'R' == 'R'
        return False
    return True
     

  # It is assumed that Point is known to have this neighbor if these functions are called
  def ReturnNeighborAbove(self, Point, PointSet):
    LandmarkLabels = []
    for LandmarkPoint in range(PointSet.__len__()):
      LandmarkLabels.append(PointSet[LandmarkPoint][0]) 
    PointIndex = LandmarkLabels.index(Point[0])
    PotentialNeighbor = PointSet[PointIndex-1]
    PotentialNeighborIndex = PointIndex - 1
    while(PotentialNeighbor[0][-1] != Point[0][-1]):
      PotentialNeighborIndex -= 1
      PotentialNeighbor = PointSet[PotentialNeighborIndex]
    # ASSERT PotentialNeighbor.label[-1] == Point.label[-1]
    return PotentialNeighbor
    
  def ReturnNeighborBelow(self, Point, PointSet):
    LandmarkLabels = []
    for LandmarkPoint in range(PointSet.__len__()):
      LandmarkLabels.append(PointSet[LandmarkPoint][0]) 
    PointIndex = LandmarkLabels.index(Point[0])
    PotentialNeighbor = PointSet[PointIndex + 1]
    PotentialNeighborIndex = PointIndex + 1
    while(PotentialNeighbor[0][-1] != Point[0][-1]):
      PotentialNeighborIndex += 1
      PotentialNeighbor = PointSet[PotentialNeighborIndex]
    # ASSERT PotentialNeighbor.label[-1] == Point.label[-1]
    return PotentialNeighbor

# Assumed to be used on undegraded LandmarkSet    
class CalculateAnglesLogic(ScriptedLoadableModuleLogic):
  
  def __init__(self):
    self.MarkupsNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    self.LandmarkLabels = []
    self.LandmarkPoints = []
    for LandmarkSet in range(self.MarkupsNodes.__len__()):
      CurrentLandmarkSet = self.MarkupsNodes.__getitem__(LandmarkSet)
      self.LandmarkLabels.append([])
      self.LandmarkPoints.append([])
      for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
        self.LandmarkLabels[LandmarkSet].append(CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint))
        self.LandmarkPoints[LandmarkSet].append(CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint,0))
    self.Angles = []
    self.MaxAngles = []
    #self.MaxVertebrae = ""
      
  def CalculateAngles(self):
    for LandmarkSet in range(self.MarkupsNodes.__len__()):
      CurrentLandmarkSet = self.MarkupsNodes.__getitem__(LandmarkSet)
      self.Angles.append([])
      print ""
      print CurrentLandmarkSet.GetName()
      for Vertebra in range(0, CurrentLandmarkSet.GetNumberOfFiducials(), 2):
        VertebraLeft = CurrentLandmarkSet.GetMarkupPointVector(Vertebra, 0)
        VertebraRight = CurrentLandmarkSet.GetMarkupPointVector(Vertebra + 1, 0)
        LtoRVector = [VertebraRight[0] - VertebraLeft[0], 0,VertebraRight[2] - VertebraLeft[2]]
        LtoRVectorLength = self.FindVectorLength(LtoRVector)
        for dim in range(3):
          LtoRVector[dim] = LtoRVector[dim] / LtoRVectorLength
        Angle = np.arccos(np.dot(LtoRVector,[1,0,0])) * (180/np.pi)
        if(VertebraLeft[2] > VertebraRight[2]):
          Angle = (-1)*Angle  
        self.Angles[LandmarkSet].append(Angle)
        print Angle
    return self.Angles
      
  # returns maximum angle between any two vertbrae, measured from the left transvserse process 
  def FindMaxCoronalAngles(self):
    self.MinLabels = []
    self.MaxLabels = []
    self.MaxAngles = []
    for LandmarkSet in range(self.MarkupsNodes.__len__()):
      CurrentLandmarkSet = self.MarkupsNodes.__getitem__(LandmarkSet)
      MinAngle = 180
      MaxAngle = -180
      MinVertebra = ''
      MaxVertebra = ''
      print CurrentLandmarkSet.GetName() + "   - Finding angles"
      for Vertebra in range(0, CurrentLandmarkSet.GetNumberOfFiducials(), 2):
        if(self.Angles[LandmarkSet][Vertebra/2] < MinAngle):
          MinVertebra = ''
          it = 0
          while not (self.LandmarkLabels[LandmarkSet][Vertebra][it] == "R" or (self.LandmarkLabels[LandmarkSet][Vertebra][it] == "L"and it > 0)):
            MinVertebra += self.LandmarkLabels[LandmarkSet][Vertebra][it]
            print "Min: " + MinVertebra
            it += 1
          MinAngle = self.Angles[LandmarkSet][Vertebra/2]
        if(self.Angles[LandmarkSet][Vertebra/2] > MaxAngle):
          MaxVertebra = ''
          it = 0
          while not (self.LandmarkLabels[LandmarkSet][Vertebra][it] == "R" or (self.LandmarkLabels[LandmarkSet][Vertebra][it] == "L" and it > 0)):
            MaxVertebra += self.LandmarkLabels[LandmarkSet][Vertebra][it]
            print "Max: " + MaxVertebra
            it += 1  
          MaxAngle = self.Angles[LandmarkSet][Vertebra/2]
      #print MinVertebra + " " + MaxVertebra
      if (MinVertebra[0] == "T" and MaxVertebra[0] == "L") or (MinVertebra[0] == MaxVertebra[0] and int(MinVertebra[1:]) < int(MaxVertebra[1:])):
        # The vertebra supposed to be inferior is on top
        Vertebrae = [MaxVertebra, MinVertebra]
        Angles = [MaxAngle, MinAngle]
        MaxVertebra = Vertebrae[1]
        MinVertebra = Vertebrae[0]
        MaxAngle = Angles[1]
        MinAngle = Angles[0]
      self.MinLabels.append(MinVertebra)
      self.MaxLabels.append(MaxVertebra)
      self.MaxAngle = MaxAngle - MinAngle
      self.MaxAngles.append(self.MaxAngle)
    self.MaxVertebrae = (self.MinLabels, self.MaxLabels)
    return (self.MaxAngles, self.MaxVertebrae)
  
  def FindVectorLength(self, Vector):
    VectorNorm = 0 # initially
    for Elem in range(Vector.__len__()):
      VectorNorm += (Vector[Elem]) ** 2
    VectorNorm = np.sqrt(VectorNorm)
    return VectorNorm
    
class DegradeTransverseProcessesTest(ScriptedLoadableModuleTest):
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
    self.test_DegradeTransverseProcesses1()

  def test_DegradeTransverseProcesses1(self):
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
    logic = DegradeTransverseProcessesLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
