import os
import unittest
import vtk, qt, ctk, slicer, numpy
from slicer.ScriptedLoadableModule import *
import logging

#
# ExtendSpine
#

class ExtendSpine(ScriptedLoadableModule):

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Extend Spine"
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

class ExtendSpineWidget(ScriptedLoadableModuleWidget):

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
    # Extend Button
    #
    self.ExtendButton = qt.QPushButton("Extend transverse process landmarks")
    self.ExtendButton.toolTip = "Extends the spine by adding landmarks for vertebrae at the same angle as the top and bottom most."
    self.ExtendButton.enabled = True
    parametersFormLayout.addRow(self.ExtendButton)
    
    #
    # Combine Button
    #
    self.CombineButton = qt.QPushButton("Combine degraded and extended data")
    self.CombineButton.toolTip = "Replaces ideal points from extended set with degraded data where the correspond"
    self.CombineButton.enabled = True
    parametersFormLayout.addRow(self.CombineButton)
    
    #
    # Save Button
    #
    self.saveExtendedPointSets = qt.QPushButton("Save extended point sets")
    self.saveExtendedPointSets.toolTip = "Saves transverse process locations to .csv file containing their names and positions."
    self.saveExtendedPointSets.enabled = True
    parametersFormLayout.addRow(self.saveExtendedPointSets)
    
    # reload button
    self.reloadButton = qt.QPushButton('Restart module')
    self.layout.addWidget(self.reloadButton)
    
    # connections
    self.ExtendButton.connect('clicked(bool)', self.onExtendButton)
    self.CombineButton.connect('clicked(bool)', self.onCombineButton)
    self.reloadButton.connect('clicked(bool)',self.onReloadButton)
    self.saveExtendedPointSets.connect('clicked(bool)', self.onSaveButton)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    #self.onSelect()

  def cleanup(self):
    pass

  #def onSelect(self):
   # self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onExtendButton(self):
    logic = ExtendSpineLogic()
    logic.ExtendSpines()
    
  def onCombineButton(self):
    logic = CombineSetsLogic()
    logic.CombineDataSets()
    
  def onReloadButton(self):
    self.MarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for PointSet in range(self.MarkupNodes.__len__()):
      if (self.MarkupNodes.__getitem__(PointSet).GetNthFiducialLabel(0)[-1] == "-C"):
        slicer.mrmlScene.RemoveNode(self.MarkupNodes.__getitem__(PointSet))
      if (self.MarkupNodes.__getitem__(PointSet).GetNthFiducialLabel(0)[-1] == "-N"):
        slicer.mrmlScene.RemoveNode(self.MarkupNodes.__getitem__(PointSet))
    slicer.util.reloadScriptedModule('ExtendSpine')
    
  # Assumes CalculateAngles of DegradeTransverseProcesses has been done
  def onSaveButton(self):
    import csv
   # OriginalDataOutput = qt.QFileDialog.getSaveFileName(0, "Unextended point sets", "", "CSV File (*.csv)")
    #ModifiedDataOutput = qt.QFileDialog.getSaveFileName(0, "Extended point sets", "", "CSV File (*.csv)")
    CombinedDataOutput = qt.QFileDialog.getSaveFileName(0, "Combined point sets", "", "CSV File (*.csv)")
    self.MarkupNodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    AnglesTable = slicer.util.getModuleGui('DegradeTransverseProcesses').children()[3].children()[1]
    #print (AnglesTable.item(0,1)[0])
    MaxAngles = []
    for InputSet in range(AnglesTable.rowCount):
      MaxAngles.append(float(AnglesTable.item(InputSet,1).data(0)))
    # with open(OriginalDataOutput, 'wb') as csvfile:
      # writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      # writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      # for MarkupsNode in range(self.MarkupNodes.__len__()):
        # CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        # if(CurrentLandmarkSet.GetName()[-1] != "-C"):
          # # If the landmark set is not an extended one
          # #writer.writerow(['Max angle:', str(MaxAngles[MarkupsNode]), CurrentLandmarkSet.GetName(), '']) 
          # writer.writerow([CurrentLandmarkSet.GetName(), '', '', '']) 
          # for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            # CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            # writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      # writer.writerow(['EOF', '', '', ''])
            
    # with open(ModifiedDataOutput, 'wb') as csvfile:
      # writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      # writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      # CompleteSetNum = 0
      # for MarkupsNode in range(self.MarkupNodes.__len__()):
        # CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        # if(CurrentLandmarkSet.GetName()[-1] == "-C"):
          # # If the landmark set is a Completed one
          # #writer.writerow(['Max angle:', str(MaxAngles[CompleteSetNum]), CurrentLandmarkSet.GetName(), '']) 
          # CompleteSetNum += 1
          # writer.writerow([CurrentLandmarkSet.GetName(), '', '', '']) 
          # for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            # CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            # writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      # writer.writerow(['EOF', '', '', ''])
      
    with open(CombinedDataOutput, 'wb') as csvfile:
      writer = csv.writer(csvfile, delimiter=',', quotechar='|')
      writer.writerow(['Landmark', 'RL', 'AP', 'SI'])
      CompleteSetNum = 0
      self.MarkupNodes = sorted(self.MarkupsNodes, key = lambda DataSet: DataSet.GetName())
      for MarkupsNode in range(self.MarkupNodes.__len__()):
        CurrentLandmarkSet = self.MarkupNodes.__getitem__(MarkupsNode)
        if(CurrentLandmarkSet.GetName()[-1] == "-N"):
          # If the landmark set is a Final one
          writer.writerow(['Max angle:', str(MaxAngles[CompleteSetNum]), CurrentLandmarkSet.GetName(), '']) 
          CompleteSetNum += 1
          for LandmarkPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
            CurrentPoint = CurrentLandmarkSet.GetMarkupPointVector(LandmarkPoint, 0)
            writer.writerow([CurrentLandmarkSet.GetNthFiducialLabel(LandmarkPoint), str(CurrentPoint[0]), str(CurrentPoint[1]), str(CurrentPoint[2])])
      writer.writerow(['EOF', '', '', ''])
    
#
# ExtendSpineLogic
#

class ExtendSpineLogic(ScriptedLoadableModuleLogic):

  def __init__(self):
    
    self.LandmarkPointSets = []      # Will contain tuples: (TrXFiducial point labels, TrXFiducial### point sets)
    self.PointSpread = 30       # Length of vector placing points from outermost vertebrae (mm)
    self.CompleteLabelList = ['T1L','T1R','T2L','T2R','T3L','T3R','T4L','T4R','T5L','T5R','T6L','T6R','T7L','T7R','T8L','T8R','T9L','T9R','T10L','T10R','T11L','T11R','T12L','T12R','L1L','L1R','L2L','L2R','L3L','L3R','L4L','L4R','L5L','L5R']
      
    
  def ExtendSpines(self):
    self.InputData = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for InputSet in range(self.InputData.__len__()):
      if self.InputData.__getitem__(InputSet).GetName()[-1] != "~" and self.InputData.__getitem__(InputSet).GetName()[-1] != "-N":
        CurrentLandmarkSet = self.InputData.__getitem__(InputSet)
        self.LandmarkPointSets.append([])
        print " "   #empty line
        print "Landmark set #" + str(InputSet)
        for InputPoint in range(CurrentLandmarkSet.GetNumberOfFiducials()):
          self.LandmarkPointSets[InputSet].append((CurrentLandmarkSet.GetNthFiducialLabel(InputPoint),CurrentLandmarkSet.GetMarkupPointVector(InputPoint,0)))
          print self.LandmarkPointSets[InputSet][InputPoint]
        
    self.ExtendSuperiorSpines()
    self.ExtendInferiorSpines()
    
    # Create new FiducialMarkupNodes for Slicer scene
    for InputSet in range(self.LandmarkPointSets.__len__()):
      NewMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewMarkupsNode.SetName(self.InputData.__getitem__(InputSet).GetName() + "-C")
      CurrentLandmarkSet = self.LandmarkPointSets[InputSet]
      for InputPoint in range(CurrentLandmarkSet.__len__()):
        CurrentLandmarkPoint = CurrentLandmarkSet[InputPoint][1]
        NewMarkupsNode.AddFiducial(CurrentLandmarkPoint[0], CurrentLandmarkPoint[1], CurrentLandmarkPoint[2])
        NewPointLabel = self.LandmarkPointSets[InputSet][InputPoint][0]
        NewMarkupsNode.SetNthFiducialLabel(InputPoint, NewPointLabel)
      slicer.mrmlScene.AddNode(NewMarkupsNode)
           
    return True
   
  def ExtendSuperiorSpines(self):
    for i, CurrentLandmarkSet in enumerate(self.LandmarkPointSets):
      TopLeftLabel = CurrentLandmarkSet[0][0]
      if(TopLeftLabel != 'T1L'):
        TopLeftLabel = self.LandmarkPointSets[i][0][0]
        TopLeftIndex = self.CompleteLabelList.index(TopLeftLabel)  
        if TopLeftIndex == 2:
          InfSupVector = [0,0,0]
          NewRightPoint = [0,0,0]
          NewLeftPoint = [0,0,0]
          for dim in range(3):
            InfSupVector[dim] = (CurrentLandmarkSet[1][1][dim] - CurrentLandmarkSet[3][1][dim])/2
            InfSupVector[dim] += (CurrentLandmarkSet[0][1][dim] - CurrentLandmarkSet[2][1][dim])/2
            NewRightPoint[dim] = CurrentLandmarkSet[1][1][dim] + InfSupVector[dim]
            NewLeftPoint[dim] = CurrentLandmarkSet[0][1][dim] + InfSupVector[dim]
          self.LandmarkPointSets[i].insert(0,(self.CompleteLabelList[TopLeftIndex-1], NewRightPoint))
          self.LandmarkPointSets[i].insert(0,(self.CompleteLabelList[TopLeftIndex-2], NewLeftPoint))
        while(TopLeftIndex > 2):
          TopLeftLabel = self.LandmarkPointSets[i][0][0]
          TopLeftIndex = self.CompleteLabelList.index(TopLeftLabel)  
          InfSupVector = [0,0,0]
          NewRightPoint = [0,0,0]
          NewLeftPoint = [0,0,0]
          for dim in range(3):
            InfSupVector[dim] = (CurrentLandmarkSet[1][1][dim] - CurrentLandmarkSet[3][1][dim])/2
            InfSupVector[dim] += (CurrentLandmarkSet[0][1][dim] - CurrentLandmarkSet[2][1][dim])/2
            NewRightPoint[dim] = CurrentLandmarkSet[1][1][dim] + InfSupVector[dim]
            NewLeftPoint[dim] = CurrentLandmarkSet[0][1][dim] + InfSupVector[dim]
          self.LandmarkPointSets[i].insert(0,(self.CompleteLabelList[TopLeftIndex-1], NewRightPoint))
          self.LandmarkPointSets[i].insert(0,(self.CompleteLabelList[TopLeftIndex-2], NewLeftPoint))

  def ExtendInferiorSpines(self):
    for i, CurrentLandmarkSet in enumerate(self.LandmarkPointSets):
      BottomRightLabel = CurrentLandmarkSet[-1][0]
      if(BottomRightLabel != "L5R"):
        BottomRightLabel = self.LandmarkPointSets[i][-1][0]
        BottomRightIndex = self.CompleteLabelList.index(BottomRightLabel) 
        if BottomRightIndex == 31:
          InfSupVector = [0,0,0]
          NewRightPoint = [0,0,0]
          NewLeftPoint = [0,0,0]
          for dim in range(3):
            InfSupVector[dim] = (CurrentLandmarkSet[-1][1][dim] - CurrentLandmarkSet[-3][1][dim])/2
            InfSupVector[dim] += (CurrentLandmarkSet[-2][1][dim] - CurrentLandmarkSet[-4][1][dim])/2
            NewRightPoint[dim] = CurrentLandmarkSet[-1][1][dim] + InfSupVector[dim]
            NewLeftPoint[dim] = CurrentLandmarkSet[-2][1][dim] + InfSupVector[dim]
          self.LandmarkPointSets[i].append((self.CompleteLabelList[BottomRightIndex+1], NewLeftPoint))
          self.LandmarkPointSets[i].append((self.CompleteLabelList[BottomRightIndex+2], NewRightPoint))
        while(BottomRightIndex < 31):
          BottomRightLabel = self.LandmarkPointSets[i][-1][0]
          BottomRightIndex = self.CompleteLabelList.index(BottomRightLabel)  
          InfSupVector = [0,0,0]
          NewRightPoint = [0,0,0]
          NewLeftPoint = [0,0,0]
          for dim in range(3):
            InfSupVector[dim] = (CurrentLandmarkSet[-1][1][dim] - CurrentLandmarkSet[-3][1][dim])/2
            InfSupVector[dim] += (CurrentLandmarkSet[-2][1][dim] - CurrentLandmarkSet[-4][1][dim])/2
            NewRightPoint[dim] = CurrentLandmarkSet[-1][1][dim] + InfSupVector[dim]
            NewLeftPoint[dim] = CurrentLandmarkSet[-2][1][dim] + InfSupVector[dim]
          self.LandmarkPointSets[i].append((self.CompleteLabelList[BottomRightIndex+1], NewLeftPoint))
          self.LandmarkPointSets[i].append((self.CompleteLabelList[BottomRightIndex+2], NewRightPoint))

   
  # returns vector pointing from symmetric neighbor to DisplacementPoint, estimates it from nearby points if symmetric neighbor is missing
  def EstimateRightLeftVector(self, DisplacementPoint, LandmarkSet, RightLeftScale):
    for LabelPoint in LandmarkSet:
      if((DisplacementPoint[0][:-1] == LabelPoint[0][:-1]) and (DisplacementPoint[0] != LabelPoint[0])):
        # The DisplacementPoint has a symmetric partner
        
        RightLeftVector = [0,0,0]
        RightLeftVectorNorm = 0
        for dim in range(3):
          RightLeftVector[dim] = DisplacementPoint[1][dim] - LabelPoint[1][dim]
          RightLeftVectorNorm += (RightLeftVector[dim]) ** 2
        LeftRightVectorNorm = numpy.sqrt(RightLeftVectorNorm)
        for dim in range(3):
          RightLeftVector[dim] = RightLeftScale * (RightLeftVector[dim] / RightLeftVectorNorm)
        #print RightLeftVector
        return RightLeftVector
    # ASSERT that DisplacementPoint has no symmetric neighbor
  
  def EstimateAntPostVector(self, DisplacementPoint, LandmarkSet, AntPostScale):
    # Need normalized RightLeftVector for cross-product
    RightLeftVector = self.EstimateRightLeftVector(DisplacementPoint, LandmarkSet, 1)
    LandmarkLabels = []
    for LandmarkPoint in range(LandmarkSet.__len__()):
      LandmarkLabels.append(LandmarkSet[LandmarkPoint][0]) 
    DisplacementPointIndex = LandmarkLabels.index(DisplacementPoint[0])
    AntPostVector = [0, 0, 0]
    InfSupVector = [0, 0, 0]
    if (not self.IsTopPoint(DisplacementPoint, LandmarkSet)):
      AbovePoint = self.ReturnNeighborAbove(DisplacementPoint, LandmarkSet)
      SupVector = [DisplacementPoint[1][0] - AbovePoint[1][0], DisplacementPoint[1][1] - AbovePoint[1][1], DisplacementPoint[1][2] - AbovePoint[1][2]]
      LengthSupVector = self.FindVectorLength(SupVector)
      for dim in range(3):
        SupVector[dim] = SupVector[dim] / LengthSupVector
      if(not self.IsBottomPoint(DisplacementPoint, LandmarkSet)):
        BelowPoint = self.ReturnNeighborBelow(DisplacementPoint, LandmarkSet)
        InfVector = [BelowPoint[1][0] - DisplacementPoint[1][0], BelowPoint[1][1] - DisplacementPoint[1][1], BelowPoint[1][2] - DisplacementPoint[1][2]]
        LengthInfVector = self.FindVectorLength(InfVector)
        for dim in range(3):
          InfVector[dim] = InfVector[dim] / LengthInfVector
          InfSupVector[dim] = (SupVector[dim] + InfVector[dim]) / 2
        AntPostVector = numpy.cross(RightLeftVector, InfSupVector)
        if(AntPostVector[1] < 0):   # if the offset is into the posterior direction (RAS = 012)
          for dim in range(3):
            AntPostVector[dim] = (-1) * AntPostVector[dim]  # Reflect the vector
        # Ensure normalization
        LengthAntPostVector = self.FindVectorLength(AntPostVector)
        for dim in range(3):
          AntPostVector[dim] = AntPostVector[dim] / LengthAntPostVector
        return AntPostVector * AntPostScale
      else:   # Displacement point is at bottom, but not top
        AntPostVector = numpy.cross(RightLeftVector, SupVector)
        if(AntPostVector[1] < 0):   # if the offset is into the posterior direction (RAS = 012)
          for dim in range(3):
            AntPostVector[dim] = (-1) * AntPostVector[dim]  # Reflect the vector
        # Ensure normalization
        LengthAntPostVector = self.FindVectorLength(AntPostVector)
        for dim in range(3):
          AntPostVector[dim] = AntPostVector[dim] / LengthAntPostVector
        return AntPostVector * AntPostScale
    else:   # Displacement point is at top, not necessarily bottom
      if(self.IsBottomPoint(DisplacementPoint, LandmarkSet)):  # DisplacementPoint is only one on its side
        print "ERROR - only one point on one side of spine, cannot compute anterior direction"
        print "   returning zero-AntPostVector"
        return [0,0,0]
      else: # DisplacementPoint is top point and not bottom point
        BelowPoint = self.ReturnNeighborBelow(DisplacementPoint, LandmarkSet)
        InfVector = [BelowPoint[1][0] - DisplacementPoint[1][0], BelowPoint[1][1] - DisplacementPoint[1][1], BelowPoint[1][2] - DisplacementPoint[1][2]]
        LengthInfVector = self.FindVectorLength(InfVector)
        for dim in range(3):
          InfVector[dim] = InfVector[dim] / LengthInfVector
        AntPostVector = numpy.cross(RightLeftVector, InfVector)
        if(AntPostVector[1] < 0):   # if the offset is into the posterior direction (RAS = 012)
          for dim in range(3):
            AntPostVector[dim] = (-1) * AntPostVector[dim]  # Reflect the vector
        # Ensure normalization
        LengthAntPostVector = self.FindVectorLength(AntPostVector)
        for dim in range(3):
          AntPostVector[dim] = AntPostVector[dim] / LengthAntPostVector
        return AntPostVector * AntPostScale
    
  def FindVectorLength(self, Vector):
    VectorNorm = 0 # initially
    for Elem in range(Vector.__len__()):
      VectorNorm += (Vector[Elem]) ** 2
    VectorNorm = numpy.sqrt(VectorNorm)
    return VectorNorm
    
    
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

class CombineSetsLogic(ScriptedLoadableModuleLogic):
  def __init__(self):
    self.CompleteLabelList = ['T1L','T1R','T2L','T2R','T3L','T3R','T4L','T4R','T5L','T5R','T6L','T6R','T7L','T7R','T8L','T8R','T9L','T9R','T10L','T10R','T11L','T11R','T12L','T12R','L1L','L1R','L2L','L2R','L3L','L3R','L4L','L4R','L5L','L5R']
    self.ExtendedDataSets = []
    self.DegradedDataSets = []
    self.CombinedDataSets = []
    self.AllDataSets = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for DataSet in range(self.AllDataSets.__len__()):
      if self.AllDataSets.__getitem__(DataSet).GetName()[-1] == "C":
        self.ExtendedDataSets.append(self.AllDataSets.__getitem__(DataSet))
      elif self.AllDataSets.__getitem__(DataSet).GetName()[-1] == "~":
        self.DegradedDataSets.append(self.AllDataSets.__getitem__(DataSet))
    self.ExtendedDataSets = sorted(self.ExtendedDataSets, key=lambda DataSet: DataSet.GetName())
    self.DegradedDataSets = sorted(self.DegradedDataSets, key=lambda DataSet: DataSet.GetName())
    
  def CombineDataSets(self):
    for PatientDataSet in range(self.ExtendedDataSets.__len__()):   # Assumes there is a degraded set for each extended one
      NewMarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewMarkupsNode.SetName(self.ExtendedDataSets.__getitem__(PatientDataSet).GetName() + "+")
      CurrentExtendedLandmarkNode = self.ExtendedDataSets.__getitem__(PatientDataSet)
      CurrentDegradedLandmarkNode = self.DegradedDataSets.__getitem__(PatientDataSet)
      
      TopOffset = self.CompleteLabelList.index(CurrentDegradedLandmarkNode.GetNthFiducialLabel(0))
      BottomOffset = self.CompleteLabelList.index(CurrentDegradedLandmarkNode.GetNthFiducialLabel(CurrentDegradedLandmarkNode.GetNumberOfFiducials() - 1))
      for ExtendedLandmark in range(TopOffset):
        CurrentExtendedLandmarkPoint = CurrentExtendedLandmarkNode.GetMarkupPointVector(ExtendedLandmark,0)
        NewMarkupsNode.AddFiducial(CurrentExtendedLandmarkPoint[0], CurrentExtendedLandmarkPoint[1], CurrentExtendedLandmarkPoint[2])
        NewMarkupsNode.SetNthFiducialLabel(ExtendedLandmark, CurrentExtendedLandmarkNode.GetNthFiducialLabel(ExtendedLandmark) + "+")
      for DegradedLandmark in range(CurrentDegradedLandmarkNode.GetNumberOfFiducials()):
        CurrentDegradedLandmarkPoint = CurrentDegradedLandmarkNode.GetMarkupPointVector(DegradedLandmark,0)
        NewMarkupsNode.AddFiducial(CurrentDegradedLandmarkPoint[0], CurrentDegradedLandmarkPoint[1], CurrentDegradedLandmarkPoint[2])
        NewMarkupsNode.SetNthFiducialLabel(TopOffset + DegradedLandmark, CurrentDegradedLandmarkNode.GetNthFiducialLabel(DegradedLandmark) + "+")
      for ExtendedLandmark in range(BottomOffset+1, CurrentExtendedLandmarkNode.GetNumberOfFiducials()):
        CurrentExtendedLandmarkPoint = CurrentExtendedLandmarkNode.GetMarkupPointVector(ExtendedLandmark,0)
        NewMarkupsNode.AddFiducial(CurrentExtendedLandmarkPoint[0], CurrentExtendedLandmarkPoint[1], CurrentExtendedLandmarkPoint[2])
        NewMarkupsNode.SetNthFiducialLabel(ExtendedLandmark, CurrentExtendedLandmarkNode.GetNthFiducialLabel(ExtendedLandmark) + "+")
      slicer.mrmlScene.AddNode(NewMarkupsNode)
 
class ExtendSpineTest(ScriptedLoadableModuleTest):
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
    self.test_ExtendSpine1()

  def test_ExtendSpine1(self):
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
    logic = ExtendSpineLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
