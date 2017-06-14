﻿import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
#from math import *
#import math

#
# PreProcessLandmarks
#

class PreProcessLandmarks(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "PreProcessLandmarks" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Scoliosis"]
    self.parent.dependencies = []
    self.parent.contributors = ["Ben Church - PerkLab, Queen's University",]
    self.parent.helpText = """
    This module can be used to extrapolate or interpolate incomplete scans of the
    transverse processes"""
    self.parent.acknowledgementText = """ """ # replace with organization, grant and thanks.

#
# PreProcessLandmarksWidget
#

class PreProcessLandmarksWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    self.reloadCollapsibleButton.collapsed = 1
    
    # Operations performed by calling self.logic.[operation] - self.SpineModel is for storing self.logic.PatientModel
    
    self.logic = None           # Use class variable to store logic to keep performing operations on same nodes

    # Store relevant anatomy models as class variables
    self.SpineModel = None
    self.LeftPatchNode = None
    self.RightPatchNode = None
    #self.PatchTag = "_Patch"      # Used to check ComboBoxes for nodes' corresponding patches
    
    self.kmWindowSize = 5
    self.PolynomialDegree = 4
    self.MaxUndos = 10
    self.NodeStack = []
    
    # Instantiate and connect widgets ...
    
    # User interface for incomplete landmark set repair
    RepairInterface = ctk.ctkCollapsibleButton()
    RepairInterface.text = "Markups node repair"
    self.layout.addWidget(RepairInterface)
    
    self.RepairInterfaceLayout = qt.QGridLayout(RepairInterface)
    self.RepairInterfaceLayout.setHorizontalSpacing(12)
    
    # Dropdown list to select MarkupsNode for repair
    self.AnatomySelector = slicer.qMRMLNodeComboBox()
    self.AnatomySelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.AnatomySelector.selectNodeUponCreation = True
    self.AnatomySelector.enabled  = True
    self.AnatomySelector.addEnabled = True
    self.AnatomySelector.noneEnabled = True
    self.AnatomySelector.removeEnabled = True
    self.AnatomySelector.renameEnabled = True
    self.AnatomySelector.toolTip = "Choose the incomplete patient landmark set to be repaired"
    self.RepairInterfaceLayout.addWidget(qt.QLabel("Node to repair -->"), 1, 0, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.AnatomySelector, 1, 1, 1, 3)
    self.AnatomySelector.setMRMLScene(slicer.mrmlScene)
    
    self.CategorizeLeftRightButton = qt.QPushButton("Categorize left-right landmarks")
    self.CategorizeLeftRightButton.toolTip = "Groups selected nodes landmarks into left and right landmark nodes used for repair operations"
    self.CategorizeLeftRightButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.CategorizeLeftRightButton, 2, 1, 1, 3)
    
    self.LeftSideSelector = slicer.qMRMLNodeComboBox()
    self.LeftSideSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.LeftSideSelector.selectNodeUponCreation = False
    self.LeftSideSelector.enabled  = False
    self.LeftSideSelector.addEnabled = False
    #self.LeftSideSelector.editEnabled = False
    self.LeftSideSelector.noneEnabled = True
    self.LeftSideSelector.removeEnabled = False
    self.LeftSideSelector.renameEnabled = False
    self.LeftSideSelector.toolTip = "Stores the left side of the patient's landmarks being repaired"
    self.LeftSideSelector.setMRMLScene(slicer.mrmlScene)

    self.RightSideSelector = slicer.qMRMLNodeComboBox()
    self.RightSideSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.RightSideSelector.selectNodeUponCreation = False
    self.RightSideSelector.enabled  = False
    self.RightSideSelector.addEnabled = False
    self.RightSideSelector.noneEnabled = True
    self.RightSideSelector.removeEnabled = True
    self.RightSideSelector.renameEnabled = False
    self.RightSideSelector.toolTip = "Stores the right side of the patient's landmarks being repaired"
    self.RightSideSelector.setMRMLScene(slicer.mrmlScene)
    
    self.RepairInterfaceLayout.addWidget(qt.QLabel(" "), 3, 1, 1, 3)
    
    self.RepairInterfaceLayout.addWidget(qt.QLabel("Unrepaired"), 4, 1, 1, 1)
    self.RepairInterfaceLayout.addWidget(qt.QLabel("Tentative patch"), 4, 3, 1, 1)
    
    self.LeftPatchSelector = slicer.qMRMLNodeComboBox()
    self.LeftPatchSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.LeftPatchSelector.selectNodeUponCreation = False
    self.LeftPatchSelector.enabled  = False
    self.LeftPatchSelector.addEnabled = False
    self.LeftPatchSelector.noneEnabled = True
    self.LeftPatchSelector.removeEnabled = True
    self.LeftPatchSelector.renameEnabled = False
    self.LeftPatchSelector.toolTip = "Stores the landmarks node with which the module will patch Unrepaired Left" 
    self.LeftPatchSelector.setMRMLScene(slicer.mrmlScene)

    self.RepairInterfaceLayout.addWidget(qt.QLabel("Left"), 5, 0, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.LeftSideSelector, 5, 1, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.LeftPatchSelector, 5, 3, 1, 1)

    self.RightPatchSelector = slicer.qMRMLNodeComboBox()
    self.RightPatchSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.RightPatchSelector.selectNodeUponCreation = False
    self.RightPatchSelector.enabled  = False
    self.RightPatchSelector.addEnabled = False
    self.RightPatchSelector.noneEnabled = True
    self.RightPatchSelector.removeEnabled = True
    self.RightPatchSelector.renameEnabled = False
    self.RightPatchSelector.toolTip = "Stores the landmarks node with which the module will patch Unrepaired Right" 
    self.RightPatchSelector.setMRMLScene(slicer.mrmlScene)
    
    self.RepairInterfaceLayout.addWidget(qt.QLabel("Right"), 6, 0, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.RightSideSelector, 6, 1, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.RightPatchSelector, 6, 3, 1, 1)
  
    self.RepairInterfaceLayout.addWidget(qt.QLabel(" "), 7, 1, 1, 3)
   
    self.RepairSideSelector = slicer.qMRMLNodeComboBox()
    self.RepairSideSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.RepairSideSelector.selectNodeUponCreation = False
    self.RepairSideSelector.enabled  = True
    self.RepairSideSelector.addEnabled = False
    self.RepairSideSelector.noneEnabled = True
    self.RepairSideSelector.removeEnabled = True
    self.RepairSideSelector.renameEnabled = False
    self.RepairSideSelector.toolTip = "Use to select the node containing the spine side to have a patch of points generated"
    self.RepairSideSelector.setMRMLScene(slicer.mrmlScene)
    self.RepairInterfaceLayout.addWidget(qt.QLabel("Make patch for:"), 8, 0, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.RepairSideSelector, 8, 1, 1, 1)
    
    self.AutoGenPatchButton = qt.QPushButton("AutoGen Patch")
    self.AutoGenPatchButton.toolTip = "Guess missing landmarks"
    self.AutoGenPatchButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.AutoGenPatchButton, 8, 2, 1, 1)
    
    self.ApplyPatchButton = qt.QPushButton("Apply Patch")
    self.ApplyPatchButton.toolTip = "Replace selected side node with itself merged with its patch"
    self.ApplyPatchButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.ApplyPatchButton, 8, 3, 1, 1)
    
    self.SubIntervalSelector = qt.QSpinBox()
    self.SubIntervalSelector.minimum = 0
    self.SubIntervalSelector.maximum = 10
    self.SubIntervalSelector.singleStep = 1
    self.RepairInterfaceLayout.addWidget(qt.QLabel("SubPatch:"), 10, 0, 1, 1)
    self.RepairInterfaceLayout.addWidget(self.SubIntervalSelector, 10, 1, 1, 1)
    
    self.AddPointToPatchButton = qt.QPushButton("Add Point")
    self.AddPointToPatchButton.toolTip = "Add one point to SubPatch"
    self.AddPointToPatchButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.AddPointToPatchButton, 10, 2, 1, 1)

    self.RemovePointButton = qt.QPushButton("Remove Point")
    self.RemovePointButton.toolTip = "Remove one point to SubPatch"
    self.RemovePointButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.RemovePointButton, 10, 3, 1, 1)    
   
    self.RepairInterfaceLayout.addWidget(qt.QLabel(" "), 11, 1, 1, 3)
   
    """
    self.UndoButton = qt.QPushButton("Undo")
    self.UndoButton.toolTip = "Undo up to" + str(self.MaxUndos) + " operations"
    self.UndoButton.enabled = True
    self.RepairInterfaceLayout.addWidget(self.UndoButton, 12, 0, 1, 4)
    """
   
    self.IdentifyOutliersButton = qt.QPushButton("Identify Outliers")
    self.IdentifyOutliersButton.toolTip = "Highlight apparent outlier points in selected side"
    self.IdentifyOutliersButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.IdentifyOutliersButton, 12, 0, 1, 2)
    
    self.RemoveOutliersButton = qt.QPushButton("Remove Outliers")
    self.RemoveOutliersButton.toolTip = "Remove points highlighted as outliers by Remove Outliers"
    self.RemoveOutliersButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.RemoveOutliersButton, 12, 2, 1, 2)
   
    self.RepairInterfaceLayout.addWidget(qt.QLabel(""), 13, 1, 1, 1)
    
    self.RepairNodeButton = qt.QPushButton("Repair landmarks node")
    self.RepairNodeButton.toolTip = "Combine left and right Unrepaired side nodes into RepairedNode"
    self.RepairNodeButton.enabled = False
    self.RepairInterfaceLayout.addWidget(self.RepairNodeButton, 14, 0, 1, 4)
    
    # Dropdown list to store output repaired MarkupsNode
    self.RepairedNodeStorage = slicer.qMRMLNodeComboBox()
    self.RepairedNodeStorage.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.RepairedNodeStorage.selectNodeUponCreation = False
    self.RepairedNodeStorage.enabled  = True
    self.RepairedNodeStorage.addEnabled = False
    self.RepairedNodeStorage.noneEnabled = True
    self.RepairedNodeStorage.removeEnabled = True
    self.RepairedNodeStorage.renameEnabled = False
    self.RepairedNodeStorage.toolTip = "Stores combined repaired left and right sides - repaired version of input repair node"
    self.RepairInterfaceLayout.addWidget(self.RepairedNodeStorage, 15, 0, 1, 3)
    self.RepairInterfaceLayout.addWidget(qt.QLabel(" <-- Repaired node"), 15, 3, 1, 1)
    self.RepairedNodeStorage.setMRMLScene(slicer.mrmlScene)
    
    # Reload module button
    self.reloadButton = qt.QPushButton('Reload module')
    self.RepairInterfaceLayout.addWidget(self.reloadButton, 16, 0, 1, 4)
    
    # Button connections
    self.CategorizeLeftRightButton.connect('clicked(bool)', self.OnCategorizeLeftRight)
    
    self.AutoGenPatchButton.connect('clicked(bool)', self.OnAutoGenPatchButton)
    self.AutoGenPatchButton.connect('clicked(bool)', self.OnLeftRightBasePatchChange)
    
    self.ApplyPatchButton.connect('clicked(bool)', self.OnApplyPatchButton)
    
    self.AddPointToPatchButton.connect('clicked(bool)', self.OnAddPointToPatchButton)
    self.RemovePointButton.connect('clicked(bool)', self.OnRemovePointButton)
    self.RemovePointButton.connect('clicked(bool)', self.OnLeftRightBasePatchChange)
    
    self.IdentifyOutliersButton.connect('clicked(bool)', self.OnIdentifyOutliersButton)
    self.RemoveOutliersButton.connect('clicked(bool)', self.OnRemoveOutliersButton)
    
    self.RepairNodeButton.connect('clicked(bool)', self.OnRepairButtonClicked)
    
    self.reloadButton.connect('clicked(bool)', self.OnReloadButton)
    
    #self.UndoButton.connect('clicked(bool)', self.OnUndoButtonClicked)
    
    # Node interface connections
    self.AnatomySelector.connect('currentNodeChanged(bool)', self.OnSelectedAnatomyChange)
    #self.LeftSideSelector.connect('currentNodeChanged(bool)', self.OnLeftRightBasePatchChange)
    #self.RightSideSelector.connect('currentNodeChanged(bool)', self.OnLeftRightBasePatchChange)
    #self.LeftPatchSelector.connect('currentNodeChanged(bool)', self.OnLeftRightBasePatchChange)
    #self.RightPatchSelector.connect('currentNodeChanged(bool)', self.OnLeftRightBasePatchChange)
    #self.RepairSideSelector.connect('currentNodeChanged(bool)', self.OnSelectedSideChange)
    self.RepairSideSelector.connect('currentNodeChanged(bool)', self.OnLeftRightBasePatchChange)
    
    #
    # Configuration panel
    #
    
    # Initialize
    ConfigInterface = ctk.ctkCollapsibleButton()
    ConfigInterface.text = "Module configuration"
    self.layout.addWidget(ConfigInterface)
    
    self.ConfigInterfaceLayout = qt.QGridLayout(ConfigInterface)
    self.ConfigInterfaceLayout.setHorizontalSpacing(12)
    
    # Configuration panel widgets
    self.kmWindowSizeSlider = qt.QSlider(0x1)
    self.kmWindowSizeSlider.toolTip = "Sets the number of points used in each iteration of k-means for right-left categorization"
    self.kmWindowSizeSlider.enabled = True
    self.kmWindowSizeSlider.setMinimum(2)
    self.kmWindowSizeSlider.setMaximum(12)
    self.kmWindowSizeSlider.value = 5
    self.kmWindowTextDisplay = qt.QLabel(str(self.kmWindowSizeSlider.value))
    self.ConfigInterfaceLayout.addWidget(qt.QLabel("kMeans Window"), 1, 0, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.kmWindowSizeSlider, 1, 1, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.kmWindowTextDisplay, 1, 2, 1, 1)
    
    self.PolyFitDegreeSlider = qt.QSlider(0x1)
    self.PolyFitDegreeSlider.toolTip = "Sets the number degrees of the polynomials fit to right and left R-L, A-P coordinates, as well as their curve-wise point frequencies"
    self.PolyFitDegreeSlider.enabled = True
    self.PolyFitDegreeSlider.setMinimum(1)
    self.PolyFitDegreeSlider.setMaximum(10)
    self.PolyFitDegreeSlider.value = 4
    self.PolyFitDegreeSliderTextDisplay = qt.QLabel(str(self.PolyFitDegreeSlider.value))
    self.ConfigInterfaceLayout.addWidget(qt.QLabel("PolyFit Degree"), 2, 0, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.PolyFitDegreeSlider, 2, 1, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.PolyFitDegreeSliderTextDisplay, 2, 2, 1, 1)
    
    self.SpecificitySlider = qt.QSlider(0x1)
    self.SpecificitySlider.toolTip = "Sets a scalar multiplier used to compare interval measures to curve statistics"
    self.SpecificitySlider.enabled = True
    self.SpecificitySlider.setMinimum(-10)
    self.SpecificitySlider.setMaximum(10)
    self.SpecificityTextDisplay = qt.QLabel(str(self.SpecificitySlider.value))
    self.ConfigInterfaceLayout.addWidget(qt.QLabel("AutoGen Specificity"), 3, 0, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.SpecificitySlider, 3, 1, 1, 1)
    self.ConfigInterfaceLayout.addWidget(self.SpecificityTextDisplay, 3, 2, 1, 1)
    
    # Configuration panel connections
    self.kmWindowSizeSlider.connect('sliderReleased()', self.OnkmWindowSliderChanged)
    self.PolyFitDegreeSlider.connect('sliderReleased()', self.OnPolyFitDegreeSliderChanged)
    self.SpecificitySlider.connect('sliderReleased()', self.OnSpecificitySliderChanged)
    
    # Add vertical spacer
    self.layout.addStretch(1)

  def cleanup(self):
    pass

  def OnSelectedAnatomyChange(self):
    SpineNode = self.AnatomySelector.currentNode()
    if SpineNode != None:
      self.logic = PreProcessLandmarksLogic(SpineNode, self.kmWindowSizeSlider.value, self.PolynomialDegree, self.SpecificitySlider.value)
      #self.SpecificitySlider.enabled = True
      self.SpineModel = self.logic.PatientModel
      self.CategorizeLeftRightButton.enabled = True
    else:
      self.CategorizeLeftRightButton.enabled = False
      
    self.UpdateColors()
      
  def OnLeftRightBasePatchChange(self):
  
    # Get node selectors' contents, even if they are None
    LeftBaseNode = self.LeftSideSelector.currentNode()
    RightBaseNode = self.RightSideSelector.currentNode()
    LeftPatchNode = self.LeftPatchSelector.currentNode()
    RightPatchNode = self.RightPatchSelector.currentNode()
    NodeForRepair = self.RepairSideSelector.currentNode()
    
    if NodeForRepair == None:   # Prevents other checks which could be None == None
      self.AutoGenPatchButton.enabled = False
      self.ApplyPatchButton.enabled = False
      self.AddPointToPatchButton.enabled = False
      self.RemovePointButton.enabled = False
      return    
    
    # Having a side node selected for patch generation allows patch generation, and then patch modification
    if NodeForRepair == LeftBaseNode:
      self.AutoGenPatchButton.enabled = True
      if LeftPatchNode!= None and LeftPatchNode.GetName() == LeftBaseNode.GetName() + "_Patch":   # IF we already have a patch for it, make a few more operations availible
        self.ApplyPatchButton.enabled = True
        self.AddPointToPatchButton.enabled = True
        self.RemovePointButton.enabled = True
        if LeftPatchNode.GetNumberOfFiducials() + NodeForRepair.GetNumberOfFiducials() > 17:
          print ""
          print "Too many points in Unrepaired and Patch node union"
          print " Can have 17 points maximum"
          print " Remove some points from the patch or unrepaired node"
          self.ApplyPatchButton.enabled = False
        
        
    if NodeForRepair == RightBaseNode:
      self.AutoGenPatchButton.enabled = True
      if RightPatchNode != None and RightPatchNode.GetName() == RightBaseNode.GetName() + "_Patch":
        self.ApplyPatchButton.enabled = True
        self.AddPointToPatchButton.enabled = True
        self.RemovePointButton.enabled = True
        if RightPatchNode.GetNumberOfFiducials() + NodeForRepair.GetNumberOfFiducials() > 17:
          print ""
          print "Too many points in Unrepaired and Patch node union"
          print " Can have 17 points maximum"
          print " Remove some points from the patch or unrepaired node"
          self.ApplyPatchButton.enabled = False
    
    # If we have a left and right side, we can combine them in a RepairedNode
    if LeftBaseNode != None and RightBaseNode != None:
      self.RepairNodeButton.enabled = True
    else:
      self.RepairNodeButton.enabled = False
    
    self.UpdateColors()
      
  def OnSelectedSideChange(self):
    NodeForRepair = self.RepairSideSelector.currentNode()
    if NodeForRepair != None:
      self.AutoGenPatchButton.enabled = True
      self.IdentifyOutliersButton.enabled = True
    else:
      self.AutoGenPatchButton.enabled = False
      self.IdentifyOutliersButton.enabled = False
    self.UpdateColors()
      
      
  def OnCategorizeLeftRight(self):
    
    # Check input data
    if self.AnatomySelector.currentNode() == None:
      print "Error - require MarkupsNode to segment into left and right sides"
      return
    
    # Following logic check may be un-necessary as it is called when Widget.AnatomySelector.currentNode() changes
    if self.logic == None:
      self.logic = PreProcessLandmarksLogic(self.AnatomySelector.currentNode(), self.kmWindowSize, self.PolynomialDegree, self.SpecificitySlider.value)
      self.SpecificitySlider.enabled = True
    
    # Get Left and Right side markups nodes
    (LeftNode, RightNode) = self.logic.PatientModel.ClassifyLeftRight()
    
    # Update scene
    PriorLeft = slicer.util.getNode(LeftNode.GetName())
    if PriorLeft != None:
      slicer.mrmlScene.RemoveNode(PriorLeft)
    
    PriorRight = slicer.util.getNode(RightNode.GetName())
    if PriorRight != None:
      slicer.mrmlScene.RemoveNode(PriorRight)    
    
    # Add new nodes to scene
    slicer.mrmlScene.AddNode(LeftNode)
    slicer.mrmlScene.AddNode(RightNode)
    
    # Update user interface
    self.LeftSideSelector.setCurrentNode(LeftNode)
    self.RightSideSelector.setCurrentNode(RightNode)
    self.UpdateColors()
    
  def OnAutoGenPatchButton(self):
    
    # Probably should necessitate continuity of data on the user's end
    #if self.logic == None:
    #  self.logic = PreProcessLandmarksLogic(self.AnatomySelector.currentNode(), self.PolynomialDegree)
    #  (LS, RS) = self.logic.CategorizeLeftRight(4)
      
    NodeToPatch = self.RepairSideSelector.currentNode()
      
    if NodeToPatch.GetName().__contains__("Left"):
      
      # Get patch
      PatchNode = self.logic.PatientModel.LeftSide.PredictAndImputeOmissions()
      
      # Update scene
      PriorNode = slicer.util.getNode(PatchNode.GetName())
      if PriorNode != None:
        slicer.mrmlScene.RemoveNode(PriorNode)
      slicer.mrmlScene.AddNode(PatchNode)
      
      # Update self class variables
      self.LeftPatchNode = PatchNode
      self.LeftPatchSelector.setCurrentNode(PatchNode)
      
    elif NodeToPatch.GetName().__contains__("Right"):
      
      # Get patch
      PatchNode = self.logic.PatientModel.RightSide.PredictAndImputeOmissions()
      
      # Update scene
      PriorNode = slicer.util.getNode(PatchNode.GetName())
      if PriorNode != None:
        slicer.mrmlScene.RemoveNode(PriorNode)
      slicer.mrmlScene.AddNode(PatchNode)
      
      # Update self class variables
      self.RightPatchNode = PatchNode
      self.RightPatchSelector.setCurrentNode(PatchNode)
    
    # Update side-independent UI 
    self.UpdateColors()
    
  def OnApplyPatchButton(self):
    
    # Get node selectors' contents, even if they are None
    NodeToPatch = self.RepairSideSelector.currentNode()
    #self.StackSaveNode(NodeToPatch)
    LeftNode = self.LeftSideSelector.currentNode()
    LeftPatchNode = self.LeftPatchSelector.currentNode()
    RightNode = self.RightSideSelector.currentNode()
    RightPatchNode = self.RightPatchSelector.currentNode()
    
    if NodeToPatch == LeftNode and LeftPatchNode != None:
      
      # Get side to patch 
      LeftSide = self.logic.PatientModel.LeftSide
      
      # Get side with patch applied
      PatchedLeftNode = LeftSide.CombineNodeWithPatch(LeftSide.MarkupsNode, LeftSide.PatchNode)
      PatchedLeftNode.SetName(LeftNode.GetName())
      
      # SpineSide.CombineNodeWithPatch does not currently update self.MarkupsNode as the function is used iteratively in searching for best subpatches
      # TODO: Fix ^
      self.logic.PatientModel.LeftSide.MarkupsNode.Copy(PatchedLeftNode)
      
      # Update scene
      slicer.mrmlScene.RemoveNode(LeftNode)
      slicer.mrmlScene.AddNode(PatchedLeftNode)
      
      # Update UI
      self.RepairSideSelector.setCurrentNode(PatchedLeftNode)
      self.LeftSideSelector.setCurrentNode(PatchedLeftNode)
    
    if NodeToPatch == RightNode and RightPatchNode != None:
      
      # Get side to patch
      RightSide = self.logic.PatientModel.RightSide
      
      # Get side with patch applied
      PatchedRightNode = RightSide.CombineNodeWithPatch(RightSide.MarkupsNode, RightSide.PatchNode)
      PatchedRightNode.SetName(RightNode.GetName())
      
      # Update scene
      slicer.mrmlScene.RemoveNode(RightNode)
      slicer.mrmlScene.AddNode(PatchedRightNode)
      
      # SpineSide.CombineNodeWithPatch does not currently update self.MarkupsNode as the function is used iteratively in searching for best subpatches
      # TODO: Fix ^
      self.logic.PatientModel.RightSide.MarkupsNode.Copy(PatchedRightNode)
      
      # Update UI
      self.RepairSideSelector.setCurrentNode(PatchedRightNode)
      self.RightSideSelector.setCurrentNode(PatchedRightNode)

    # Disable the button, we used the patch and need a new one
    self.ApplyPatchButton.enabled = False
      
    self.UpdateColors()
  
  def OnkmWindowSliderChanged(self):
    NewKMeansWindowSize = self.kmWindowSizeSlider.value
    self.kmWindowTextDisplay.setText(str(NewKMeansWindowSize))
    if self.logic != None:
      self.logic.PatientModel.kmWindowSize = NewKMeansWindowSize
    return
  
  def OnPolyFitDegreeSliderChanged(self):
    NewPolyFitDegree = self.PolyFitDegreeSlider.value
    self.PolyFitDegreeSliderTextDisplay.setText(str(NewPolyFitDegree))
    if self.logic != None:
      self.logic.PatientModel.PolyFitDegree = NewPolyFitDegree
      self.logic.PatientModel.RightSide.PolyFitDegree = NewPolyFitDegree
      self.logic.PatientModel.LeftSide.PolyFitDegree = NewPolyFitDegree
    return
  
  def OnSpecificitySliderChanged(self):
    NewOmissionSpecificity = self.SpecificitySlider.value
    self.SpecificityTextDisplay.setText(str(NewOmissionSpecificity))
    if self.logic != None:
      self.logic.PatientModel.OmissionDetectionSpecificity = NewOmissionSpecificity
      self.logic.PatientModel.RightSide.OmissionDetectionSpecificity = NewOmissionSpecificity
      self.logic.PatientModel.LeftSide.OmissionDetectionSpecificity = NewOmissionSpecificity
    #print str(self.SpecificitySlider.value)
    return
  
  def OnAddPointToPatchButton(self):
  
    # Get sub-patch to add point to 
    SubIntervalNumber = self.SubIntervalSelector.value
    
    # Get the patch node having a point added
    NodeUnderRepair = self.RepairSideSelector.currentNode()
    OldPatchName = NodeUnderRepair.GetName() + "_Patch"
    OldPatchNode = slicer.util.getNode(OldPatchName)
    
    if NodeUnderRepair.GetNumberOfFiducials() + OldPatchNode.GetNumberOfFiducials() >= 17:
      print ""
      print "Cannot add any more points to patch"
      print " Patch and unrepaired node together would have over 17 points"
      print " Remove some points from the patch or repair node first"
      print " Adding no points"
      return
    
    if OldPatchName.__contains__("Left"):
      
      # Get patch with point added
      SideUnderRepair = self.logic.PatientModel.LeftSide
      NewPatchNode = SideUnderRepair.AddPointToSubPatch(SubIntervalNumber)
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldPatchNode)
      slicer.mrmlScene.AddNode(NewPatchNode)
      
      # Update UI
      self.LeftPatchSelector.setCurrentNode(NewPatchNode)
      
    else:
      
      # Get patch with point added
      SideUnderRepair = self.logic.PatientModel.RightSide
      NewPatchNode = SideUnderRepair.AddPointToSubPatch(SubIntervalNumber)
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldPatchNode)
      slicer.mrmlScene.AddNode(NewPatchNode)
      
      # Update UI
      self.RightPatchSelector.setCurrentNode(NewPatchNode)

    self.UpdateColors() 
    
    return
    
  def OnRemovePointButton(self):
    
    # Get sub-patch to add point to 
    SubIntervalNumber = self.SubIntervalSelector.value
    NodeUnderRepair = self.RepairSideSelector.currentNode()
    OldPatchName = NodeUnderRepair.GetName() + "_Patch"
    OldPatchNode = slicer.util.getNode(OldPatchName)
    
    if OldPatchName.__contains__("Left"):
      
      # Get patch with point added
      SideUnderRepair = self.logic.PatientModel.LeftSide
      NewPatchNode = SideUnderRepair.RemovePointFromSubPatch(SubIntervalNumber)
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldPatchNode)
      slicer.mrmlScene.AddNode(NewPatchNode)
      
      # Update UI
      self.LeftPatchSelector.setCurrentNode(NewPatchNode)
      
    else:
    
      # Get patch with point added
      SideUnderRepair = self.logic.PatientModel.RightSide
      NewPatchNode = SideUnderRepair.RemovePointFromSubPatch(SubIntervalNumber)
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldPatchNode)
      slicer.mrmlScene.AddNode(NewPatchNode)
      
      # Update UI
      self.RightPatchSelector.setCurrentNode(NewPatchNode)
     
    self.UpdateColors() 
    
    return
    
  def OnIdentifyOutliersButton(self):
    
    # Get node to be searched for outliers
    OldWorkingNode = self.RepairSideSelector.currentNode()
    
    if OldWorkingNode.GetName().__contains__("Left"):
      # Get node with suspected outliers highlighted (unselected)
      NewWorkingNode = self.logic.PatientModel.LeftSide.IdentifyOutliers()
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldWorkingNode)
      #OldWorkingNode.Copy(NewWorkingNode)
      slicer.mrmlScene.AddNode(NewWorkingNode)
      
      # Update UI
      self.RepairSideSelector.setCurrentNode(NewWorkingNode)
      self.LeftSideSelector.setCurrentNode(NewWorkingNode)
      
    else:
      # Get node with suspected outliers highlighted (unselected)
      NewWorkingNode = self.logic.PatientModel.RightSide.IdentifyOutliers()
      
      # Update scene
      slicer.mrmlScene.RemoveNode(OldWorkingNode)
      #OldWorkingNode.Copy(NewWorkingNode)
      slicer.mrmlScene.AddNode(NewWorkingNode)
      
      # Update UI
      self.RepairSideSelector.setCurrentNode(NewWorkingNode)
      self.RightSideSelector.setCurrentNode(NewWorkingNode)

    self.UpdateColors()
    
    return
    
  def OnRemoveOutliersButton(self):
    
    return
    
  def OnRepairButtonClicked(self):  # Combines self.LeftSideSelector and self.RightSideSelector nodes into a node stored in self.RepairedNodeStorage
  
    # Get left and right combined into repaired node
    RepairedNode = self.logic.PatientModel.CombineRepairedSides()
    
    # Impose naming convention required by ModelToPatientRegistration extension
    RepairedNode = self.logic.ApplyNamingConvention(RepairedNode)
    
    # Update mrmlScene
    OldRepairedNode = self.RepairedNodeStorage.currentNode()
    if OldRepairedNode != None:
      slicer.mrmlScene.RemoveNode(OldRepairedNode)
    slicer.mrmlScene.AddNode(RepairedNode)
    
    # Update UI
    self.RepairedNodeStorage.setCurrentNode(RepairedNode)
    self.UpdateColors()
    
    return
    
  def StackSaveNode(self, LatestNode):
    # Create node into which LatestNode can be copied
    NodeTypeStr = LatestNode.GetClassName()
    InstantiationStr = "StorageNode = slicer." + NodeTypeStr + "()"
    exec InstantiationStr
    StorageNode.Copy(LatestNode)
    #StorageNode = LatestNode
    #slicer.mrmlScene.RemoveNode(StorageNode)
    
    # Add this node to self.NodeStack
    if len(self.NodeStack) >= self.MaxUndos:
      self.NodeStack.__delitem__(0)
    self.NodeStack.append(StorageNode)
    """
    else:
      self.NodeStack.
      for i, Node in enumerate(self.NodeStack[1:], start=1):
        self.NodeStack[i-1] = Node
      self.NodeStack[-1] = StorageNode
    """
  def OnUndoButtonClicked(self):
    if len(self.NodeStack) == 0:
      print "Nothing to undo"
      return
      
    LastNode = self.NodeStack.pop()
    
    RollbackNode = slicer.util.getNode(LastNode.GetName())
    if RollbackNode == None:
      print "Error - Could not find last node from stack"
      return
    if LastNode.GetClassName() != RollbackNode.GetClassName():
      print "Error - Undo stack top and Rollback node types don't match"
      return
    else:
      #RollbackNode = LastNode
      slicer.mrmlScene.RemoveNode(RollbackNode)
      slicer.mrmlScene.AddNode(LastNode)
      RollbackNode.Copy(LastNode)
      slicer.mrmlScene.AddNode(RollbackNode)
      slicer.mrmlScene.RemoveNode(LastNode)
    
    self.UpdateColors()
    
  def UpdateColors(self):
    # Input set
    MainNode = self.AnatomySelector.currentNode()
    if MainNode != None:
      MainNode.GetDisplayNode().SetSelectedColor(1,0,1)
    
    # Both unrepaired sides
    (LeftUnrep, RightUnrep) = (self.LeftSideSelector.currentNode(), self.RightSideSelector.currentNode())
    if LeftUnrep != None:
      LeftUnrep.GetDisplayNode().SetSelectedColor(0,0,1)
    if RightUnrep != None:
      RightUnrep.GetDisplayNode().SetSelectedColor(0,0,1)
      
    # Node selected for patching
    NodeToPatch = self.RepairSideSelector.currentNode()
    if NodeToPatch != None:
      NodeToPatch.GetDisplayNode().SetSelectedColor(1,0,0)
      
    # Selected patches
    (LeftPatchNode, RightPatchNode) = (self.LeftPatchSelector.currentNode(), self.RightPatchSelector.currentNode())
    if LeftPatchNode != None:
      LeftPatchNode.GetDisplayNode().SetSelectedColor(1,1,0)
    if RightPatchNode != None:
      RightPatchNode.GetDisplayNode().SetSelectedColor(1,1,0)    
    
    # Repaired node
    RepairedNode = self.RepairedNodeStorage.currentNode()
    if RepairedNode != None:
      RepairedNode.GetDisplayNode().SetSelectedColor(0,1,0)
    
  def OnReloadButton(self):
    self.logic = None
    #slicer.mrmlScene.Clear(0)
    self.SpecificitySlider.enabled = False
    print str(slicer.moduleNames.PreProcessLandmarks) + " reloaded"
    slicer.util.reloadScriptedModule(slicer.moduleNames.PreProcessLandmarks)
    
#
# PreProcessLandmarksLogic
#

class PreProcessLandmarksLogic(ScriptedLoadableModuleLogic):
  class PatientTransverseProcesses:
    def __init__(self, parent, Node, kmWindowSize, PolyFitDegree, OmissionDetectionSpecificity):
      self.ParentLogic = parent
    
      self.MarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      self.MarkupsNode.Copy(Node)
      self.MarkupsNode.SetAndObserveDisplayNodeID(Node.GetDisplayNodeID())
      self.LabelsCoords = [(self.MarkupsNode.GetNthFiducialLabel(i), self.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      
      self.kmWindowSize = kmWindowSize
      self.PolyFitDegree = PolyFitDegree
      self.OmissionDetectionSpecificity = OmissionDetectionSpecificity
      
      
    def SortPointsVertically(self):
      self.MarkupsNode.RemoveAllMarkups()
      self.LabelsCoords = sorted(self.LabelsCoords, key=lambda Tup: -1*Tup[1][2])
      for i, Markup in enumerate(self.LabelsCoords):
        self.MarkupsNode.AddFiducialFromArray(Markup[1])
        self.MarkupsNode.SetNthFiducialLabel(i, Markup[0])
      
    def ClassifyLeftRight(self):  

      self.SortPointsVertically()
      SortedPointsLeftVotes = len(self.LabelsCoords) * [0]
      SortedPointsRightVotes = len(self.LabelsCoords) * [0]
      for i in range(0, len(self.LabelsCoords)-self.kmWindowSize+1):
        MarkupsWindow = self.LabelsCoords[i:i+self.kmWindowSize]
        #print MarkupsWindow
        NormalizedWindow = self.AnisotrpoicNormalization(MarkupsWindow)
        (KmLabels, KmCentroids) = self.KMeans(NormalizedWindow, 2)

        # If KmLabel == 0 indicates a left-side point
        if KmCentroids[0][0] < KmCentroids[1][0]:
          for j, Label in enumerate(KmLabels):
            if Label == 0:
              SortedPointsLeftVotes[i + j] = SortedPointsLeftVotes[i + j] + 1
            else:
              SortedPointsRightVotes[i + j] = SortedPointsRightVotes[i + j] + 1
        else: # If KmLabel == 0 indicates a right-side point
          for j, Label in enumerate(KmLabels):
            #print i, j
            if Label == 0:
              SortedPointsRightVotes[i + j] = SortedPointsRightVotes[i + j] + 1
            else:
              SortedPointsLeftVotes[i + j] = SortedPointsLeftVotes[i + j] + 1
      
      NewLeftSide = slicer.vtkMRMLMarkupsFiducialNode()
      NewLeftSide.SetName(self.MarkupsNode.GetName() + "_Left")
      
      NewRightSide = slicer.vtkMRMLMarkupsFiducialNode()
      NewRightSide.SetName(self.MarkupsNode.GetName() + "_Right")
      
      for i, UnclassifiedPoint in enumerate(self.LabelsCoords):
        #OriginalPointIndex = self.LabelsCoords.index(UnclassifiedPoint)
        #OriginalPoint = self.LabelsCoords[OriginalPointIndex]
        if SortedPointsLeftVotes[i] > SortedPointsRightVotes[i]:
          NewLeftSide.AddFiducialFromArray(UnclassifiedPoint[1])
          NewLeftSide.SetNthFiducialLabel(NewLeftSide.GetNumberOfFiducials()-1, UnclassifiedPoint[0] + '_Left')
        else:
          NewRightSide.AddFiducialFromArray(UnclassifiedPoint[1])
          NewRightSide.SetNthFiducialLabel(NewRightSide.GetNumberOfFiducials()-1, UnclassifiedPoint[0] + '_Right')

      self.LeftSide = self.ParentLogic.SpineSide(self, NewLeftSide, self.PolyFitDegree, self.OmissionDetectionSpecificity)
      self.RightSide = self.ParentLogic.SpineSide(self, NewRightSide, self.PolyFitDegree, self.OmissionDetectionSpecificity)
      
      return (self.LeftSide.MarkupsNode, self.RightSide.MarkupsNode)
      
    def AnisotrpoicNormalization(self, PointSet):
      # Start by finding top and bottom points of spine for verticle normalization
      SetRight = -10000
      SetLeft = 10000
      SetFront = -10000
      SetBack = 10000
      SetTop = -10000
      SetBottom = 10000

      for Point in PointSet:
        Coords  = Point[1]
        if Coords[0] < SetLeft:
          SetLeft = Coords[0]
        if Coords[0] > SetRight:
          SetRight = Coords[0]
        if Coords[1] < SetBack:
          SetBack = Coords[1]
        if Coords[1] > SetFront:
          SetFront = Coords[1]
        if Coords[2] > SetTop:
          SetTop = Coords[2]
        if Coords[2] < SetBottom:
          SetBottom = Coords[2]
        
      SetHeight = SetTop - SetBottom
      SetWidth = SetRight - SetLeft
      SetDepth = SetFront - SetBack
      
      # (Re) initialize normalized point list
      NormalizedPoints = len(PointSet) * [0]
      
      # Normalize S-I dimension to R-L scale
      for i, Point in enumerate(PointSet):
        Coords = Point[1]
        NormalizedPoints[i] = np.array([Coords[0], (Coords[1]) * (SetWidth / (SetDepth * 3.0)), (Coords[2]) * (SetWidth / (SetHeight * 2.0))])
      return NormalizedPoints
    
    def ShouldStopKMeans(self, oldCentroids, Centroids, iterations):
      MaxIterations = 500
      StoppingDelta = 0.05
      #print oldCentroids, Centroids
      if iterations > MaxIterations: return True
      if iterations == 0:
        return False
      for C in range(len(oldCentroids)):
        #print oldCentroids[C], Centroids[C]
        if np.linalg.norm(oldCentroids[C] - Centroids[C]) < StoppingDelta:
          return True
      return False

    def GetKMeansLabels(self, DataSet, Centroids, Labels):
      for i, Coords in enumerate(DataSet):
        #PointCentroidDistance = np.linalg.norm(Coords - Centroids[0])
        minDist = 1000000
        Labels[i] = 0
        for j, CentroidVector in enumerate(Centroids):
          #print Coords, CentroidVector#, np.linalg.norm(Coords - CentroidVector)
          PointCentroidDistance = np.linalg.norm(Coords - CentroidVector)
          if PointCentroidDistance < minDist:
            minDist = PointCentroidDistance
            Labels[i] = j
      return Labels

    def GetCentroids(self, DataSet, Labels, k):
      Centroids = []
      #print Labels
      for C in range(k):    # For each centroid
        Centroid = np.random.uniform() * np.ones(len(DataSet[0]))    # Each centroid with as many dimensions as the data
        
        for i, Coords in enumerate(DataSet): # Take each data point contributing to the centroid into consideration
          if Labels[i] == C:                    # if it belongs to the current centroid
            for dim in range(len(Centroid)):
              Centroid[dim] += Coords[dim]
              
        for dim in range(len(Centroid)):
          Centroid[dim] = Centroid[dim] / np.count_nonzero(Labels == C)
        Centroids.append(Centroid)
      return Centroids

    def KMeans(self, DataSet, k=2):   # Expects DataSet as list of (Label, np.array[R,A,S]) tuples
      DataSetLabels = np.zeros(len(DataSet))
      for i in range(len(DataSetLabels)):
        DataSetLabels[i] = int(round(np.random.uniform()))
      # Initialize centroids
      numFeatures = len(DataSet[0])
      Centroids = k * [0]
      # Try initializing one centroid on each side
      Centroids[0] = np.array([min([i[0] for i in DataSet]),0,0])
      Centroids[1] = np.array([max([i[0] for i in DataSet]),0,0])
      
    # Initialize book keeping variables
      iterations = 0
      oldCentroids = k * [0]
      for Cent in range(k):
        oldCentroids[Cent] = np.array(numFeatures * [np.random.uniform()])

      # Run the k-means algorithm
      while not self.ShouldStopKMeans(oldCentroids, Centroids, iterations):
        oldCentroids = Centroids
        iterations += 1
        
        DataSetLabels = self.GetKMeansLabels(DataSet, Centroids, DataSetLabels)
        Centroids = self.GetCentroids(DataSet, DataSetLabels, k)
        #print Centroids
      return (DataSetLabels, Centroids)

    def CombineRepairedSides(self):
      
      # Retrieve and sort LabelCoords of Left and Right sides
      LeftLabelCoords = [(self.LeftSide.MarkupsNode.GetNthFiducialLabel(i), self.LeftSide.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.LeftSide.MarkupsNode.GetNumberOfFiducials())]
      SortedLeftLabelCoords = sorted(LeftLabelCoords, key=lambda Tup: -1*Tup[1][2])
      RightLabelsCoords = [(self.RightSide.MarkupsNode.GetNthFiducialLabel(j), self.RightSide.MarkupsNode.GetMarkupPointVector(j,0)) for j in range(self.RightSide.MarkupsNode.GetNumberOfFiducials())]
      SortedRightLabelCoords = sorted(RightLabelsCoords, key=lambda Tup: -1*Tup[1][2])
      
      # Consolidate Left,Right LabelCoords - in that order
      # Assumes that left and right sides start at the same vertebra and each point corresponds to one on the other side
      CombinedLabelCoords = []
      for i, (LeftPoint, RightPoint) in enumerate(zip(LeftLabelCoords, RightLabelsCoords)):
        CombinedLabelCoords.append(LeftPoint)
        CombinedLabelCoords.append(RightPoint)
      #CombinedLabelCoords = LeftLabelCoords + RightLabelsCoords
      
      # Create node to populate with markups using CombinedLabelCoordsLabelCoords
      RepairedNode = slicer.vtkMRMLMarkupsFiducialNode()
      RepairedNode.SetName(self.MarkupsNode.GetName() + "_Repaired")
      
      # Populate RepairedNode with markups points
      for i, LabelCoord in enumerate(CombinedLabelCoords):
        RepairedNode.AddFiducialFromArray(LabelCoord[1])
        RepairedNode.SetNthFiducialLabel(i, LabelCoord[0])

      return RepairedNode
      
  class SpineSide:
    def __init__(self, parent, Node, PolyFitDegree, OmissionDetectionSpecificity):
      self.ParentPatient = parent
      self.PointsPerPolynomialCurve = 500
      
      self.MarkupsNode = slicer.vtkMRMLMarkupsFiducialNode()
      self.MarkupsNode.Copy(Node)
      
      self.PatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      self.PatchNode.SetName(parent.MarkupsNode.GetName() + "_Patch")
      
      # Each element corresponds to an interval of self.MarkupsNode identified as having omissions,
      self.SubPatchesPointCounts = []      # the value represents the number of omissions to be estimated in self.PatchNode
      
      self.PolyFitDegree = PolyFitDegree
      self.OmissionDetectionSpecificity = OmissionDetectionSpecificity
      
      #self.OutlierIdentificationVotes = np.zeros(self.MarkupsNode.GetNumberOfFiducials())
      #self.ImputationErrorImprovementThreshold = 0.05
     
    def OrderPointsSuperiorToInferior(self, Node):
      #print  "Sorting node " + Node.GetName() + " landmarks"
      LabelsCoords = [(Node.GetNthFiducialLabel(i), Node.GetMarkupPointVector(i,0)) for i in range(Node.GetNumberOfFiducials())]
      Node.RemoveAllMarkups()
      LabelsCoords = sorted(LabelsCoords, key=lambda Tup: -1*Tup[1][2])
      for i, (Label, Coord) in enumerate(LabelsCoords):
        Node.AddFiducialFromArray(Coord)
        Node.SetNthFiducialLabel(i, Label)
     
    def CombineNodeWithPatch(self, BaseNode, PatchNode):
      # Returns a node containing the original markups points, and those from the patch - used to check patch fit improvements
      
      TentativeNode = slicer.vtkMRMLMarkupsFiducialNode()
      OriginalLabelPoints = [(BaseNode.GetNthFiducialLabel(i), BaseNode.GetMarkupPointVector(i,0)) for i in range(BaseNode.GetNumberOfFiducials())]
      PatchLabelPoints = [(PatchNode.GetNthFiducialLabel(i), PatchNode.GetMarkupPointVector(i,0)) for i in range(PatchNode.GetNumberOfFiducials())]
      AllLabelPoints = OriginalLabelPoints + PatchLabelPoints
      
      # Also sort points
      AllLabelPoints = sorted(AllLabelPoints, key=lambda LabelPoint: -1*LabelPoint[1][2])
      
      for i, (Label, Coords) in enumerate(AllLabelPoints):
        TentativeNode.AddFiducialFromArray(Coords)
        TentativeNode.SetNthFiducialLabel(i, Label)
        
      return TentativeNode
     
    def CoordsPolyFit(self, Node):
      #Coords = [self.MarkupsNode.GetMarkupPointVector(i,0) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      #print Coords
      #LabelsCoords = [(self.MarkupsNode.GetNthFiducialLabel(i), self.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      LabelsCoords = [(Node.GetNthFiducialLabel(i), Node.GetMarkupPointVector(i,0)) for i in range(Node.GetNumberOfFiducials())]
      R = [LabelsCoords[i][1][0] for i in range(len(LabelsCoords))]
      A = [LabelsCoords[i][1][1] for i in range(len(LabelsCoords))]
      S = [LabelsCoords[i][1][2] for i in range(len(LabelsCoords))]
      sSpace = np.linspace(S[0], S[-1], self.PointsPerPolynomialCurve)

      # The degrees of these polynomials should not be hard coded
      S_R_FitCoefs = np.polyfit(S, R, self.PolyFitDegree)
      S_A_FitCoefs = np.polyfit(S, A, self.PolyFitDegree)
      
      SrPolynomial = np.poly1d(S_R_FitCoefs)
      SaPolynomial = np.poly1d(S_A_FitCoefs)
  
      return (SrPolynomial, SaPolynomial)
      
    def GetCurvewiseInterpointDistances(self, Node, SrPolynomial, SaPolynomial):
      #Coords = [self.MarkupsNode.GetMarkupPointVector(i,0) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      LabelsCoords = [(Node.GetNthFiducialLabel(i), Node.GetMarkupPointVector(i,0)) for i in range(Node.GetNumberOfFiducials())]
      R = [LabelsCoords[i][1][0] for i in range(len(LabelsCoords))]
      A = [LabelsCoords[i][1][1] for i in range(len(LabelsCoords))]
      S = [LabelsCoords[i][1][2] for i in range(len(LabelsCoords))]
      sSpace = np.linspace(S[0], S[-1], self.PointsPerPolynomialCurve)
     
      # Distances along the polynomial to each landmark in the given dimension
      PointDistances = []
      #priorS = sSpace[0]
      sIndex = 1
      
      # Find points in sSpace corresponding to landmarks in both R and A dimension
      # Using "sSpace" for each polynomial's independent variable allows the use of the S coordinate to recognize when we've reached a point along the polynomials
      CurveDistance = 0
      for i, Landmark in enumerate(zip(R[:-1],A[:-1])):
        CurrentIntervalLength = 0
        while sIndex < self.PointsPerPolynomialCurve and sSpace[sIndex] > S[i+1]:
          PriorS = sSpace[sIndex-1]
          CurrentS = sSpace[sIndex]
          PriorR = SrPolynomial(PriorS)
          CurrentR = SrPolynomial(CurrentS)
          PriorA = SaPolynomial(PriorS)
          CurrentA = SaPolynomial(CurrentS)
          CurveIncrementDistance = np.sqrt(((CurrentS - PriorS)**2) + ((CurrentR - PriorR)**2) + ((CurrentA - PriorA)**2)) 
          CurrentIntervalLength += CurveIncrementDistance
          CurveDistance += CurveIncrementDistance
          sIndex += 1
        PointDistances.append((CurveDistance, CurrentIntervalLength))
          
      return PointDistances
      
    def FrequencyPolyFit(self, IntervalData):
      IntervalDistances = [IntervalData[i][1] for i in range(len(IntervalData))]
      
      # Try fitting polynomial to relationship between total distance travelled down the spine to the inter-landmark distances
      CumulativeDistance = 0
      CumulativeCurveDistances = []
      for i, Interval in enumerate(IntervalDistances):
        CumulativeCurveDistances.append(CumulativeDistance + (Interval/2.0))
        CumulativeDistance += Interval
      
      CurveSpace = np.linspace(0, CumulativeDistance, self.PointsPerPolynomialCurve)
      FrequencyCoeffs = np.polyfit(CumulativeCurveDistances, IntervalDistances, self.PolyFitDegree)
      FrequencyPolynomial = np.poly1d(FrequencyCoeffs)
      
      return FrequencyPolynomial
      
    def GetPolyfitErrors(self, Data, Polynomial):
      FitErrors = []
      for Datum in Data:
        PolyPrediction = Polynomial(Datum[0])
        PredictionError = Datum[1] - PolyPrediction
        FitErrors.append(PredictionError)
      return FitErrors
  
    def GetPolyfitRMS(self, FitErrors):
      if len(FitErrors) == 0:
        print "Computing RMS error of empty error set - returning 0"
        return 0
      else:
        SumSquaredError = sum([(Error)**2 for Error in FitErrors])
        RMS = np.sqrt(SumSquaredError)
        #print RMS
        return RMS
       
    def GetNonparametricSkewness(self, Data):   # 
      Mean = np.mean(Data)
      Median = np.median(Data)
      StdDev = np.std(Data)
      Skewness = (Mean - Median) / StdDev
      return Skewness
       
    def IdentifyOmissionsFromLocalDistances(self, IntervalIndex, IntervalData):  # Not currently used, other Identifiers work better so far
      IntervalLength = IntervalData[IntervalIndex][1]
      
      # The first interval is a boundary, having only one neighboring interval
      if IntervalIndex == 0:
        NextIntervalLength = IntervalData[IntervalIndex+1][1]
        if NextIntervalLength < 0.75*IntervalLength:
          return True
        else:
          return False
          
      # The last interval is also a boundary condiditon
      if IntervalIndex == len(IntervalData)-1:
        PriorIntervalLength = IntervalData[IntervalIndex-1][1]
        if PriorIntervalLength < 0.75*IntervalLength:
          return True
        else:
          return False
          
      #else:
      NextIntervalLength = IntervalData[IntervalIndex+1][1]
      PriorIntervalLength = IntervalData[IntervalIndex-1][1]
      if NextIntervalLength < 0.75*IntervalLength or PriorIntervalLength < 0.75*IntervalLength:
        return True
      else:
        return False
   
    def IdentifyOmissionsFromIntervalLengthFit(self, IntervalIndex, FitErrors):
      AbsErrors = [abs(E) for E in FitErrors]
      IntervalError = FitErrors[IntervalIndex]
      MeanAbsError = np.mean(AbsErrors)
      ErrorStdDev = np.std(FitErrors)
      
      # Try skewness instead of standard deviation - 
      IntervalLengthSkewness = self.GetNonparametricSkewness(FitErrors)
      
      # The first interval is a boundary, having only one neighboring interval
      if IntervalIndex == 0:
        if abs(IntervalError) > MeanAbsError + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
          
      # The last interval is also a boundary condiditon
      elif IntervalIndex == len(FitErrors)-1:
        if abs(IntervalError) > MeanAbsError + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
          
      else:
        PriorIntervalError = FitErrors[IntervalIndex-1]
        NextIntervalError = FitErrors[IntervalIndex+1]
        if (np.sign(IntervalError) != np.sign(PriorIntervalError) and abs(PriorIntervalError)!=max(AbsErrors)) and (np.sign(IntervalError) != np.sign(NextIntervalError) and abs(NextIntervalError)!=max(AbsErrors)) and (abs(IntervalError) > MeanAbsError):
          return True
        if abs(IntervalError) > MeanAbsError + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
    
    def IdentifyOmissionsFromLocalIntervalLengths(self, IntervalIndex, IntervalData):
      IntervalLengths = [IntervalData[i][1] for i in range(len(IntervalData))]
      IntervalLengthMean = np.mean(IntervalLengths)
      IntervalLengthStd = np.std(IntervalLengths)
      
      # Try skewness instead of standard deviation - 
      IntervalLengthSkewness = self.GetNonparametricSkewness(IntervalLengths)
      
      CurrentIntervalLength = IntervalLengths[IntervalIndex]
      
      # Boundary condition checks leave the possibility that two intervals with omissions neighbor
      
      # Superior-most boundary condition - only one neighbor, below
      if IntervalIndex == 0:
        NextIntervalLength = IntervalLengths[IntervalIndex+1]
        if CurrentIntervalLength > NextIntervalLength + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
        
      # Inferior-most boundary condition - only one neighbor, above
      if IntervalIndex == len(IntervalData)-1:
        PriorIntervalLength = IntervalLengths[IntervalIndex-1]
        if CurrentIntervalLength > PriorIntervalLength + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
      
      else:         # Interval is somewhere in the middle, with two neighbors
        NextIntervalLength = IntervalLengths[IntervalIndex+1]
        PriorIntervalLength = IntervalLengths[IntervalIndex-1]
        if CurrentIntervalLength > NextIntervalLength + (self.OmissionDetectionSpecificity*IntervalLengthSkewness) or CurrentIntervalLength > PriorIntervalLength + (self.OmissionDetectionSpecificity*IntervalLengthSkewness):
          return True
        else:
          return False
     
    def EstimateOmissions(self, Node, IntervalIndex, IntervalData):   # Finds and adds the optimum number of points to a given interval to minimize the entire curve's frequency fit RMS
      # Returns CountEstimate - If CountEstimate == 0 and this method was called, an infinite loop may result from omissions being identified and never fixed
      LabelsCoords = [(Node.GetNthFiducialLabel(i), Node.GetMarkupPointVector(i,0)) for i in range(Node.GetNumberOfFiducials())]
      S = [LabelsCoords[i][1][2] for i in range(len(LabelsCoords))]
      
      # Initialize original measures for comparison
      (OriginalSr, OriginalSa) = self.CoordsPolyFit(Node)
      OriginalFreqPoly = self.FrequencyPolyFit(IntervalData)
      OriginalFitErrors = self.GetPolyfitErrors(IntervalData, OriginalFreqPoly)
      BestRMS = self.GetPolyfitRMS(OriginalFitErrors)
      BestIntervalFitError = abs(OriginalFitErrors[IntervalIndex])
      #print OriginalFitErrors
      
      print ""
      print " About to test imputation - ", self.MarkupsNode.GetName(), " - Interval ", IntervalIndex
      
      # Try one imputation, as boundary, to get initial tentative improvement
      CountEstimate = 1
      print "Try ", CountEstimate, " imputation"

      SubPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      """
      if self.UseLinearInterpolation:
        PointAboveInterval = Node.GetMarkupPointVector(IntervalIndex,0)
        PointBelowInterval = Node.GetMarkupPointVector(IntervalIndex,0)
        DirectionVector = [PointBelowInterval[dim]-PointAboveInterval[dim] for dim in range(3)]
        DirectionUnitVector = [DirectionVector[dim]/np.linalg.norm(DirectionVector) for dim in range(3)]
        for NewPoint in range(CountEstimate):
          SubPatchNode.AddFiducialFromArray([(NewPoint+1) * PointAboveInterval[dim] for dim in range(3))
      """
      #else:
      #TentativePatchNode.Copy(self.MarkupsNode)
      # Create SubPatchNode and add points by interpolating Node characteristic polynomials
      NewSLocations = np.linspace(S[IntervalIndex], S[IntervalIndex+1], CountEstimate+2)[1:-1]
      for x in NewSLocations: 
        SubPatchNode.AddFiducialFromArray([OriginalSr(x), OriginalSa(x), x])
        #TentativePatchNode.AddFiducialFromArray([OriginalSr(x), OriginalSa(x), x])
              
      # Create Node containing MarkupsNode combined with tentative subpatches to see how the number of imputations improves fit
      TentativeRepairNode = slicer.vtkMRMLMarkupsFiducialNode()
      TentativeRepairNode.Copy(self.CombineNodeWithPatch(Node, SubPatchNode))
      TentativeRepairNode.SetName("TentPat")
      
      # Create TentativeSide object for operations needed to reasses fit error
      TentativeSide = self.ParentPatient.ParentLogic.SpineSide(self.ParentPatient, TentativeRepairNode, self.PolyFitDegree, self.OmissionDetectionSpecificity)
      TentativeSide.OrderPointsSuperiorToInferior(TentativeSide.MarkupsNode)
      
      # Use said operations to update measures
      (TentSrPoly, TentSaPoly) = TentativeSide.CoordsPolyFit(TentativeSide.MarkupsNode)
      #TentInterpointData = TentativeSide.GetCurvewiseInterpointDistances(OriginalSr, OriginalSa)
      TentInterpointData = TentativeSide.GetCurvewiseInterpointDistances(TentativeSide.MarkupsNode, TentSrPoly, TentSaPoly)
      TentFreqPolynomial = TentativeSide.FrequencyPolyFit(TentInterpointData)
      TentFitErrors = TentativeSide.GetPolyfitErrors(TentInterpointData, TentFreqPolynomial)
      #TentFitErrors = TentativeSide.GetPolyfitErrors(TentInterpointData, OriginalFreqPoly)
      TentRMS = TentativeSide.GetPolyfitRMS(TentFitErrors)
      TentIntervalFitError = np.mean([abs(Error) for Error in TentFitErrors[IntervalIndex:IntervalIndex+CountEstimate+1]])
      
      # Try 2 or more imputations in the interval, see if the freq. fit improves
      
      while TentIntervalFitError - BestIntervalFitError < 0:    # While the latest addition improved the frequency fit over the original polynomial
        #print ""
        #print TentIntervalFitError - BestIntervalFitError
        #print TentFitErrors
        #print ""
        CountEstimate += 1            # Try adding one more point
        print "Try ", CountEstimate, " imputations"
        BestRMS = TentRMS
        BestIntervalFitError = TentIntervalFitError

        SubPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
        """
        if self.UseLinearInterpolation:
          PointAboveInterval = Node.GetMarkupPointVector(IntervalIndex,0)
          PointBelowInterval = Node.GetMarkupPointVector(IntervalIndex,0)
          DirectionVector = [PointBelowInterval[dim]-PointAboveInterval[dim] for dim in range(3)]
          DirectionUnitVector = [DirectionVector[dim]/np.linalg.norm(DirectionVector) for dim in range(3)]
          for NewPoint in range(CountEstimate):
            SubPatchNode.AddFiducialFromArray([(NewPoint+1) * PointAboveInterval[dim] for dim in range(3))
        else:
        """
        #TentativePatchNode.Copy(self.MarkupsNode)
        # Create SubPatchNode and add points by interpolating Node characteristic polynomials
        NewSLocations = np.linspace(S[IntervalIndex], S[IntervalIndex+1], CountEstimate+2)[1:-1]
        for x in NewSLocations: 
          SubPatchNode.AddFiducialFromArray([OriginalSr(x), OriginalSa(x), x])
          #TentativePatchNode.AddFiducialFromArray([OriginalSr(x), OriginalSa(x), x])
      
        TentativeRepairNode.Copy(self.CombineNodeWithPatch(TentativeRepairNode, SubPatchNode))
        TentativeRepairNode.SetName("TentRepair")
        TentativeSide = self.ParentPatient.ParentLogic.SpineSide(self.ParentPatient, TentativeRepairNode, self.PolyFitDegree, self.OmissionDetectionSpecificity)
        TentativeSide.OrderPointsSuperiorToInferior(TentativeSide.MarkupsNode)
        (TentSrPoly, TentSaPoly) = TentativeSide.CoordsPolyFit(TentativeSide.MarkupsNode)
        TentInterpointData = TentativeSide.GetCurvewiseInterpointDistances(TentativeSide.MarkupsNode, TentSrPoly, TentSaPoly)
        #TentInterpointData = TentativeSide.GetCurvewiseInterpointDistances(OriginalSr, OriginalSa)
        TentFreqPolynomial = TentativeSide.FrequencyPolyFit(TentInterpointData)
        TentFitErrors = TentativeSide.GetPolyfitErrors(TentInterpointData, TentFreqPolynomial)
        #TentFitErrors = TentativeSide.GetPolyfitErrors(TentInterpointData, OriginalFreqPoly)
        TentRMS = TentativeSide.GetPolyfitRMS(TentFitErrors)
        TentIntervalFitError = np.mean([abs(Error) for Error in TentFitErrors[IntervalIndex:IntervalIndex+CountEstimate+1]])
        
      # ASSERT TentRMS > BestRMS, because CountEstimate is one too high
      CountEstimate -= 1
      
      print " Done testing imputation - ", self.MarkupsNode.GetName(), " - Interval ", IntervalIndex
      print ""
      print "Frequency fit optimized by imputing ", CountEstimate, " points in interval ", IntervalIndex, " - ", self.MarkupsNode.GetName()
      #print CountEstimate
      
      # Create node with optimum number of imputations

      SubPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      SubPatchNode.SetName(self.MarkupsNode.GetName() + "_SubPatch-" + str(IntervalIndex))
      
      """
      if self.UseLinearInterpolation:
          PointAboveInterval = Node.GetMarkupPointVector(IntervalIndex,0)
          PointBelowInterval = Node.GetMarkupPointVector(IntervalIndex,0)
          DirectionVector = [PointBelowInterval[dim]-PointAboveInterval[dim] for dim in range(3)]
          DirectionUnitVector = [DirectionVector[dim]/np.linalg.norm(DirectionVector) for dim in range(3)]
          for NewPoint in range(CountEstimate):
            SubPatchNode.AddFiducialFromArray([(NewPoint+1) * PointAboveInterval[dim] for dim in range(3))
      else:
      """
      NewSLocations = np.linspace(S[IntervalIndex], S[IntervalIndex+1], CountEstimate+2)[1:-1]
      for i, x in enumerate(NewSLocations): 
        SubPatchNode.AddFiducialFromArray([OriginalSr(x), OriginalSa(x), x])
        NewLabel = str("SubPatch_" + str(len(self.SubPatchesPointCounts)) + "-Point_" + str(i))
        print ""
        SubPatchNode.SetNthFiducialLabel(i, NewLabel)

      #self.OrderPointsSuperiorToInferior(SubPatchNode)

      return SubPatchNode
      
    def PredictAndImputeOmissions(self):
      self.SubPatchesPointCounts = []  # Reinitialize
      # SkipIntervals only used when looping over the side a number of times, currently not implemented
      SkipIntervals = np.zeros(self.MarkupsNode.GetNumberOfFiducials()-1)   # Used to keep track of which intervals were tested for omissions - IdentifyOmissionsFromIntervalLengthFit may indicate omissions while EstimateOmissions might impute none
      
      SubPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      PatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      
      TentativeRepairNode = slicer.vtkMRMLMarkupsFiducialNode()
      TentativeRepairNode.Copy(self.MarkupsNode)
        
      (SrPolynomial, SaPolynomial) = self.CoordsPolyFit(self.MarkupsNode)
      CurvewiseInterpointData = self.GetCurvewiseInterpointDistances(self.MarkupsNode, SrPolynomial, SaPolynomial)
      CurvewiseFrequencyPolynomial = self.FrequencyPolyFit(CurvewiseInterpointData)
      FrequencyFitErrors = self.GetPolyfitErrors(CurvewiseInterpointData, CurvewiseFrequencyPolynomial)

      print ""
      print "About to look for intervals with omissions  - ", self.MarkupsNode.GetName()
      ImputingPoints = True
      #while ImputingPoints:
      ImputingPoints = False
        
      for TentativeInterval, (IntervalDatum, FitError) in enumerate(zip(CurvewiseInterpointData, FrequencyFitErrors), start=sum(self.SubPatchesPointCounts)):
        # For each interval in the TentativeRepairNode
        CorrespondingOriginalInterval = TentativeInterval - sum(self.SubPatchesPointCounts)
        if not SkipIntervals[TentativeInterval] and (self.IdentifyOmissionsFromIntervalLengthFit(CorrespondingOriginalInterval, FrequencyFitErrors) or self.IdentifyOmissionsFromLocalIntervalLengths(CorrespondingOriginalInterval, CurvewiseInterpointData)):
          print ""
          print "Interval " + str(CorrespondingOriginalInterval) + " apparently missing points - ", self.MarkupsNode.GetName()
          ImputingPoints = True
          #print CorrespondingOriginalInterval, TentativeInterval
          
          # Find optimum number of points to impute in interval i
          SubPatchNode = self.EstimateOmissions(TentativeRepairNode, CorrespondingOriginalInterval, CurvewiseInterpointData)
          
          # Keep track of interval changes
          NumImputations = SubPatchNode.GetNumberOfFiducials()
          if NumImputations > 0:
            self.SubPatchesPointCounts.append(NumImputations)
          SkipIntervals = np.insert(SkipIntervals, TentativeInterval, np.ones(NumImputations+1))      # Don't check intervals we've already identified as incomplete - infinite loop
          
          """
          # Update node tracking overall patch improvement
          TentativeRepairNode = self.CombineNodeWithPatch(TentativeRepairNode, SubPatchNode)
          
          # Update node characteristics
          (SrPolynomial, SaPolynomial) = self.CoordsPolyFit(TentativeRepairNode)
          CurvewiseInterpointData = self.GetCurvewiseInterpointDistances(TentativeRepairNode, SrPolynomial, SaPolynomial)
          CurvewiseFrequencyPolynomial = self.FrequencyPolyFit(CurvewiseInterpointData)
          FrequencyFitErrors = self.GetPolyfitErrors(CurvewiseInterpointData, CurvewiseFrequencyPolynomial)
          #break               # Start again from the top of the newly shaped (frequency fit) spine
          """
          PatchNode = self.CombineNodeWithPatch(PatchNode, SubPatchNode)

        SubPatchNode.RemoveAllMarkups()
      print "Done looking for intervals with omissions - ", self.MarkupsNode.GetName()

      self.PatchNode.Copy(PatchNode)
      self.PatchNode.SetName(self.MarkupsNode.GetName() + "_Patch")
      
      return self.PatchNode
   
    def AddPointToSubPatch(self, SubPatchIndex):    # Updates self.PatchNode to include on more point in interval corresponding to SubPatchIndex
      if SubPatchIndex > len(self.SubPatchesPointCounts):
        print "Non-existent patch interval - cannot add point"
        return self.PatchNode
      
      # Identify top and bottom points of SubPatch so they and the points they bound can be modified
      SubPatchIntervalStartIndex = sum(self.SubPatchesPointCounts[:SubPatchIndex])
      SubPatchIntervalEndIndex = sum(self.SubPatchesPointCounts[:SubPatchIndex]) + self.SubPatchesPointCounts[SubPatchIndex] - 1
      
      TopSubPatchPoint = (self.PatchNode.GetNthFiducialLabel(SubPatchIntervalStartIndex), self.PatchNode.GetMarkupPointVector(SubPatchIntervalStartIndex,0))
      BottomSubPatchPoint = (self.PatchNode.GetNthFiducialLabel(SubPatchIntervalEndIndex), self.PatchNode.GetMarkupPointVector(SubPatchIntervalEndIndex,0))
      
      # Find last point before, and first point after SubPatch's interval
      LabelsCoords = [(self.MarkupsNode.GetNthFiducialLabel(i), self.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      for i, LabelCoord in enumerate(LabelsCoords[:-1]):
        if LabelsCoords[i+1][1][2] < BottomSubPatchPoint[1][2]:
          BoundingPointBelow = LabelCoord
          break
        
      LabelsCoords.reverse()
      for i, LabelCoord in enumerate(LabelsCoords[:-1]):
        if LabelsCoords[i+1][1][2] > TopSubPatchPoint[1][2]:
          BoundingPointAbove = LabelCoord
          break
      
      OriginalPatchPoints = [(self.PatchNode.GetNthFiducialLabel(i), self.PatchNode.GetMarkupPointVector(i,0)) for i in range(self.PatchNode.GetNumberOfFiducials())]
      NewPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewPatchNode.SetName(self.PatchNode.GetName())
      
      # Copy points above interval having addition, if any
      for OriginalPoint in OriginalPatchPoints[:SubPatchIntervalStartIndex]:
        NewPatchNode.AddFiducialFromArray(OriginalPoint[1])
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, OriginalPoint[0])      
      
      
      # Populate interval with one more point than before
      (SrPolynomial, SaPolynomial) = self.CoordsPolyFit(self.MarkupsNode)
      NewSiLocations = np.linspace(BoundingPointAbove[1][2], BoundingPointBelow[1][2], SubPatchIntervalEndIndex - SubPatchIntervalStartIndex + 4)[1:-1]
      NewPointCoords = [[SrPolynomial(S), SaPolynomial(S), S] for S in NewSiLocations]
      for i, NewPoint in enumerate(NewPointCoords):
        NewPatchNode.AddFiducialFromArray(NewPoint)
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, "SubPatch-"+str(SubPatchIndex)+"_Point-"+str(i))
      
      # Copy points inferior to interval having addition, if any
      for OriginalPoint in OriginalPatchPoints[SubPatchIntervalEndIndex+1:]:
        NewPatchNode.AddFiducialFromArray(OriginalPoint[1])
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, OriginalPoint[0])
      
      # Update self attributes
      self.SubPatchesPointCounts[SubPatchIndex] += 1
      self.PatchNode.Copy(NewPatchNode)
      
      return self.PatchNode
      
    def RemovePointFromSubPatch(self, SubPatchIndex):
      if SubPatchIndex > len(self.SubPatchesPointCounts):
        print "Non-existent patch interval - no points to remove"
        return self.PatchNode
      
      # Identify top and bottom points of SubPatch so they and the points they bound can be modified
      SubPatchIntervalStartIndex = sum(self.SubPatchesPointCounts[:SubPatchIndex])
      SubPatchIntervalEndIndex = sum(self.SubPatchesPointCounts[:SubPatchIndex]) + self.SubPatchesPointCounts[SubPatchIndex] - 1
      
      TopSubPatchPoint = (self.PatchNode.GetNthFiducialLabel(SubPatchIntervalStartIndex), self.PatchNode.GetMarkupPointVector(SubPatchIntervalStartIndex,0))
      BottomSubPatchPoint = (self.PatchNode.GetNthFiducialLabel(SubPatchIntervalEndIndex), self.PatchNode.GetMarkupPointVector(SubPatchIntervalEndIndex,0))
      
      # Find last point before, and first point after SubPatch's interval
      LabelsCoords = [(self.MarkupsNode.GetNthFiducialLabel(i), self.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      for i, LabelCoord in enumerate(LabelsCoords[:-1]):
        if LabelsCoords[i+1][1][2] < BottomSubPatchPoint[1][2]:
          BoundingPointBelow = LabelCoord
          break
        
      LabelsCoords.reverse()
      for i, LabelCoord in enumerate(LabelsCoords[:-1]):
        if LabelsCoords[i+1][1][2] > TopSubPatchPoint[1][2]:
          BoundingPointAbove = LabelCoord
          break
      
      OriginalPatchPoints = [(self.PatchNode.GetNthFiducialLabel(i), self.PatchNode.GetMarkupPointVector(i,0)) for i in range(self.PatchNode.GetNumberOfFiducials())]
      NewPatchNode = slicer.vtkMRMLMarkupsFiducialNode()
      NewPatchNode.SetName(self.PatchNode.GetName())
      
      # Copy points superior to the interval having deletions, if any
      for OriginalPoint in OriginalPatchPoints[:SubPatchIntervalStartIndex]:
        NewPatchNode.AddFiducialFromArray(OriginalPoint[1])
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, OriginalPoint[0])      
      
      # Populate interval with one point fewer than before
      NewSiLocations = np.linspace(BoundingPointAbove[1][2], BoundingPointBelow[1][2], SubPatchIntervalEndIndex - SubPatchIntervalStartIndex + 2)[1:-1]
      (SrPolynomial, SaPolynomial) = self.CoordsPolyFit(self.MarkupsNode)
      NewPointCoords = [[SrPolynomial(S), SaPolynomial(S), S] for S in NewSiLocations]
      for i, NewPoint in enumerate(NewPointCoords):
        NewPatchNode.AddFiducialFromArray(NewPoint)
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, "SubPatch-"+str(SubPatchIndex)+"_Point-"+str(i))
      
      # Keep track of sub patches with self.SubPatchesPointCounts
      self.SubPatchesPointCounts[SubPatchIndex] -= 1
          
      # Copy points inferior to the interval having deletions, if any
      for OriginalPoint in OriginalPatchPoints[SubPatchIntervalEndIndex+1:]:
        NewPatchNode.AddFiducialFromArray(OriginalPoint[1])
        NewPatchNode.SetNthFiducialLabel(NewPatchNode.GetNumberOfFiducials()-1, OriginalPoint[0])
      
      if self.SubPatchesPointCounts[SubPatchIndex] == 0:
        self.SubPatchesPointCounts.__delitem__(SubPatchIndex)
        PatchPointsBeforeIndexShift = sum(self.SubPatchesPointCounts[:SubPatchIndex])
        # Adjust the subsequent Subpatch indices
        for j, SubPatch in enumerate(self.SubPatchesPointCounts[SubPatchIndex:]):
          for i, SubPatchPointNumber in enumerate(range(SubPatch)):
            PointsAlreadyRelabelled = sum(self.SubPatchesPointCounts[SubPatchIndex:j]) + i
            NewPatchNode.SetNthFiducialLabel(PatchPointsBeforeIndexShift+SubPatchPointNumber+PointsAlreadyRelabelled, "SubPatch-"+str(SubPatchIndex+j)+"_Point-"+str(i))
      
      self.PatchNode.Copy(NewPatchNode)
      
      return self.PatchNode
   
    def IdentifyOutliers(self):     # self.MarkupsNode remains the same EXCEPT that points identified as outliers are UNSELECTED
      #from np import math
      LabelsCoords = [(self.MarkupsNode.GetNthFiducialLabel(i), self.MarkupsNode.GetMarkupPointVector(i,0)) for i in range(self.MarkupsNode.GetNumberOfFiducials())]
      
      # Polynomial curve fit to see which points deviate the most
      (SrPolynomial, SaPolynomial) = self.CoordsPolyFit(self.MarkupsNode)
      
      RlFitSqErrors = []
      ApFitSqErrors = []
      print " R-L  - A-P"
      print "" 
      
      print "Polynomial fit errors"
      for PointIndex in range(self.MarkupsNode.GetNumberOfFiducials()):
        CurrentPointCoords = self.MarkupsNode.GetMarkupPointVector(PointIndex,0)
        RlFitSqErrors.append((CurrentPointCoords[0] - SrPolynomial(CurrentPointCoords[2]))**2)
        ApFitSqErrors.append((CurrentPointCoords[1] - SaPolynomial(CurrentPointCoords[2]))**2)
        print RlFitSqErrors[-1], "  ", ApFitSqErrors[-1]
   
      print "Mean: ", np.mean(RlFitSqErrors), " - ", np.mean(ApFitSqErrors) 
      print "StdDev: ", np.std(RlFitSqErrors), "  - ", np.std(ApFitSqErrors)
      print ""
   
      # Try local inter-landmark vector direction changes
      RsAngles = []
      AsAngles = []
      sVector = [0.0,0.0,1.0]    # For angle measurment with cross-product
      
      print "Inter-landmark angles (relative to vertical)"
      for i, LabelCoords in enumerate(LabelsCoords[1:], start=1):
        #InterLandmarkVector = [0.0,0.0,0.0]
        RsVector = [(LabelsCoords[i-1][1][0] - LabelCoords[1][0]), 0.0, (LabelsCoords[i-1][1][2] - LabelCoords[1][2])]
        AsVector = [0.0, (LabelsCoords[i-1][1][1] - LabelCoords[1][1]), (LabelsCoords[i-1][1][2] - LabelCoords[1][2])]
        RsDot = np.dot(RsVector, sVector) / np.linalg.norm(RsVector)
        AsDot = np.dot(AsVector, sVector) / np.linalg.norm(AsVector)
        RsAngle = np.math.acos(RsDot) * np.sign(RsVector[0])
        AsAngle = np.math.acos(AsDot) * np.sign(AsVector[1])
        RsAngles.append(RsAngle)
        AsAngles.append(AsAngle)
        print RsAngle, "  ", AsAngle
        
      print "Mean (of abs): ", np.mean([abs(Angle) for Angle in RsAngles]), "  - ", np.mean([abs(Angle) for Angle in AsAngles])
      print "StdDev: ", np.std(RsAngles), " - ", np.std(AsAngles)
      
      # Check angle/fit criteria for outlier identification
      for i, ((RlFitSqError, ApFitSqError),(RsAngle, AsAngle)) in enumerate(zip(zip(RlFitSqErrors, ApFitSqErrors), zip(RsAngles, AsAngles))[1:-1], start=1):
        CurrentRsAngle = RsAngles[i]
        if np.sign(CurrentRsAngle) != np.sign(RsAngles[i-1]) and np.sign(CurrentRsAngle) != np.sign(RsAngles[i+1]):
          # Outlier detected
          print "Outlier detected at ", self.MarkupsNode.GetNthFiducialLabel(i-1), " of ", self.MarkupsNode.GetName()
          self.OutlierIdentificationVotes[i] += 1
          self.MarkupsNode.SetNthFiducialSelected(i,0)    # Currently, unselect outliers
      
      print self.OutlierIdentificationVotes
      print ""
      
      return self.MarkupsNode
   
  def __init__(self, Markups, kmWindowSize, PolyFitDegree, Specificity):
    #  Instantiation of PatientModel here performs initial left-right classification into self.PatientModel.LeftSide etc...
    #self.kmWindowSize = kmWindowSize
    #self.PolyDegree = PolyDegree
    #self.Specificity = Specificity
    self.PatientModel = self.PatientTransverseProcesses(self, Markups, kmWindowSize, PolyFitDegree, Specificity)
   
  def ApplyNamingConvention(self, MarkupsNode):   # Renames points in MarkupsNode according to convention used by ModelToPatientRegistration extension
    # List names expected by ModelToPatientRegistration
    ExpectedLabels = ['T1L','T1R','T2L','T2R','T3L','T3R','T4L','T4R','T5L','T5R','T6L','T6R','T7L','T7R','T8L','T8R','T9L','T9R',\
      'T10L','T10R','T11L','T11R','T12L','T12R','L1L','L1R','L2L','L2R','L3L','L3R','L4L','L4R','L5L','L5R']
      
    # Get names currently used by MarkupsNode
    #IndicesLabels = [(i, MarkupsNode.GetNthFiducialLabel(i)) for i in range(MarkupsNode.GetNumberOfFiducials())]
    Labels = [MarkupsNode.GetNthFiducialLabel(i) for i in range(MarkupsNode.GetNumberOfFiducials())]
    
    NewNode = slicer.vtkMRMLMarkupsFiducialNode()
    NewNode.SetName(MarkupsNode.GetName())
    for Coord in [MarkupsNode.GetMarkupPointVector(i,0) for i in range(MarkupsNode.GetNumberOfFiducials())]:
      NewNode.AddFiducialFromArray(Coord)
    
    # Find some label in MarkupsNode which contains corresponds to some ExpectedLabel
    for i, ExpectedLabel in enumerate(ExpectedLabels):
      # Check each label of MarkupsNode against all ExpectedLabels
      for j, Label in enumerate(Labels):
        if Label.__contains__(ExpectedLabel):
          if j > i: # Check compatability of point set with naming convention
            print "Error - Tried to name more vertebrae than currently recognized"
            print " Does scan include cervical spine?"
            print " Maybe point sets were over-extrapolated."
            print " Leaving names unmodified."
            return
          else:
            #CurrentLabel = ExpectedLabels[i-j]
            for k, ExpectedLabel in enumerate(ExpectedLabels[i:MarkupsNode.GetNumberOfFiducials()+2+i]):
              #print ExpectedLabel
              NewNode.SetNthFiducialLabel(k, ExpectedLabel)
            return NewNode
              
    return
   
class PreProcessLandmarksTest(ScriptedLoadableModuleTest):
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
    self.test_PreProcessLandmarks1()

  def test_PreProcessLandmarks1(self):
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
    logic = PreProcessLandmarksLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')