from __main__ import vtk, qt, ctk, slicer

#
# ModelToPatientRegistration
#

class ModelToPatientRegistration:
  def __init__(self, parent):
    parent.title = "ModelToPatientRegistration"
    parent.categories = ["Landmark Registration"]
    parent.dependencies = []
    parent.contributors = ["Ben Church - Queen's University, PerkLab"]
    parent.helpText = """
    This scripted loadable module was initially created for the B-spline (?)
    spinal model registration project.
    """
    parent.acknowledgementText = """ """ 

#
# ModelToPatientRegistrationWidget
#

class ModelToPatientRegistrationWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qtQVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parnet.show()

  def setup(self):
    # Another panel for modifying landmark sets, deleting points and such
    LandmarkModification = ctk.ctkCollapsibleButton()
    LandmarkModification.text = "Landmark Set Modification"
    self.layout.addWidget(LandmarkModification)
    LandmarkModificationLayout = qt.QFormLayout(LandmarkModification)
    LandmarkModification.collapsed = True
    
    self.ToModifySelector = slicer.qMRMLNodeComboBox()
    self.ToModifySelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.ToModifySelector.selectNodeUponCreation = False
    self.ToModifySelector.enabled = True
    self.ToModifySelector.addEnabled = True
    self.ToModifySelector.noneEnabled = True
    self.ToModifySelector.removeEnabled = True
    self.ToModifySelector.renameEnabled = True
    self.ToModifySelector.toolTip = "Choose which landmark set to perform operation on"
    LandmarkModificationLayout.addRow("Input node:", self.ToModifySelector)
    self.ToModifySelector.setMRMLScene(slicer.mrmlScene)
    
    self.ModifiedStorageNode = slicer.qMRMLNodeComboBox()
    self.ModifiedStorageNode.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.ModifiedStorageNode.selectNodeUponCreation = True
    self.ModifiedStorageNode.enabled = True
    self.ModifiedStorageNode.addEnabled = True
    self.ModifiedStorageNode.noneEnabled = True
    self.ModifiedStorageNode.removeEnabled = True
    self.ModifiedStorageNode.renameEnabled = True
    self.ModifiedStorageNode.toolTip = "Choose markups node to contain modified landmarks"
    LandmarkModificationLayout.addRow("Output node:", self.ModifiedStorageNode)
    self.ModifiedStorageNode.setMRMLScene(slicer.mrmlScene)
    
    self.RemoveThirdLandmarks = qt.QCheckBox("Remove every third vertebra's landmarks")
    self.RemoveThirdLandmarks.toolTip = "Creates markups node with one third of original node's vertebrae removed"
    self.RemoveThirdLandmarks.enabled = True
    LandmarkModificationLayout.addRow(self.RemoveThirdLandmarks)
    self.RemoveThirdLandmarks.connect('toggled(bool)', self.onReconfigure)
    
    self.RemoveTwoThirdsLandmarks = qt.QCheckBox("Remove two of every three vertebrae's landmarks")
    self.RemoveTwoThirdsLandmarks.toolTip = "Creates markups node with one third of original node's vertebrae remaining"
    self.RemoveTwoThirdsLandmarks.enabled = True
    LandmarkModificationLayout.addRow(self.RemoveTwoThirdsLandmarks)
    self.RemoveTwoThirdsLandmarks.connect('toggled(bool)', self.onReconfigure)
    
    self.RemoveHalfLandmarks = qt.QCheckBox("Remove every other vertebra's landmarks")
    self.RemoveHalfLandmarks.toolTip = "Creates markups node with half the inputs vertebra's landmarks removed"
    self.RemoveHalfLandmarks.enabled = True
    LandmarkModificationLayout.addRow(self.RemoveHalfLandmarks)
    self.RemoveHalfLandmarks.connect('toggled(bool)', self.onReconfigure)
    
    self.ModifyLandmarkSet = qt.QPushButton("Modify landmark set")
    self.ModifyLandmarkSet.toolTip = "Run chosen operation on chosen landmark set"
    LandmarkModificationLayout.addRow(self.ModifyLandmarkSet)
    
    self.ToModifySelector.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.ToModifySelector))
    self.ModifiedStorageNode.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.ModifiedStorageNode))
    self.ModifyLandmarkSet.connect('clicked(bool)', self.OnModifyMarkupsNode)
    
    # Set up user interface panel
    moduleInterface = ctk.ctkCollapsibleButton()
    moduleInterface.text = "Model to patient registration"
    self.layout.addWidget(moduleInterface)
    moduleInterfaceLayout = qt.QFormLayout(moduleInterface)
    
    # Input data dropdown list selector
    self.FromInputSelector = slicer.qMRMLNodeComboBox()
    self.FromInputSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.FromInputSelector.selectNodeUponCreation = True
    self.FromInputSelector.enabled  = True
    self.FromInputSelector.addEnabled = True
    self.FromInputSelector.noneEnabled = False
    self.FromInputSelector.removeEnabled = True
    self.FromInputSelector.renameEnabled = True
    self.FromInputSelector.toolTip = "Choose the landmark list to be registered to the other"
    moduleInterfaceLayout.addRow("From-landmark input list:", self.FromInputSelector)
    self.FromInputSelector.setMRMLScene(slicer.mrmlScene)

    # Input data dropdown list selector
    self.ToInputSelector = slicer.qMRMLNodeComboBox()
    self.ToInputSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.ToInputSelector.selectNodeUponCreation = True
    self.ToInputSelector.enabled  = True
    self.ToInputSelector.addEnabled = True
    self.ToInputSelector.noneEnabled = False
    self.ToInputSelector.removeEnabled = True
    self.ToInputSelector.renameEnabled = True
    self.ToInputSelector.toolTip = "Choose the landmark list the other will be registered to"
    moduleInterfaceLayout.addRow("To-landmark input list:", self.ToInputSelector)
    self.ToInputSelector.setMRMLScene(slicer.mrmlScene)

    # Output continer storage selector
    self.FromOutputSelector = slicer.qMRMLNodeComboBox()
    self.FromOutputSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.FromOutputSelector.selectNodeUponCreation = True
    self.FromOutputSelector.enabled  = True
    self.FromOutputSelector.addEnabled = True
    self.FromOutputSelector.noneEnabled = False
    self.FromOutputSelector.removeEnabled = True
    self.FromOutputSelector.renameEnabled = True
    self.FromOutputSelector.toolTip = "Select a vtkMRMLMarkupsFiducialNode to store landmarks to be registered to the other"
    moduleInterfaceLayout.addRow("From-landmark output storage:", self.FromOutputSelector)
    self.FromOutputSelector.setMRMLScene(slicer.mrmlScene)
    
    # Output continer storage selector
    self.ToOutputSelector = slicer.qMRMLNodeComboBox()
    self.ToOutputSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode",]
    self.ToOutputSelector.selectNodeUponCreation = True
    self.ToOutputSelector.enabled  = True
    self.ToOutputSelector.addEnabled = True
    self.ToOutputSelector.noneEnabled = False
    self.ToOutputSelector.removeEnabled = True
    self.ToOutputSelector.renameEnabled = True
    self.ToOutputSelector.toolTip = "Select a vtkMRMLMarkupsFiducialNode to store landmarks the other will be registered to"
    moduleInterfaceLayout.addRow("To-landmark output storage:", self.ToOutputSelector)
    self.ToOutputSelector.setMRMLScene(slicer.mrmlScene)
    
    # Scaling factors option
    self.UseLocalFactors = qt.QCheckBox("Compute scaling factors individually")
    self.UseLocalFactors.toolTip = "Uses information near to each landmark to compute individual scaling factors, rather than one global scaling factor"
    self.UseLocalFactors.enabled = True
    moduleInterfaceLayout.addRow(self.UseLocalFactors)
    
    self.UseAverageAnatomicScaling = qt.QCheckBox("Use average of individual scaling factors")
    self.UseAverageAnatomicScaling.toolTip = "Uses average anatomic scaling factor for uniform scaling of model to patient models"
    self.UseAverageAnatomicScaling.enabled = True
    moduleInterfaceLayout.addRow(self.UseAverageAnatomicScaling)
    
    self.UseVerticalAnatomicScaling = qt.QCheckBox("Use vertical inter-landmark distances to scale anchor offsets")
    self.UseVerticalAnatomicScaling.toolTip = "One possible way of scaling the model spine to the patient's"
    self.UseVerticalAnatomicScaling.enabled = True
    moduleInterfaceLayout.addRow(self.UseVerticalAnatomicScaling)
    
    self.UseScalingPoints = qt.QCheckBox("Use scale points nodes to calculate anchor point offset magnitude")
    self.UseScalingPoints.toolTip = "Uses superior and inferior scaling points placed on vertebral bodies (in seperate nodes from TrX) for anchor point offset scaling"
    self.UseScalingPoints.enabled = True
    moduleInterfaceLayout.addRow(self.UseScalingPoints)
    
    # Find-list-intersection button
    self.GeneratePoints = qt.QPushButton("Generate registration points")
    self.GeneratePoints.toolTip = "Finds intersection of model and patient anatomic landmark sets, then supplements them with anchor points."
    if(self.FromInputSelector.currentNode and self.ToInputSelector.currentNode and (self.FromInputSelector.currentNode != self.ToInputSelector.currentNode)):
      if(self.FromOutputSelector.currentNode and self.ToOutputSelector.currentNode and (self.FromOutputSelector.currentNode != self.ToOutputSelector.currentNode)):
        self.GeneratePoints.enabled = True
      else:
        self.GeneratePoints.enabled = False
    else:
      self.GeneratePoints.enabled = False
    moduleInterfaceLayout.addRow(self.GeneratePoints)
    
    # Reload module button
    self.reloadButton = qt.QPushButton('Reload module')
    moduleInterfaceLayout.addRow(self.reloadButton)
    
    self.FromInputSelector.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.FromInputSelector))
    self.ToInputSelector.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.ToInputSelector))
    self.FromOutputSelector.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.FromOutputSelector))
    self.ToOutputSelector.connect("currentNodeChanged(vtkMRMLMarkupsFiducialNode*)", self.onSelect(self.ToOutputSelector))
    self.GeneratePoints.connect('clicked(bool)', self.onGeneratePoints)
    self.reloadButton.connect('clicked(bool)', self.onReloadButton)
    
  def onSelect(self, ComboBox):
    ComboBox.setMRMLScene(slicer.mrmlScene)
    #logic = ModelToPatientRegistrationLogic(self.FromInputSelector.currentNode, self.ToInputSelector.currentNode, self.FromOutputSelector.currentNode, self.ToOutputSelector.currentNode)
    
  def onReconfigure(self):
    Check = [self.RemoveThirdLandmarks.checkState(), self.RemoveHalfLandmarks.checkState(), self.RemoveTwoThirdsLandmarks.checkState()]
    if Check.count(2) > 1:
      print "One operation may be selected at a time."
      return
    if Check.count(2) == 0:
      print "Select one operation"
      return

    
  def OnModifyMarkupsNode(self):
    logic = MarkupsNodeModificationLogic(self.ToModifySelector.currentNode(), self.ModifiedStorageNode.currentNode())
    if self.RemoveThirdLandmarks.checkState():
      logic.RemoveThirdVertebrae()
    elif self.RemoveHalfLandmarks.checkState():
      logic.RemoveHalfVertebrae()
    elif self.RemoveTwoThirdsLandmarks.checkState():
      logic.RemoveTwoThirdsVertebrae()
    
  def onGeneratePoints(self):
    logic = ModelToPatientRegistrationLogic(self.FromInputSelector.currentNode(), self.ToInputSelector.currentNode(), self.FromOutputSelector.currentNode(), self.ToOutputSelector.currentNode())
    #Reanme run and pass parameters
    logic.AnchorPointSets()
  
  def cleanup(self):
    pass
    
  def onReloadButton(self):
    if(slicer.util.getNode('ModelRegistrationPointsNode') != None):
      slicer.mrmlScene.RemoveNode(slicer.util.getNode('ModelRegistrationPointsNode'))
    if(slicer.util.getNode('PatientRegistrationPointsNode') != None):  
      slicer.mrmlScene.RemoveNode(slicer.util.getNode('PatientRegistrationPointsNode'))
    slicer.util.reloadScriptedModule(slicer.moduleNames.ModelToPatientRegistration)

#
# ModelToPatientRegistrationLogic
#

class ModelToPatientRegistrationLogic:
  
  def __init__(self, OriginalModelPoints, OriginalPatientPoints, ModelRegistrationPointsNode, PatientRegistrationPointsNode):
    
    self.OriginalPatientPoints = OriginalPatientPoints
    self.OriginalModelPoints = OriginalModelPoints
    self.ModelRegistrationPointsNode = ModelRegistrationPointsNode
    self.PatientRegistrationPointsNode = PatientRegistrationPointsNode
    
    """
    # Unit vectors
    self.PatientSupInfVectorsLeft = []
    self.PatientRightLeftVectorsLeft = []
    self.PatientAntPostVectorsLeft = []
    self.PatientSupInfVectorsRight = []
    self.PatientRightLeftVectorsRight = []
    self.PatientAntPostVectorsRight = []
    self.ModelSupInfVectorsLeft = []
    self.ModelRightLeftVectorsLeft = []
    self.ModelAntPostVectorsLeft = []
    self.ModelSupInfVectorsRight = []
    self.ModelRightLeftVectorsRight = []
    self.ModelAntPostVectorsRight = []
    
    # Scaling factors
    self.PatientScalingFactorsLeft = []
    self.PatientScalingFactorsRight = []
    self.ModelScalingFactorsLeft = []
    self.ModelScalingFactorsRight = []

    
    # Input and output (anchor) points
    self.ModelRegistrationPoints = []
    self.ModelRegistrationPointsLeft = []
    self.ModelRegistrationPointsRight = []
    self.PatientRegistrationPoints = []
    self.PatientRegistrationPointsLeft = []
    self.PatientRegistrationPointsRight = []
    self.NamesIntersection = []

    
    # Scale points in (Label, [R,A,S]) format
    self.PatientScalePoints = []
    self.ModelScalePoints = []
    """
    # Reinitialize
    slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().RemoveAllMarkups()
    slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().RemoveAllMarkups()
    
    self.UseVertebraWiseScaling = slicer.modules.ModelToPatientRegistrationWidget.UseLocalFactors.isChecked()
    self.UseAverageScaling = slicer.modules.ModelToPatientRegistrationWidget.UseAverageAnatomicScaling.isChecked()
    self.UseVerticalScaling = slicer.modules.ModelToPatientRegistrationWidget.UseVerticalAnatomicScaling.isChecked()
    self.UseScalingPoints = slicer.modules.ModelToPatientRegistrationWidget.UseScalingPoints.isChecked()
    if self.UseScalingPoints:
      self.PatientScalePointsNode = slicer.util.getNode('ScalePoints' + self.OriginalPatientPoints.GetName()[-3:])
      self.AllModelScalePointsNode = slicer.util.getNode('ModelScalePoints')
    self.GlobalVertebralScalingFactor = 30       # Distance (mm) in anterior direction to duplicate fiducial points for registration accuracy
  
  class Patient:
    def __init__(self, Parent, LandmarksNode, RegistrationPointsNode, ScalePointsNode=None):    # Parent is ModelToPatientRegistrationLogic class
      import numpy as np
      # Sort landmarks, left to right, superior to inferior, into a new node
      self.RawLandmarksNode = LandmarksNode     # 'Raw' meaning may be missing points, or be out of order
      self.LandmarksNode = self.OrderRawLandmarksIntoRegistration()
      self.ScalePointsNode = ScalePointsNode
      
      # Make sure landmarks have valid labels
      if not self.CheckNamingConvention():
        print "Naming convention not followed - Patient object initialized but left empty"
        return
      
      self.AnchorPointsNode = slicer.vtkMRMLMarkupsFiducialNode()
      self.AnchorPointsNode.SetName("PatientAnchorPoints" + LandmarksNode.GetName()[-3:])
      
      self.AnchorOffsetScaleFactors = np.zeros(self.LandmarksNode.GetNumberOfFiducials())
      self.AnchorOffsetDirectionVectors = []
      #for P in range(len(self.AnchorOffsetScaleFactors)):
      #  self.AnchorOffsetDirectionVectors.append(np.zeros(3))
      
      Parent.ComputeOffsetUnitVectors(self)
      
      self.RegistrationPointsNode = RegistrationPointsNode                 # Will contain sorted landmarks and anchor points, in that order
      
    def CheckNamingConvention(self):       # Returns bool
      NamesValid = True
      
      ValidTrPNames = ['T1L','T1R','T2L','T2R','T3L','T3R','T4L','T4R','T5L','T5R','T6L','T6R','T7L','T7R','T8L','T8R',
        'T9L','T9R','T10L','T10R','T11L','T11R','T12L','T12R','L1L','L1R','L2L','L2R','L3L','L3R','L4L','L4R','L5L','L5R']
      for i in range(self.LandmarksNode.GetNumberOfFiducials()-1):
        if self.LandmarksNode.GetNthFiducialLabel(i) not in ValidTrPNames:
          print "TrP naming convention not followed."
          print " Name " + PatientScalePoint[0] + " not allowed"
          NamesValid = False
          return NamesValid
          
      if self.ScalePointsNode != None:
        ValidScalePointNames = ['T1S','T1I','T2S','T2I','T3S','T3I','T4S','T4I','T5S','T5I','T6S','T6I','T7S','T7I','T8S','T8I',
          'T9S','T9I','T10S','T10I','T11S','T11I','T12S','T12I','L1S','L1I','L2S','L2I','L3S','L3I','L4S','L4I','L5S','L5I']
        for i in range(self.ScalePointsNode.GetNumberOfFiducials()-1):
          if self.ScalePointsNode.GetNthFiducialLabel(i) not in ValidScalePointNames:
            print "ScalePoints naming convention not followed."
            print " Name " + self.ScalePointsNode.GetNthFiducialLabel(i) + " not allowed"
            NamesValid = False
            return NamesValid
            
      return NamesValid
      
    def OrderRawLandmarksIntoRegistration(self):              # Returns vtkMRMLMarkupsFiducialNode
      RawThoracicPoints = []   # Will contain landmark points, represented as (Label, [R, A, S]) for sorting
      RawLumbarPoints = []
      for i in range(self.RawLandmarksNode.GetNumberOfFiducials()):
        if self.RawLandmarksNode.GetNthFiducialLabel(i)[0] == "T":
          RawThoracicPoints.append((self.RawLandmarksNode.GetNthFiducialLabel(i), self.RawLandmarksNode.GetMarkupPointVector(i, 0)))
        else:
          RawLumbarPoints.append((self.RawLandmarksNode.GetNthFiducialLabel(i), self.RawLandmarksNode.GetMarkupPointVector(i, 0)))
      SortedThoracicPoints = sorted(RawThoracicPoints, key=lambda tup: int(tup[0][1:-1]))  # Sort by labels alphabetically
      SortedLumbarPoints = sorted(RawLumbarPoints, key=lambda tup: int(tup[0][1:-1]))  # Sort by labels alphabetically
      SortedLandmarks =  SortedThoracicPoints + SortedLumbarPoints
      SortedLandmarksNode = slicer.vtkMRMLMarkupsFiducialNode()
      # Might not need name, not adding this node to scene?
      #SortedLandmarksNode.SetName(self.RawLandmarksNode.GetName() + "Sorted")
      for i, Landmark in enumerate(SortedLandmarks):
        SortedLandmarksNode.AddFiducialFromArray(Landmark[1])
        SortedLandmarksNode.SetNthFiducialLabel(i, Landmark[0])
      return SortedLandmarksNode
  
  class Model:
    def __init__(self, Parent, Patient, FullLandmarksNode, RegistrationPointsNode, ScalePointsNode=None):
      import numpy as np
      self.FullLandmarksNode = FullLandmarksNode    # Contains all landmarks in UsLandmarks_Atlas, not trimmed to correspond to patient
      self.FullScalePointsNode = ScalePointsNode
      (self.LandmarksNode, self.ScalePointsNode) = self.EstablishCorrespondence(Patient)
      
      self.AnchorPointsNode = slicer.vtkMRMLMarkupsFiducialNode()
      self.AnchorPointsNode.SetName("ModelAnchorPoints" + Patient.RawLandmarksNode.GetName()[-3:])
      
      self.AnchorOffsetScaleFactors = np.zeros(self.LandmarksNode.GetNumberOfFiducials())
      self.AnchorOffsetDirectionVectors = []
      #for P in range(len(self.AnchorOffsetScaleFactors)):
       # self.AnchorOffsetDirectionVectors.append(np.zeros(3))
      
      Parent.ComputeOffsetUnitVectors(self)
      
      self.RegistrationPointsNode = RegistrationPointsNode
      
    def EstablishCorrespondence(self, Patient):     # Returns (vtkMRMLMarkupsFiducialNode(), vtkMRMLMarkupsFiducialNode()) ---- (LandmarksNode, ScalePointsNode)
      CorrespondingLandmarksNode = slicer.vtkMRMLMarkupsFiducialNode()
      # Might not need to set name, not adding to scene
      #CorrespondingLandmarksNode.SetName("aaa")
      PatientLandmarkPointCounter = 0
      for P in range(self.FullLandmarksNode.GetNumberOfFiducials()):
        if Patient.LandmarksNode.GetNthFiducialLabel(PatientLandmarkPointCounter) == self.FullLandmarksNode.GetNthFiducialLabel(P):
          CorrespondingLandmarksNode.AddFiducialFromArray(self.FullLandmarksNode.GetMarkupPointVector(P,0))
          CorrespondingLandmarksNode.SetNthFiducialLabel(PatientLandmarkPointCounter, Patient.LandmarksNode.GetNthFiducialLabel(PatientLandmarkPointCounter))
          PatientLandmarkPointCounter += 1
        
      CorrespondingScalePointsNode = slicer.vtkMRMLMarkupsFiducialNode()
      PatientScalePointCounter = 0
      for P in range(self.FullScalePointsNode.GetNumberOfFiducials()):
        if Patient.ScalePointsNode.GetNthFiducialLabel(PatientScalePointCounter) == self.FullScalePointsNode.GetNthFiducialLabel(P):
          CorrespondingScalePointsNode.AddFiducialFromArray(self.FullScalePointsNode.GetMarkupPointVector(P, 0))
          CorrespondingScalePointsNode.SetNthFiducialLabel(PatientLandmarkPointCounter, self.FullScalePointsNode.GetNthFiducialLabel(P))
          PatientScalePointCounter += 1
      #print CorrespondingLandmarksNode
      return (CorrespondingLandmarksNode, CorrespondingScalePointsNode)
      
  def AnchorPointSets(self):
    import math, numpy
    
    # Anatomy models, patient and average model
    self.PatientAnatomy = self.Patient(self, self.OriginalPatientPoints, self.PatientRegistrationPointsNode, self.PatientScalePointsNode)
    self.ModelAnatomy = self.Model(self, self.PatientAnatomy, self.OriginalModelPoints, self.ModelRegistrationPointsNode, self.AllModelScalePointsNode)
    
    #self.ComputeOffsetUnitVectors(self.ModelAnatomy)
    #self.ComputeOffsetUnitVectors(self.PatientAnatomy)
    
    if self.UseVertebraWiseScaling or self.UseAverageScaling:
      self.ComputeLocalAnatomicScalingFactors()
      
    if self.UseScalingPoints:
      self.ComputeScalingFactorsFromScalePoints(self.ModelAnatomy)
      self.ComputeScalingFactorsFromScalePoints(self.PatientAnatomy)
    
    self.AddAnchorPoints(self.ModelAnatomy)
    self.AddAnchorPoints(self.PatientAnatomy)
    
    """
    # Construct a point set at patient-spine space of points with correspondence
    self.AnchorPatientSpine()

    # Construct a point set at average-spine model of points with correspondence
    self.AnchorModelSpine()
    
    
    if(len(self.NamesIntersection) == 0):   # For whatever reason, no points in common could be found between the CT and atlas
      print "Warning - intersection of CT and atlas points is the empty set."
      print "   Maybe points are named incorrectly."
      print "   Correspondence lists not genereated."
      return
    
    print "\n Intersection of CT and atlas landmarks:"
    for i, name in enumerate(self.NamesIntersection[:-1]):
      if((name[:-1] == self.NamesIntersection[i+1][:-1]) and (name[-1] != self.NamesIntersection[i+1][-1])):
        print name + "  " + self.NamesIntersection[i+1]
          
      else:
        if(i>0 and ((name[:-1] + "L" not in self.NamesIntersection) or (name[:-1] + "R" not in self.NamesIntersection))):
          print name

    # This should fix bug where last landmark is not printed if it has no partner
    if(self.NamesIntersection[-1][:-1] != self.NamesIntersection[-2][:-1]): # and (self.NamesIntersection[-1][-1] != self.NamesIntersection))
      print self.NamesIntersection[-1]
    return    
  """
  """
  def CheckNamingConvention(self, Anatomy):

    self.PatientPointNames = [x[0] for x in self.PatientRegistrationPoints]
    for PatientScalePoint in self.PatientRegistrationPoints:
      if PatientScalePoint[0] not in self.ValidTrPNames:

      if(self.PatientPointNames.count(PatientScalePoint[0]) != 1):
        print "Duplicated point present in input."
        NamesValid = False
        return NamesValid
    
    if self.UseScalingPoints:
      for PatientScalePoint in self.PatientScalePoints:
        if PatientScalePoint[0] not in self.ValidScalePointNames:
          print "Scale point naming convention not followed."
          print " Name " + PatientScalePoint[0] + " not allowed."
          NamesValid = False
          return NamesValid
    return NamesValid
    

  def OrderPoints(self):
    self.ThoracicPoints = []
    self.LumbarPoints = []
    for PatientPoint in self.PatientRegistrationPoints:
      if PatientPoint[0][0] == "T":
        self.ThoracicPoints.append((PatientPoint[0][1:], PatientPoint[1], int(PatientPoint[0][1:-1])))
      else:
        self.LumbarPoints.append((PatientPoint[0][1:], PatientPoint[1], int(PatientPoint[0][1:-1])))
        
    self.ThoracicPoints.sort(key=lambda tup: tup[2])
    self.LumbarPoints.sort(key=lambda tup: tup[2])
    for i, ThoracicPoint in enumerate(self.ThoracicPoints):
      self.ThoracicPoints[i] = ("T" + ThoracicPoint[0], ThoracicPoint[1])
    for i, LumbarPoint in enumerate(self.LumbarPoints):
      self.LumbarPoints[i] = ("L" + LumbarPoint[0], LumbarPoint[1])
    
    self.PatientRegistrationPoints = self.ThoracicPoints + self.LumbarPoints
    self.CtNames = [x[0] for x in self.PatientRegistrationPoints]
    return

    
  def EstablishCorrespondence(self):
    for i, PatientPoint in enumerate(self.PatientRegistrationPoints):
      if(PatientPoint[0][-1] == "L"):
        self.PatientRegistrationPointsLeft.append(PatientPoint)
      else:
        self.PatientRegistrationPointsRight.append(PatientPoint)
      slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().AddFiducial(PatientPoint[1][0], PatientPoint[1][1], PatientPoint[1][2])
      slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().SetNthFiducialLabel(i, PatientPoint[0])
  
    for i in range(self.OriginalModelPoints.GetNumberOfFiducials()):
      CurrentLabel = self.OriginalModelPoints.GetNthFiducialLabel(i)
      CurrentPoint = self.OriginalModelPoints.GetMarkupPointVector(i,0)
      if(CurrentLabel[1] == "0"):
        if((CurrentLabel[0]+CurrentLabel[2:]) in self.CtNames):
          self.NamesIntersection.append((CurrentLabel[0]+CurrentLabel[2:]))
          self.ModelRegistrationPoints.append((CurrentLabel[0]+CurrentLabel[2:], [CurrentPoint[0], CurrentPoint[1], CurrentPoint[2]]))
          slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
          slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().GetNumberOfFiducials() - 1, self.ModelRegistrationPoints[-1][0])
          if(CurrentLabel[-1] == "L"):
            self.ModelRegistrationPointsLeft.append((CurrentLabel[0]+CurrentLabel[2:], CurrentPoint))
          else:
            self.ModelRegistrationPointsRight.append((CurrentLabel[0]+CurrentLabel[2:], CurrentPoint))
      else:
        if(CurrentLabel in self.CtNames):
          self.NamesIntersection.append(CurrentLabel)
          self.ModelRegistrationPoints.append((CurrentLabel,[CurrentPoint[0], CurrentPoint[1], CurrentPoint[2]]))
          slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
          slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().GetNumberOfFiducials() - 1, CurrentLabel)
          if(CurrentLabel[-1] == "L"):
            self.ModelRegistrationPointsLeft.append((CurrentLabel,[CurrentPoint[0], CurrentPoint[1], CurrentPoint[2]]))
          else:
            self.ModelRegistrationPointsRight.append((CurrentLabel,[CurrentPoint[0], CurrentPoint[1], CurrentPoint[2]]))
            
    if self.UseScalingPoints:
      TopScalePointStart = self.ValidScalePointNames.index(self.PatientScalePoints[0][0])
      for i in range(TopScalePointStart, TopScalePointStart + len(self.PatientScalePoints)):
        CurrentLabel = self.AllModelScalePointsNode.GetNthFiducialLabel(i)
        CurrentPoint = self.AllModelScalePointsNode.GetMarkupPointVector(i, 0)
        self.ModelScalePoints.append((CurrentLabel, CurrentPoint))
  """
  
  def ComputeOffsetUnitVectors(self, Anatomy):
    import numpy as np
    # Treat top (superior-most) boundary condition
    # This method can also search through the rest of the spine, eliminating the boundary condition
    LeftTopSearch = 0
    #print Anatomy.LandmarksNode
    while Anatomy.LandmarksNode.GetNthFiducialLabel(LeftTopSearch)[-1] != "L":
      LeftTopSearch += 1
    TopLeftPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(LeftTopSearch), Anatomy.LandmarksNode.GetMarkupPointVector(LeftTopSearch, 0))
    SecondLeftSearch = LeftTopSearch + 1
    while Anatomy.LandmarksNode.GetNthFiducialLabel(SecondLeftSearch)[-1] != "L":
      SecondLeftSearch += 1
    SecondLeftPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(SecondLeftSearch), Anatomy.LandmarksNode.GetMarkupPointVector(SecondLeftSearch, 0))
    
    RightTopSearch = 0
    while Anatomy.LandmarksNode.GetNthFiducialLabel(RightTopSearch)[-1] != "R":
      RightTopSearch += 1
    TopRightPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(RightTopSearch), Anatomy.LandmarksNode.GetMarkupPointVector(RightTopSearch, 0))
    SecondRightSearch = RightTopSearch + 1
    while Anatomy.LandmarksNode.GetNthFiducialLabel(SecondRightSearch)[-1] != "R":
      SecondRightSearch += 1
    SecondRightPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(SecondRightSearch), Anatomy.LandmarksNode.GetMarkupPointVector(SecondRightSearch, 0))
    
    print (TopLeftPoint[0], TopRightPoint[0])
    
    # Top-end boundary requires only Top and Middle points
    CurrentAxialVector = [0, 0, 0]
    CurrentLateralVector = [0, 0, 0]
    for dim in range(3):
      CurrentAxialVector[dim] = ((TopLeftPoint[1][dim] - SecondLeftPoint[1][dim]) + (TopRightPoint[1][dim] - SecondRightPoint[1][dim])) / 2.0
      CurrentLateralVector[dim] = TopRightPoint[1][dim] - TopLeftPoint[1][dim]
    
    CurrentOffsetVector = np.cross(CurrentAxialVector, CurrentLateralVector)
    OffsetVectorNorm = np.linalg.norm(CurrentOffsetVector)
    for dim in range(3):
      CurrentOffsetVector[dim] = CurrentOffsetVector[dim] / OffsetVectorNorm
    print " " + str(CurrentOffsetVector)
    
    # Append OffsetVector twice, once for each top point ASSUMES FIRST TWO POINTS IN LANDMARKS ARE TOP LEFT AND TOP RIGHT - defeats purpose of above search
    # Consider having one offset vector per vertebra
    Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
    Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
    
    # ASSERT top-end boundary condition dealt with - offset direction calculated
    
    for Point in range(2, Anatomy.LandmarksNode.GetNumberOfFiducials() - 2, 2):  # Assumes first point in set is left points - skips boundary conditions
      # ASSUMES ALL POINTS ARE PRESENT
      CurrentLeftPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(Point), Anatomy.LandmarksNode.GetMarkupPointVector(Point, 0))
      
      LeftIndexAbove = self.FindNegihborPointAbove(Anatomy, (CurrentLeftPoint[0], Point))
      LeftPointAbove = (Anatomy.LandmarksNode.GetNthFiducialLabel(LeftIndexAbove), Anatomy.LandmarksNode.GetMarkupPointVector(LeftIndexAbove, 0))
      
      LeftIndexBelow = self.FindNegihborPointBelow(Anatomy, (CurrentLeftPoint[0], Point))
      LeftPointBelow = (Anatomy.LandmarksNode.GetNthFiducialLabel(LeftIndexBelow), Anatomy.LandmarksNode.GetMarkupPointVector(LeftIndexBelow, 0))
      
      NeighborIndexBeside = self.FindNegihborPointBeside(Anatomy, (CurrentLeftPoint[0], Point))
      CurrentRightPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(NeighborIndexBeside), Anatomy.LandmarksNode.GetMarkupPointVector(NeighborIndexBeside, 0))
      
      RightIndexAbove = self.FindNegihborPointAbove(Anatomy, (CurrentRightPoint[0], NeighborIndexBeside))
      RightPointAbove = (Anatomy.LandmarksNode.GetNthFiducialLabel(RightIndexAbove), Anatomy.LandmarksNode.GetMarkupPointVector(RightIndexAbove, 0))
      
      RightIndexBelow = self.FindNegihborPointBelow(Anatomy, (CurrentRightPoint[0], NeighborIndexBeside))
      RightPointBelow = (Anatomy.LandmarksNode.GetNthFiducialLabel(RightIndexBelow), Anatomy.LandmarksNode.GetMarkupPointVector(RightIndexBelow, 0))
      
      print (CurrentLeftPoint[0], CurrentRightPoint[0])
      
      TopAxialVector = [0, 0, 0]
      BottomAxialVector = [0, 0, 0]
      AverageAxialVector = [0, 0, 0]
      CurrentLateralVector = [0, 0, 0]
      for dim in range(3):
        TopAxialVector[dim] = ((LeftPointAbove[1][dim] - CurrentLeftPoint[1][dim]) + (RightPointAbove[1][dim] - CurrentRightPoint[1][dim])) / 2.0
        BottomAxialVector[dim] = ((CurrentLeftPoint[1][dim] - LeftPointBelow[1][dim]) + (CurrentRightPoint[1][dim] - RightPointBelow[1][dim])) / 2.0
        AverageAxialVector[dim] = (TopAxialVector[dim] + BottomAxialVector[dim]) / 2.0
        CurrentLateralVector[dim] = CurrentRightPoint[1][dim] - CurrentLeftPoint[1][dim]
      
      CurrentOffsetVector = np.cross(AverageAxialVector, CurrentLateralVector)
      OffsetVectorNorm = np.linalg.norm(CurrentOffsetVector)
      
      for dim in range(3):
        CurrentOffsetVector[dim] = CurrentOffsetVector[dim] / OffsetVectorNorm
      print " " + str(CurrentOffsetVector)

      # Append OffsetVector twice, once for each top point ASSUMES FIRST TWO POINTS IN LANDMARKS ARE TOP LEFT AND TOP RIGHT - defeats purpose of above search
      # Consider having one offset vector per vertebra
      Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
      Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
    
    
    LeftBottomSearch = Anatomy.LandmarksNode.GetNumberOfFiducials() - 1
    while Anatomy.LandmarksNode.GetNthFiducialLabel(LeftBottomSearch)[-1] != "L":
      LeftBottomSearch -= 1
    LeftBottomPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(LeftBottomSearch), Anatomy.LandmarksNode.GetMarkupPointVector(LeftBottomSearch, 0))
    
    SecondLastLeftIndex = self.FindNegihborPointAbove(Anatomy, (LeftBottomPoint[0], LeftBottomSearch))
    SecondLastLeftPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(SecondLastLeftIndex), Anatomy.LandmarksNode.GetMarkupPointVector(SecondLastLeftIndex, 0))
    
    RightBottomSearch = Anatomy.LandmarksNode.GetNumberOfFiducials() - 1
    while Anatomy.LandmarksNode.GetNthFiducialLabel(RightBottomSearch)[-1] != "R":
      RightBottomSearch -= 1
    RightBottomPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(RightBottomSearch), Anatomy.LandmarksNode.GetMarkupPointVector(RightBottomSearch, 0))
    
    SecondLastRightIndex = self.FindNegihborPointAbove(Anatomy, (RightBottomPoint[0], RightBottomSearch))
    SecondLastRightPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(SecondLastRightIndex), Anatomy.LandmarksNode.GetMarkupPointVector(SecondLastRightIndex, 0))
    
    print (LeftBottomPoint[0], RightBottomPoint[0])
    
    CurrentAxialVector = [0, 0, 0]
    CurrentLateralVector = [0, 0, 0]
    for dim in range(3):
      CurrentAxialVector[dim] = ((SecondLastLeftPoint[1][dim] - LeftBottomPoint[1][dim]) + (SecondLastRightPoint[1][dim] - RightBottomPoint[1][dim])) / 2.0
      CurrentLateralVector[dim] = SecondLastRightPoint[1][dim] - SecondLastLeftPoint[1][dim]
    
    CurrentOffsetVector = np.cross(CurrentAxialVector, CurrentLateralVector)
    OffsetVectorNorm = np.linalg.norm(CurrentOffsetVector)
    for dim in range(3):
      CurrentOffsetVector[dim] = CurrentOffsetVector[dim] / OffsetVectorNorm
    print " " + str(CurrentOffsetVector)  
      
    Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
    Anatomy.AnchorOffsetDirectionVectors.append(CurrentOffsetVector)
      
  """
    for PatientPoint in self.PatientRegistrationPointsLeft:
      self.PatientSupInfVectorsLeft.append(self.ComputeSupInfVector(PatientPoint, self.PatientRegistrationPointsLeft))
      self.PatientRightLeftVectorsLeft.append(self.ComputeRightLeftVector(PatientPoint, self.PatientRegistrationPointsRight))
      self.PatientAntPostVectorsLeft.append(numpy.cross(self.PatientRightLeftVectorsLeft[-1], self.PatientSupInfVectorsLeft[-1]))
    for PatientPoint in self.PatientRegistrationPointsRight:
      self.PatientSupInfVectorsRight.append(self.ComputeSupInfVector(PatientPoint, self.PatientRegistrationPointsRight))
      self.PatientRightLeftVectorsRight.append(self.ComputeRightLeftVector(PatientPoint, self.PatientRegistrationPointsLeft))
      self.PatientAntPostVectorsRight.append(numpy.cross(self.PatientRightLeftVectorsRight[-1], self.PatientSupInfVectorsRight[-1]))
    for ModelPoint in self.ModelRegistrationPointsLeft:
      self.ModelSupInfVectorsLeft.append(self.ComputeSupInfVector(ModelPoint, self.ModelRegistrationPointsLeft))
      self.ModelRightLeftVectorsLeft.append(self.ComputeRightLeftVector(ModelPoint, self.ModelRegistrationPointsRight))
      self.ModelAntPostVectorsLeft.append(numpy.cross(self.ModelRightLeftVectorsLeft[-1], self.ModelSupInfVectorsLeft[-1]))
    for ModelPoint in self.ModelRegistrationPointsRight:
      self.ModelSupInfVectorsRight.append(self.ComputeSupInfVector(ModelPoint, self.ModelRegistrationPointsRight))
      self.ModelRightLeftVectorsRight.append(self.ComputeRightLeftVector(ModelPoint, self.ModelRegistrationPointsLeft))
      self.ModelAntPostVectorsRight.append(numpy.cross(self.ModelRightLeftVectorsRight[-1], self.ModelSupInfVectorsRight[-1]))
  """
  def FindNegihborPointBelow(self, Anatomy, Point):   # Point is of the form ('Label', index) --- returns PointSearchIndex, indicating neighbor location
    PointSearchIndex = Point[1] + 1                   # Start at Point location, for speed
    while Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchIndex)[-1] != Point[0][-1]:
      #print Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchIndex)
      if PointSearchIndex >= Anatomy.LandmarksNode.GetNumberOfFiducials():
        print "ERROR - could not find neighbor below " + Point[0]
        print " Returning nothing - results invalid"
        return
      PointSearchIndex += 1
    return PointSearchIndex
    
  def FindNegihborPointAbove(self, Anatomy, Point):   # Point is of the form ('Label', index) --- returns PointSearchIndex, indicating neighbor location
    PointSearchIndex = Point[1] - 1                   # Start at Point location, for speed
    while Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchIndex)[-1] != Point[0][-1]:
      #print (Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchIndex), Point)
      if PointSearchIndex < 0:
        print "ERROR - could not find neighbor above " + Point[0]
        print " Returning nothing - results invalid"
        return
      PointSearchIndex -= 1
    return PointSearchIndex
    
  def FindNegihborPointBeside(self, Anatomy, Point):
    PointSearchUpIndex = Point[1] - 1       # Search out from start point, up and down, for speed
    PointSearchDownIndex = Point[1] + 1
    #print "SearchUp: " + str((Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchUpIndex)[:-1], Point[0][:-1]))
    #print "SearchDown" + str((Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchDownIndex)[:-1], Point[0][:-1]))
    while Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchUpIndex)[:-1] != Point[0][:-1] and Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchDownIndex)[:-1] != Point[0][:-1]:
      #print "SearchUp: " + str((Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchUpIndex), Point[0]))
      #print "SearchDown: " + str((Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchDownIndex), Point[0]))
      if PointSearchUpIndex > 0:
        PointSearchUpIndex -= 1
      if PointSearchDownIndex < Anatomy.LandmarksNode.GetNumberOfFiducials()-1:
        PointSearchDownIndex += 1
      #print "SearchUpIndex: " + str(PointSearchUpIndex)
      #print "SearchDownIndex: " + str(PointSearchDownIndex)
    if Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchUpIndex)[:-1] == Point[0][:-1]:
      return PointSearchUpIndex
    if Anatomy.LandmarksNode.GetNthFiducialLabel(PointSearchDownIndex)[:-1] == Point[0][:-1]:
      return PointSearchDownIndex
    print "ERROR - could not find symmetric partner for " + Point[0]
    print " Returning nothing - results invalid"
    return
  
  """
  def ComputeSupInfVector(self, Point, PointSet):
    import math
    SupInfVector = [0, 0, 0]
    
    if Point[0] == PointSet[0][0]:
      # Boundary conition, top point
      for dim in range(3):
        SupInfVector[dim] += PointSet[1][1][dim] - Point[1][dim]
    elif Point[0] == PointSet[-1][0]:
      # Boundary condition, bottom point
      for dim in range(3):
        SupInfVector[dim] += Point[1][dim] - PointSet[-2][1][dim]
    else:
      # Non-bounary, expecting a point above and below
      SearchIterator = 1
      while Point[0] != PointSet[SearchIterator][0]: # Search by label until point is found
        SearchIterator += 1
        if SearchIterator > len(PointSet):
          print "Error - could not compute SupInfVector for point " + Point[0]
          print " Point not found in registration points - Returning 0-vector."
          return [0,0,0]
      for dim in range(3):
        # Use average since we have points above and below
        SupInfVector[dim] += ((Point[1][dim] - PointSet[SearchIterator - 1][1][dim]) + (PointSet[SearchIterator + 1][1][dim] - Point[1][dim])) / 2.0
   
    return SupInfVector
    
  def ComputeRightLeftVector(self, Point, ParallelPointSet):
    import math
    RightLeftVecor = [0, 0, 0]
    SearchIterator = 0
    
    while Point[0][:-1] != ParallelPointSet[SearchIterator][0][:-1]:    # Compare labels from beginning, slow
      SearchIterator += 1
      if SearchIterator > len(ParallelPointSet):
        print "Error - could not find point symmetric to point " + Point[0]
        print " Point not found in registration points - Returning 0-vector."
        return [0,0,0]
    for dim in range(3):
      RightLeftVecor[dim] = ParallelPointSet[SearchIterator][1][dim] - Point[1][dim]
      
    if Point[0][-1] == "R":
      # Ensure vectors are always Left to Right, convention for easy crossing
      for i, Coord in enumerate(RightLeftVecor):
        RightLeftVecor[i] = -1 * Coord
        
    return RightLeftVecor
 """
  def ComputeLocalAnatomicScalingFactors(self):
    import numpy
    if self.UseVerticalScaling:
      # Uses lengths of SupInfVectors as anatomic scaling factors
      SumPatientScalingFactors = 0
      SumModelScalingFactors = 0
      for i in range(len(self.PatientRegistrationPointsLeft)):
        self.PatientScalingFactorsLeft.append(numpy.linalg.norm(self.PatientSupInfVectorsLeft[i]))
        SumPatientScalingFactors += self.PatientScalingFactorsLeft[-1]
      for i in range(len(self.PatientRegistrationPointsRight)):
        self.PatientScalingFactorsRight.append(numpy.linalg.norm(self.PatientSupInfVectorsRight[i]))
        SumPatientScalingFactors += self.PatientScalingFactorsRight[-1]
      for i in range(len(self.ModelRegistrationPointsLeft)):
        self.ModelScalingFactorsLeft.append(numpy.linalg.norm(self.ModelSupInfVectorsLeft[i]))
        SumModelScalingFactors += self.ModelScalingFactorsLeft[-1]
      for i in range(len(self.ModelRegistrationPointsRight)):
        self.ModelScalingFactorsRight.append(numpy.linalg.norm(self.ModelSupInfVectorsRight[i]))
        SumModelScalingFactors += self.ModelScalingFactorsRight[-1]
      
      if self.UseAverageScaling:
        AveragePatientScalingFactor = SumPatientScalingFactors / (len(self.PatientScalingFactorsLeft) + len(self.PatientScalingFactorsRight))
        AverageModelScalingFactor = SumModelScalingFactors / (len(self.ModelScalingFactorsLeft) + len(self.ModelScalingFactorsRight))
        for i in range(len(self.PatientScalingFactorsLeft)):
          self.PatientScalingFactorsLeft[i] = AveragePatientScalingFactor
        for i in range(len(self.PatientScalingFactorsRight)):
          self.PatientScalingFactorsRight[i] = AveragePatientScalingFactor
        for i in range(len(self.ModelScalingFactorsLeft)):
          self.ModelScalingFactorsLeft[i] = AverageModelScalingFactor
        for i in range(len(self.ModelScalingFactorsRight)):
          self.ModelScalingFactorsRight[i] = AverageModelScalingFactor
        
      for i in range(len(self.PatientScalingFactorsLeft)):
        # The multiplication by * (AveragePatientScalingFactor / AverageModelScalingFactor) is the VSF, to factor in the relative lengths of model/patient
        self.PatientScalingFactorsLeft[i] = self.PatientScalingFactorsLeft[i] * (self.PatientScalingFactorsLeft[i] / self.ModelScalingFactorsLeft[i])
      for i in range(len(self.PatientScalingFactorsRight)):
        self.PatientScalingFactorsRight[i] = self.PatientScalingFactorsRight[i] * (self.PatientScalingFactorsRight[i] / self.ModelScalingFactorsRight[i])
 
  def ComputeScalingFactorsFromScalePoints(self, Anatomy):
    # Computes offset magnitudes as distance between CurrentVertebraAveragePostPoint and CurrentVertebraAverageAntPoint
    # Assumes all TrPs have symmetric partners, and all corresponding scale points are present
    import numpy as np
    
    CurrentVertebraAverageAntPoint = [0, 0, 0]
    CurrentVertebraAveragePostPoint = [0, 0, 0]
    AntPostScaleVector = [0, 0, 0]
    for (VertebraAnt, VertebraPost) in zip(range(0, Anatomy.ScalePointsNode.GetNumberOfFiducials() - 1, 2), range(0, Anatomy.LandmarksNode.GetNumberOfFiducials() - 1, 2)):
      SupScalePoint = (Anatomy.ScalePointsNode.GetNthFiducialLabel(VertebraAnt), Anatomy.ScalePointsNode.GetMarkupPointVector(VertebraAnt, 0))
      InfScalePoint = (Anatomy.ScalePointsNode.GetNthFiducialLabel(VertebraAnt + 1), Anatomy.ScalePointsNode.GetMarkupPointVector(VertebraAnt + 1, 0))
      LeftTrP = (Anatomy.LandmarksNode.GetNthFiducialLabel(VertebraPost), Anatomy.LandmarksNode.GetMarkupPointVector(VertebraPost, 0))
      RightTrP = (Anatomy.LandmarksNode.GetNthFiducialLabel(VertebraPost + 1), Anatomy.LandmarksNode.GetMarkupPointVector(VertebraPost + 1, 0))
      for dim in range(3):
        CurrentVertebraAverageAntPoint[dim] = (SupScalePoint[1][dim] + InfScalePoint[1][dim]) / 2.0
        CurrentVertebraAveragePostPoint[dim] = (LeftTrP[1][dim] + RightTrP[1][dim]) / 2.0
        AntPostScaleVector[dim] = CurrentVertebraAverageAntPoint[dim] - CurrentVertebraAveragePostPoint[dim]
      #print np.linalg.norm(AntPostScaleVector)
      Anatomy.AnchorOffsetScaleFactors[VertebraPost] = (np.linalg.norm(AntPostScaleVector)) # Once for the left anchor point
      Anatomy.AnchorOffsetScaleFactors[VertebraPost + 1] = (np.linalg.norm(AntPostScaleVector)) # And once for the right
    #print Anatomy.AnchorOffsetScaleFactors
  
  def AddAnchorPoints(self, Anatomy):
    # Assumes Anatomy's anchor offset magnitudes and directions are instantiated
    print ""
    print Anatomy
    print Anatomy.AnchorOffsetDirectionVectors
    print " "
    print Anatomy.AnchorOffsetScaleFactors
    print " "
    
    #print Anatomy.LandmarksNode
    # Add original landmarks to RegistrationPoints
    for LandmarkIndex in range(Anatomy.LandmarksNode.GetNumberOfFiducials()):
      OriginalLandmarkPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(LandmarkIndex), Anatomy.LandmarksNode.GetMarkupPointVector(LandmarkIndex, 0))
      Anatomy.RegistrationPointsNode.AddFiducialFromArray(OriginalLandmarkPoint[1])
      Anatomy.RegistrationPointsNode.SetNthFiducialLabel(LandmarkIndex, OriginalLandmarkPoint[0])

    # Add anchor points to RegistrationPoints
    for LandmarkIndex in range(Anatomy.LandmarksNode.GetNumberOfFiducials()):
      OriginalLandmarkPoint = (Anatomy.LandmarksNode.GetNthFiducialLabel(LandmarkIndex), Anatomy.LandmarksNode.GetMarkupPointVector(LandmarkIndex, 0))
      CorrespondingAnchorCoord = [0, 0, 0]
      for dim in range(3):
        CorrespondingAnchorCoord[dim] = OriginalLandmarkPoint[1][dim] + (Anatomy.AnchorOffsetDirectionVectors[LandmarkIndex][dim] * Anatomy.AnchorOffsetScaleFactors[LandmarkIndex])
      CorrespondingAnchorPoint = (OriginalLandmarkPoint[0] + "A", CorrespondingAnchorCoord)
      Anatomy.RegistrationPointsNode.AddFiducialFromArray(CorrespondingAnchorCoord)
      Anatomy.RegistrationPointsNode.SetNthFiducialLabel(LandmarkIndex, CorrespondingAnchorPoint[0])
  
  def AnchorPatientSpine(self):
    import numpy
    LeftIt = 0
    RightIt = 0
    for i in range(len(self.PatientRegistrationPoints)):
      if self.PatientRegistrationPoints[i][0][-1] == "L":
        AnchorPoint = (self.PatientRegistrationPoints[i][0] + "A", [0, 0, 0])
        AntPostVectorNorm = numpy.linalg.norm(self.PatientAntPostVectorsLeft[LeftIt])
        for dim in range(3):
          AnchorPoint[1][dim] = self.PatientRegistrationPoints[i][1][dim] + (self.PatientAntPostVectorsLeft[LeftIt][dim] * self.PatientScalingFactorsLeft[LeftIt] / AntPostVectorNorm)
        slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().AddFiducial(AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2])
        slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().GetNumberOfFiducials() - 1, AnchorPoint[0])
        LeftIt += 1
      else: #self.PatientRegistrationPoints[i][0][-1] == "R":
        AnchorPoint = (self.PatientRegistrationPoints[i][0] + "A", [0, 0, 0])
        AntPostVectorNorm = numpy.linalg.norm(self.PatientAntPostVectorsRight[RightIt])
        for dim in range(3):
          AnchorPoint[1][dim] = self.PatientRegistrationPoints[i][1][dim] + (self.PatientAntPostVectorsRight[RightIt][dim] * self.PatientScalingFactorsRight[RightIt] / AntPostVectorNorm)
        slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().AddFiducial(AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2])
        slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.ToOutputSelector.currentNode().GetNumberOfFiducials() - 1, AnchorPoint[0])
        RightIt += 1
      self.PatientRegistrationPoints.append((AnchorPoint[0] + "A", [AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2]]))
  
  def AnchorModelSpine(self):
    import numpy
    LeftIt = 0
    RightIt = 0
    for i in range(len(self.ModelRegistrationPoints)):
      if self.ModelRegistrationPoints[i][0][-1] == "L":
        AnchorPoint = (self.ModelRegistrationPoints[i][0] + "A", [0, 0, 0])
        AntPostVectorNorm = numpy.linalg.norm(self.ModelAntPostVectorsLeft[LeftIt])
        for dim in range(3):
          AnchorPoint[1][dim] = self.ModelRegistrationPoints[i][1][dim] + (self.ModelAntPostVectorsLeft[LeftIt][dim] * self.ModelScalingFactorsLeft[LeftIt] / AntPostVectorNorm)
        slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().AddFiducial(AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2])
        slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().GetNumberOfFiducials() - 1, AnchorPoint[0])
        LeftIt += 1
      else: #self.ModelRegistrationPoints[i][0][-1] == "R":
        AnchorPoint = (self.ModelRegistrationPoints[i][0] + "A", [0, 0, 0])
        AntPostVectorNorm = numpy.linalg.norm(self.ModelAntPostVectorsRight[RightIt])
        for dim in range(3):
          AnchorPoint[1][dim] = self.ModelRegistrationPoints[i][1][dim] + (self.ModelAntPostVectorsRight[RightIt][dim] * self.ModelScalingFactorsRight[RightIt] / AntPostVectorNorm)
        slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().AddFiducial(AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2])
        slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().SetNthFiducialLabel(slicer.modules.ModelToPatientRegistrationWidget.FromOutputSelector.currentNode().GetNumberOfFiducials() - 1, AnchorPoint[0])
        RightIt += 1
      self.ModelRegistrationPoints.append((AnchorPoint[0] + "A", [AnchorPoint[1][0], AnchorPoint[1][1], AnchorPoint[1][2]]))
 
class ModelToPatientRegistrationTest:
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
    self.test_ModelToPatientRegistration1()

  def test_ModelToPatientRegistration1(self):
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
    logic = CtToAtlasRegistrationLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')

class MarkupsNodeModificationLogic:
  def __init__(self, Input, Output):
    self.MarkupsNode = Input
    self.OutputNode = Output
    
  def RemoveThirdVertebrae(self):
    print "Creating landmark set every third vertebra removed"
    for it, i in enumerate(range(0, self.MarkupsNode.GetNumberOfFiducials(), 2)):
      if it % 3 == 0 or it % 3 == 1:
        CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i,0)
        self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
        self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i))
        CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i+1,0)
        self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
        self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i+1))
        print "Landmark " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-2) + " and " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1) + " created."
       
  def RemoveHalfVertebrae(self):
    print "Creating landmark set with half vertebrae removed"
    for i in range(0, self.MarkupsNode.GetNumberOfFiducials(), 4):
      CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i,0)
      self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
      self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i))
      CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i+1,0)
      self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
      self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i+1))
      print "Landmark " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-2) + " and " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1) + " created."
  
  def RemoveTwoThirdsVertebrae(self):
    print "Creating landmark set with two-thirds of vertebrae removed"
    for i in range(0, self.MarkupsNode.GetNumberOfFiducials(), 6):
      CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i,0)
      self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
      self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i))
      CurrentPoint = self.MarkupsNode.GetMarkupPointVector(i+1,0)
      self.OutputNode.AddFiducial(CurrentPoint[0], CurrentPoint[1], CurrentPoint[2])
      self.OutputNode.SetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1, self.MarkupsNode.GetNthFiducialLabel(i+1))
      print "Landmark " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-2) + " and " + self.OutputNode.GetNthFiducialLabel(self.OutputNode.GetNumberOfFiducials()-1) + " created."

  