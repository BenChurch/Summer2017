import sys
import csv
#from PyQt5.QtGui import QIcon
from PyQt5 import QtCore as qc
from PyQt5 import QtGui as qg
from PyQt5 import QtWidgets as qw
#from PyQt5.QtWidgets import QApplication, QWidget

#from PyQt4

####### User defined operation variables #######
# Relative path to InterfaceWindowIcon
IconPath = "./IndustryIcon.png"

####### Interface class #######
class InterfaceWindow(qw.QWidget):
  def __init__(self):
    # Start application
    self.InterfaceApp = qw.QApplication(sys.argv)
    
    # Generate window in application
    super().__init__()
    
    self.InitInterfaceData()
   
   # Initialize window display attributes
    self.InitDisplayAttr()    
    
    self.show()
    sys.exit(self.InterfaceApp.exec_())
    return 
  
  def InitDisplayAttr(self):
    
    self.setWindowTitle('Industry Analysis')
    self.setWindowIcon(qg.QIcon(IconPath))
    
    self.setGeometry(600,400,500,300)
    
    return True
  
  def InitInterfaceData(self):
    # Update list of sovereign states from wikipedia
    self.StateSelectionBox = qw.QComboBox(self)
    self.StateSelectionBox.setEditable(True)
    StateNameStrings = self.GetStateNames()
    self.StateSelectionBox.addItems(StateNameStrings)
    
    self.StateSelectionBox.setGeometry(0,0,140,30)
    self.StateSelectionBox.activated[str].connect(self.OnStateSelectionBoxChanged)
    
    # QLabel to recieve state name
    self.SelectedStateLabel = qw.QLabel(self.StateSelectionBox.currentText(), self)
    self.SelectedStateLabel.setGeometry(170, 0, 140, 30)
    
    return True
  
  def OnStateSelectionBoxChanged(self, StateName):
    self.SelectedStateLabel.setText(self.StateSelectionBox.currentText())
  
    return True
  
  def GetStateNames(self):
    from pyquery import PyQuery as pq
    
    Page = pq(url="https://en.wikipedia.org/wiki/List_of_sovereign_states", opeer=lambda url, **kw: urlopen(url).read())

    AllTables = Page('table')#.filter('.sortable wikitable')
    TableOfStates = AllTables(".sortable.wikitable")
    RowsOfTable = TableOfStates('tr')
    DataOfRows = RowsOfTable('td')
    
    StateNameStrings = []
    for Datum in RowsOfTable.items('span'):
      if Datum.attr("id") == "Other_states":    # Cut off the weird ones
        break
        
      if Datum.attr['id']:
        StateNameStrings.append(Datum.attr['id'])
    
    return StateNameStrings
    
####### Main program #######

MyInterface = InterfaceWindow()

#MyInterface.ReadInStateNames()
