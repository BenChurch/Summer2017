import sys
import csv
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtGui import QIcon

####### User defined operation variables #######
AllStatesHtmlFileName = "AllStates.html"

class Interface(QWidget):
  def __init__(self):
    # Update list of sovereign states from wikipedia
    return 
  
  def WriteBatForDataUpdate(self):
    HtmlGrabberBat = open("WikiHtmlGrabber.bat", "w")
    HtmlGrabberBat.write("curl https://en.wikipedia.org/wiki/List_of_sovereign_states > " + AllStatesHtmlFileName)
    return 
  
  def ExecBatForDataUpdate(self):
    from subprocess import Popen
    p = Popen(r"WikiHtmlGrabber.bat")
    stdout, stderr = p.communicate()
  
####### Main program #######
MyInterface = Interface()
MyInterface.WriteBatForDataUpdate()
MyInterface.ExecBatForDataUpdate()