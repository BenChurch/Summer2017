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
    
  def ReadInWikiTable(self):
    from pyquery import PyQuery as pq
    from lxml import etree
    import urllib
    #d = pq("<html></html>")
    #d = pq(etree.fromstring("<html></html>"))
    d = pq(url="https://en.wikipedia.org/wiki/List_of_sovereign_states", opeer=lambda url, **kw: urlopen(url).read())
    d('table').filter('.sortable wikitable')
    #d = pq(filename=AllStatesHtmlFileName, encoding="UTF-8")
    
  
####### Main program #######
MyInterface = Interface()
#MyInterface.WriteBatForDataUpdate()
#MyInterface.ExecBatForDataUpdate()
MyInterface.ReadInWikiTable()