# Open qt window to get user input
import sys, PyQt5
from PyQt5.QtWidgets import QApplication, QWidget

WindowApplication = QApplication(sys.argv)
Window = QWidget()
Window.resize(320,240)
Window.setWindowTitle("Wiki Industry Analysis")
Window.show()

sys.exit(WindowApplication.exec_())

# Generate and execute WikiHtmlGrabber.bat file with console-line curl-commands


# Read in scrapped wiki page and find table


# Perform analysis

