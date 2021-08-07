#!/usr/bin/env python

import os 
import sys 
import logging
from PyQt5 import QtCore, QtWidgets,QtGui

#============================= Function =================================
##logging info
DEBUG="" #change it when debugging
logFormat = "%(asctime)s [%(levelname)s] %(message)s"
level = "DEBUG" if DEBUG != "" else "INFO"
logging.basicConfig( stream=sys.stderr, level=level, format=logFormat )
#========================================================================

class Screenshot(QtWidgets.QWidget):
    def __init__(self):
        super(Screenshot, self).__init__()
        self.screenshotLabel = QtWidgets.QLabel()
        self.screenshotLabel.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding)
        self.screenshotLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.screenshotLabel.setMinimumSize(240, 160)
        self.createOptionsGroupBox()
        self.createButtonsLayout()
        mainLayout = QtWidgets.QVBoxLayout()
        mainLayout.addWidget(self.screenshotLabel)
        mainLayout.addWidget(self.optionsGroupBox)
        mainLayout.addLayout(self.buttonsLayout)
        self.setLayout(mainLayout)
        self.shootScreen()
        self.delaySpinBox.setValue(0)
        self.setWindowTitle("Screenshot")
        self.resize(300, 200)
        
    def resizeEvent(self, event):
        scaledSize = self.originalPixmap.size()
        scaledSize.scale(self.screenshotLabel.size(), QtCore.Qt.KeepAspectRatio)
        if not self.screenshotLabel.pixmap() or scaledSize != self.screenshotLabel.pixmap().size():
            self.updateScreenshotLabel()
            
    def newScreenshot(self):
        if self.hideThisWindowCheckBox.isChecked():
            self.hide()
        self.newScreenshotButton.setDisabled(True)
        QtCore.QTimer.singleShot(self.delaySpinBox.value() * 1000,
                self.shootScreen)
                
    def saveScreenshot(self):
        format = 'png'
        initialPath = QtCore.QDir.currentPath() + "/untitled." + format
        fileName,filetype = QtWidgets.QFileDialog.getSaveFileName(self, "Save As",
                initialPath,
                "%s Files (*.%s);;All Files (*)" % (format.upper(), format))
        if fileName:
            self.originalPixmap.save(fileName, format)
            logging.info(f"Screenshot has been saved to '{os.path.abspath(fileName)}'")
           
    def shootScreen(self):
        if self.delaySpinBox.value() != 0:
            QtWidgets.qApp.beep()
        # Garbage collect any existing image first.
        self.originalPixmap = None
       
        screen= QtWidgets.QApplication.primaryScreen()
        self.originalPixmap = screen.grabWindow(QtWidgets.QApplication.desktop().winId())
       
        self.updateScreenshotLabel()
        self.newScreenshotButton.setDisabled(False)
        if self.hideThisWindowCheckBox.isChecked():
            self.show()
            
    def updateCheckBox(self):
        if self.delaySpinBox.value() == 0:
            self.hideThisWindowCheckBox.setDisabled(True)
        else:
            self.hideThisWindowCheckBox.setDisabled(False)
            
    def createOptionsGroupBox(self):
        self.optionsGroupBox = QtWidgets.QGroupBox("Options")
        self.delaySpinBox = QtWidgets.QSpinBox()
        self.delaySpinBox.setSuffix(" s")
        self.delaySpinBox.setMaximum(60)
        self.delaySpinBox.valueChanged.connect(self.updateCheckBox)
        self.delaySpinBoxLabel = QtWidgets.QLabel("Screenshot Delay:")
        self.hideThisWindowCheckBox = QtWidgets.QCheckBox("Hide This Window")
        optionsGroupBoxLayout = QtWidgets.QGridLayout()
        optionsGroupBoxLayout.addWidget(self.delaySpinBoxLabel, 0, 0)
        optionsGroupBoxLayout.addWidget(self.delaySpinBox, 0, 1)
        optionsGroupBoxLayout.addWidget(self.hideThisWindowCheckBox, 1, 0, 1, 2)
        self.optionsGroupBox.setLayout(optionsGroupBoxLayout)
        
    def createButtonsLayout(self):
        self.newScreenshotButton = self.createButton("New Screenshot",
                self.newScreenshot)
        self.saveScreenshotButton = self.createButton("Save Screenshot",
                self.saveScreenshot)
        self.quitScreenshotButton = self.createButton("Quit", self.close)
        self.buttonsLayout = QtWidgets.QHBoxLayout()
        self.buttonsLayout.addStretch()
        self.buttonsLayout.addWidget(self.newScreenshotButton)
        self.buttonsLayout.addWidget(self.saveScreenshotButton)
        self.buttonsLayout.addWidget(self.quitScreenshotButton)
    def createButton(self, text, member):
        button = QtWidgets.QPushButton(text)
        button.clicked.connect(member)
        return button
    def updateScreenshotLabel(self):
        self.screenshotLabel.setPixmap(self.originalPixmap.scaled(
                self.screenshotLabel.size(), QtCore.Qt.KeepAspectRatio,
                QtCore.Qt.SmoothTransformation))


if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    screenshot = Screenshot()
    screenshot.show()
    sys.exit(app.exec_())
