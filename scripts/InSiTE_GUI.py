#!/usr/bin/python

import sys

from PyQt5.QtWidgets import (QWidget, QLabel, QHBoxLayout,QVBoxLayout,
                             QComboBox, QApplication, QCheckBox)


class Example(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        filetype = QHBoxLayout() #filetype container for dropdown menu and label
        ticks= QHBoxLayout() #ticks menu for check boxes
        vbox = QVBoxLayout() # vbox overall layout
        vbox.addLayout(filetype)
        vbox.addLayout(ticks)
        self.setContentsMargins(20, 20, 20, 20)
        self.setLayout(vbox)
        self.move(300, 300)
        self.setWindowTitle('InSiTE')

        self.lbl = QLabel('FASTQ', self) #starting label
        combo = QComboBox(self) #dropdown menu
        combo.addItem('FASTQ')
        combo.addItem('FASTA')
        combo.addItem('SAM')
        combo.addItem('CSV')
        combo.activated[str].connect(self.onActivated) #when dropdown menu is changed

        filetype.addWidget(combo)
        filetype.setSpacing(20)
        filetype.addStretch(1)
        filetype.addWidget(self.lbl)

        box = QCheckBox("Awesome?", self)
        box.stateChanged.connect(self.clickBox)
        ticks.addWidget(box)
        #self.b.move(20, 40)
        #self.b.resize(320, 40)
        self.show()
    def clickBox(self, state):

        if state == QtCore.Qt.Checked:
            print('Checked')
        else:
            print('Unchecked')
    def onActivated(self, text):
        self.lbl.setText('new text')#text)
        self.lbl.adjustSize()
        input_type=text

def main():
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

