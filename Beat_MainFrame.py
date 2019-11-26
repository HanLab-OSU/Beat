import sys
import os
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication,  QWidget, QPushButton, QLabel, QLineEdit, \
    QFileDialog,QRadioButton, QTextBrowser
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout,QHBoxLayout,QComboBox
from Beat_mac import Beat

class Beat_MainFrame(QWidget):
    def __init__(self):
        super().__init__()
        self.title = "Base Editing Analysis Tool (BEAT)"
        self.setStyleSheet("QLabel {font: 14pt Comic Sans MS}")
        self.setStyleSheet("QGroupBox {font: 14pt Arial}")
        self.left = 500
        self.top = 200
        self.width = 300
        self.height = 250
        self.init_window()


    def init_window(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        v_main = QVBoxLayout()
        v_main.addWidget(self.create_path_reader())
        v_main.addWidget(self.create_arg_reader())
        self.run_btn = QPushButton("RUN")
        self.run_btn.resize(100,50)
        self.run_btn.clicked.connect(self.run_script)
        self.browser = QTextBrowser()
        v_main.addWidget(self.browser)
        v_main.addWidget(self.run_btn)
        self.setLayout(v_main)
        self.show()

    def create_arg_reader(self):
        argbox = QGroupBox("Please enter sequence parameters")
        v_box = QVBoxLayout()

        h_layout_spa = QHBoxLayout()
        spa_label = QLabel("Spacer: ")
        self.spa_input = QLineEdit(self)
        h_layout_spa.addWidget(spa_label)
        h_layout_spa.addWidget(self.spa_input)

        h_layout_pos_base = QHBoxLayout()
        pos_label = QLabel("Posotion: ")
        self.pos_input = QLineEdit(self)
        change_label = QLabel("Base_Change: ")
        # self.change_input = QLineEdit(self)
        self.change_input = QComboBox(self)
        self.change_input.addItems(["AG", "TC", "GA", "CT", "AC", "AT", "GC", "GT", "CA", "CG", "TA", "TG"])
        h_layout_pos_base.addWidget(pos_label)
        h_layout_pos_base.addWidget(self.pos_input)
        h_layout_pos_base.addWidget(change_label)
        h_layout_pos_base.addWidget(self.change_input)

        v_box.addLayout(h_layout_spa)
        v_box.addLayout(h_layout_pos_base)

        argbox.setLayout(v_box)
        return argbox

    def create_path_reader(self):
        pathbox = QGroupBox("Please enter the file(s) for analysis   Renzhi.Han@osumc.edu")
        v_box = QVBoxLayout()

        h_radiobtn = QHBoxLayout()

        self.rad_single = QRadioButton("Single File")
        self.rad_single.setChecked(1)
        self.rad_single.toggled.connect(self.enable_input)

        self.rad_muti = QRadioButton("Batch")
        self.rad_muti.toggled.connect(self.enable_input)

        self.rad_csv = QRadioButton("Batch with different spacer(.csv)")
        self.rad_csv.toggled.connect(self.disable_input)

        h_radiobtn.addWidget(self.rad_single)
        h_radiobtn.addWidget(self.rad_muti)
        h_radiobtn.addWidget(self.rad_csv)

        h_box = QHBoxLayout()
        self.path_input = QLineEdit(self)
        btn = QPushButton("Browse")
        btn.clicked.connect(self.get_file)
        h_box.addWidget(self.path_input)
        h_box.addWidget(btn)

        v_box.addLayout(h_radiobtn)
        v_box.addLayout(h_box)
        pathbox.setLayout(v_box)
        return pathbox

    def get_file(self):
        if self.rad_single.isChecked():
            path = QFileDialog.getOpenFileName(filter='.ab1 file(*.ab1)')
            self.path_input.setText(path[0])
            self.browser.append('----'*15)
            self.browser.append("Calculating the sequencing data of " + self.beat.filename + "...")
            self._update_timer = QtCore.QTimer()
            self._update_timer.timeout.connect(self.run_script)
            self._update_timer.start(10000) # milliseconds            
        if self.rad_muti.isChecked():
            path = QFileDialog.getExistingDirectory()
            self._path = path
            self.path_input.setText(self._path)
            self.browser.append('----'*15)
            self.browser.append("Calculating the sequencing data of ")
            for files in os.listdir(path):
                if files[-4:] == '.ab1':
                    self.browser.append(files)            
            self._update_timer = QtCore.QTimer()
            self._update_timer.timeout.connect(self.run_script)
            self._update_timer.start(10000) # milliseconds             
        if self.rad_csv.isChecked():
            path = QFileDialog.getOpenFileName(filter='csv file(*.csv *.xls, *.xlsx)')
            self.path_input.setText(path[0])
            self.browser.append('----'*15)
            self.browser.append("Run the batch analysis with " + path + ' ...')

    def run_script(self):
        try:
            self.run_btn.setEnabled(False)
            path = self.path_input.text()
            self.beat = Beat(path)
            print("running...")
            
            # self.browser.append('----'*15)
##            if 'result' not in os.listdir('./'):
##                self.browser.append('result folder didn\'t find, creating result under: ' + os.curdir)
##                os.mkdir('result')

            if self.rad_single.isChecked():
                spacer = self.spa_input.text()
                position = int(self.pos_input.text())
                base_change = self.change_input.currentText()
                self.beat.getpathfilename(self.beat.path)
                self.beat.spacer = spacer
                self.beat.position = position
                self.beat.change_base = base_change
                print(self.beat.path, self.beat.filename)                
                self.beat.single_analysis()

            if self.rad_muti.isChecked():
                spacer = self.spa_input.text()
                position = int(self.pos_input.text())
                base_change = self.change_input.currentText()
                self.beat.spacer = spacer
                self.beat.position = position
                self.beat.change_base = base_change
##                self.browser.append("Calculating the sequencing data of: ")
##                for files in os.listdir(self.beat.path):
##                    if files[-4:] == '.ab1':
##                        self.browser.append(files)
                self.beat.batch_analysis(path)

            if self.rad_csv.isChecked():
                print(path)                
                self.beat.csv_analysis(path)

            self.browser.append('----'*15)
            self.browser.append("Saved image and .xslx files in "+ self.beat.path +'/result')
            self.browser.append("Done!")
            self._update_timer.stop()
            
        except Exception as err:
            print(str(err))
            self.run_btn.setEnabled(True)

        self.run_btn.setEnabled(True)

    def get_path(self,path):
        while (1):
            if path[-1] in ['/', '\\']:
                return path
            else:
                path = path[:-1]

    def disable_input(self):
        self.spa_input.setEnabled(False)
        self.pos_input.setEnabled(False)
        self.change_input.setEnabled(False)

    def enable_input(self):
        self.spa_input.setEnabled(True)
        self.pos_input.setEnabled(True)
        self.change_input.setEnabled(True)

if __name__ == '__main__':
    App = QApplication(sys.argv)
    window = Beat_MainFrame()
    sys.exit(App.exec_())
