# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn

from HazardsFinder import Ui_widget
import suspectscreening_1
import suspectscreening_2
import nontagetscreening_1
import nontagetscreening_2
from datetime import datetime
import os
import pandas as pd
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QTableWidgetItem
from PyQt5.QtWidgets import QDesktopWidget
from PyQt5.QtGui import QCursor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from PyQt5.Qt import QThread

class ThreadOne(QThread):
    def __init__(self):
        super().__init__()
        self.mzml_path_1 = ''
        self.current_path=''
        self.standard_database=''
        self.ms1_error=''
        self.rt_error=''
        self.isotope_sc=''
        self.SIM_sc=''
    def run(self):
        os.chdir(self.current_path)
        suspectscreening_1.screening(self.mzml_path_1, self.current_path, self.standard_database, self.ms1_error, self.rt_error,
                  self.isotope_sc, self.SIM_sc)

class ThreadTwo(QThread):
    def __init__(self):
        super().__init__()
        self.mzml_path_1 = ''
        self.current_path=''
        self.online_database=''
        self.literature_database=''
        self.ms1_error=''
        self.rt_error=''
        self.isotope_sc = ''
        self.SIM_online=''
        self.frgnum_literature=''

    def run(self):
        os.chdir(self.current_path)
        suspectscreening_2.screening(self.mzml_path_1, self.current_path, self.online_database, self.literature_database,self.ms1_error, self.rt_error,self.isotope_sc, self.SIM_online, self.frgnum_literature)



class ThreadThree(QThread):
    def __init__(self):
        super().__init__()
        self.mzml_path_2 = ''
        self.current_path=''
        self.tps_database=''
        self.peaklistpath_1=''
        self.ms1_error=''
        self.rt_error=''
        self.isotope_sc = ''
        self.frgnum_cfmid=''

    def run(self):
        os.chdir(self.current_path)

        nontagetscreening_1.screening(self.mzml_path_2, self.current_path, self.tps_database, self.peaklistpath_1,
                                      self.ms1_error, self.rt_error, self.isotope_sc, self.frgnum_cfmid)


class ThreadFour(QThread):
    def __init__(self):
        super().__init__()
        self.mzml_path_2 = ''
        self.current_path=''
        self.cfi_database=''
        self.peaklistpath_1=''
        self.ms2_error=''
        self.inten_filter=''


    def run(self):
        os.chdir(self.current_path)
        nontagetscreening_2.screening(self.mzml_path_2, self.current_path, self.cfi_database, self.peaklistpath_1,self.ms2_error,self.inten_filter)
class DrawFunction(object):
    def Image_In_Process(self, feature):
        f_num=feature.index_no
        try:
            f_name=feature.Name[f_num]
            f_tRtheor = feature.tRtheor[f_num]
            f_tRexper=feature.tRMS[f_num]
            f_EIC = feature.EIC[f_num]
            f_MSMStheor = feature.MSMStheor[f_num]
            f_MSMSexper = feature.MSMSexper[f_num]
            f_mztheor=feature.mztheor[f_num]
            f_sim=feature.SIM[f_num]
            f_EIC=f_EIC.split(';')
            inten = []
            rt = []
            for w in range(len(f_EIC)):
                rt.append((float(f_EIC[w].split(' ')[0]))/60)
                inten.append(float(f_EIC[w].split(' ')[1]))

            feature.fig1.clear()
            ax1 = feature.fig1.add_subplot(111)
            # 调整图像大小
            ax1.cla()
            rt_max_now = max(rt)
            rt_min_now = min(rt)
            ax1.plot(rt, inten, linewidth=0.8, color='black')

            ax1.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

            ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax1.set_ylim((0, None))
            ax1.set_xlim((rt_min_now-2, rt_max_now+2))
            ax1.set_xlabel('Time (min)',fontsize=9)
            ax1.set_ylabel('Relative Intensity',fontsize=9)
            plt.xticks(size=9)
            plt.yticks(size=9)

            text = 'Index: ' + str(f_num + 1)
            ax1.text((rt_min_now + rt_max_now) / 2, 1.1 * max(inten), text, horizontalalignment='center',
                     verticalalignment='top',
                     color='red',
                     fontsize=8)

            text='Name: '+ str(f_name)
            ax1.text(rt_min_now-1, 0.98*max(inten), text, horizontalalignment='center', verticalalignment='top', color='blue',
                    fontsize=8)
            text = 'm/ztheor: ' + str(f_mztheor)
            ax1.text(rt_min_now-1, 0.9* max(inten), text, horizontalalignment='center', verticalalignment='top',
                     color='blue',
                     fontsize=8)
            text ='tRtheor: ' + str(f_tRtheor) + 'min'
            ax1.text(rt_min_now-1, 0.82* max(inten), text, horizontalalignment='center', verticalalignment='top',
                     color='blue',
                     fontsize=8)

            text = 'tR: ' + str(f_tRexper) + 'min'
            ax1.text(0.5*(rt_min_now+rt_max_now)+0.6, max(inten), text, horizontalalignment='center', verticalalignment='top',
                     color='red',
                     fontsize=8)

            feature.fig1.subplots_adjust(left=0.2, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
            feature.canvas1.draw()


            f_MSMSexper=f_MSMSexper.split(';')
            X1 = []
            Y1 = []

            for l in range(len(f_MSMSexper)):
                X1.append(float((f_MSMSexper[l].split(' ')[0])))
                Y1.append((float((f_MSMSexper[l].split(' ')[1])) * 100))

            f_MSMStheor = f_MSMStheor.split(';')
            X2 = []
            Y2 = []

            for l in range(len(f_MSMStheor)):
                X2.append(float((f_MSMStheor[l].split(' ')[0])))
                Y2.append((float((f_MSMStheor[l].split(' ')[1])) * (-100)))

            feature.fig2.clear()
            ax2 = feature.fig2.add_subplot(111)
            # 调整图像大小
            ax2.cla()
            msms_max = max(X1)
            ax2.bar(X1, Y1, width=msms_max / 150, color="red")
            ax2.bar(X2, Y2, width=msms_max / 150, color="green")
            ax2.plot([50, msms_max+10], [-0.05, -0.05], color="black", linewidth=1)
            ax2.set_ylim((-110, 120))
            ax2.set_xlim((50, msms_max+10))

            index_max = [x[0] for x in sorted(enumerate(Y1), key=lambda x: x[1])[-4:]]
            for y, x in enumerate(X1):
                if y in index_max:
                    plt.text(x, Y1[y] + 0.5, "%s" % round(x, 4), ha="center", va="bottom",fontsize=7)
            ax2.set_xlabel('m/z',fontsize=9)
            plt.xticks(size=9)
            plt.yticks(size=9)

            text = 'SIM: ' + str(f_sim)
            ax2.text((min(X1)+msms_max)/2, 130,text, horizontalalignment='center',verticalalignment='top',color='red',fontsize=8)

            feature.fig2.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
            feature.canvas2.draw()
        except:
            try:
                f_name = feature.Name[f_num]
                f_tRtheor = feature.tRtheor[f_num]
                f_tRexper = feature.tRMS[f_num]
                f_EIC = feature.EIC[f_num]
                f_MSMSexper = feature.MSMSexper[f_num]
                f_mztheor = feature.mztheor[f_num]
                f_sim_num = feature.SIM_num[f_num]
                f_same_Fig=feature.Same_Fig[f_num]
                f_same_Fig=f_same_Fig.strip(',').split(',')

                f_EIC = f_EIC.split(';')
                inten = []
                rt = []
                for w in range(len(f_EIC)):
                    rt.append((float(f_EIC[w].split(' ')[0])) / 60)
                    inten.append(float(f_EIC[w].split(' ')[1]))

                feature.fig1.clear()
                ax1 = feature.fig1.add_subplot(111)
                # 调整图像大小
                ax1.cla()
                rt_max_now = max(rt)
                rt_min_now = min(rt)

                ax1.plot(rt, inten, linewidth=0.8, color='black')

                ax1.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

                ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                ax1.set_ylim((0, None))
                ax1.set_xlim((rt_min_now - 2, rt_max_now + 2))
                ax1.set_xlabel('Time (min)', fontsize=9)
                ax1.set_ylabel('Relative Intensity', fontsize=9)
                plt.xticks(size=9)
                plt.yticks(size=9)

                text = 'Index: ' + str(f_num + 1)
                ax1.text((rt_min_now+rt_max_now)/2, 1.1 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='red',
                         fontsize=8)

                text = 'Name: ' + str(f_name)
                ax1.text(rt_min_now - 1, 0.98 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='blue',
                         fontsize=8)
                text = 'm/ztheor: ' + str(f_mztheor)
                ax1.text(rt_min_now - 1, 0.9 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='blue',
                         fontsize=8)
                text = 'tRtheor: ' + str(f_tRtheor) + 'min'
                ax1.text(rt_min_now - 1, 0.82 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='blue',
                         fontsize=8)

                text = 'tR: ' + str(f_tRexper) + 'min'
                ax1.text(0.5 * (rt_min_now + rt_max_now) + 0.6, max(inten), text, horizontalalignment='center',
                         verticalalignment='top',
                         color='red',
                         fontsize=8)

                feature.fig1.subplots_adjust(left=0.2, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
                feature.canvas1.draw()


                f_MSMSexper = f_MSMSexper.split(';')
                X1 = []
                Y1 = []
                X2 = []
                Y2 = []
                for l in range(len(f_MSMSexper)):
                    X1.append(float((f_MSMSexper[l].split(' ')[0])))
                    Y1.append((float((f_MSMSexper[l].split(' ')[1])) * 100))
                    for m in range(len(f_same_Fig)):
                        if abs(float((f_MSMSexper[l].split(' ')[0]))-float(f_same_Fig[m]))<0.01:
                            X2.append(float(f_MSMSexper[l].split(' ')[0]))
                            Y2.append(float(f_MSMSexper[l].split(' ')[1]) * 100)

                feature.fig2.clear()
                ax2 = feature.fig2.add_subplot(111)
                # 调整图像大小
                ax2.cla()
                msms_max = max(X1)
                ax2.bar(X1, Y1, width=msms_max / 150, color="dimgray")
                ax2.bar(X2, Y2, width=msms_max / 150, color="red")
                ax2.set_ylim((0, 120))
                ax2.set_xlim((50, msms_max + 10))
                index_max = [x[0] for x in sorted(enumerate(Y1), key=lambda x: x[1])[-4:]]
                for y, x in enumerate(X2):
                    if y in index_max:
                        plt.text(x, Y2[y] + 0.5, "%s" % round(x, 4), ha="center", va="bottom", fontsize=7)
                ax2.set_xlabel('m/z', fontsize=9)
                plt.xticks(size=9)
                plt.yticks(size=9)

                text = 'Same fragments: ' + str(f_sim_num)
                ax2.text((min(X1)+msms_max)/2, 130, text, horizontalalignment='center', verticalalignment='top', color='red',
                         fontsize=8)

                feature.fig2.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
                feature.canvas2.draw()
            except:
                f_tRexper = feature.tR[f_num]
                f_EIC = feature.EIC[f_num]
                f_MSMSexper = feature.MSMSexper[f_num]
                f_mzexper = feature.mzexper[f_num]
                f_cfi_num = feature.cfi_num[f_num]
                f_EIC = f_EIC.split(';')
                inten = []
                rt = []
                for w in range(len(f_EIC)):
                    rt.append((float(f_EIC[w].split(' ')[0])) / 60)
                    inten.append(float(f_EIC[w].split(' ')[1]))

                feature.fig1.clear()
                ax1 = feature.fig1.add_subplot(111)
                # 调整图像大小
                ax1.cla()
                rt_max_now = max(rt)
                rt_min_now = min(rt)
                ax1.plot(rt, inten, linewidth=0.8, color='black')

                ax1.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

                ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                ax1.set_ylim((0, None))
                ax1.set_xlim((rt_min_now - 2, rt_max_now + 2))
                ax1.set_xlabel('Time (min)', fontsize=9)
                ax1.set_ylabel('Relative Intensity', fontsize=9)
                plt.xticks(size=9)
                plt.yticks(size=9)

                text = 'Index: ' + str(f_num + 1)
                ax1.text((rt_min_now + rt_max_now) / 2, 1.1 * max(inten), text, horizontalalignment='center',
                         verticalalignment='top',
                         color='red',
                         fontsize=8)

                text = 'm/z: ' + str(f_mzexper)
                ax1.text(rt_min_now - 1, 0.98 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='blue',
                         fontsize=8)
                text = 'tR: ' + str(f_tRexper) + 'min'
                ax1.text(rt_min_now - 1, 0.9 * max(inten), text, horizontalalignment='center', verticalalignment='top',
                         color='blue',
                         fontsize=8)


                feature.fig1.subplots_adjust(left=0.2, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
                feature.canvas1.draw()

                f_MSMSexper = f_MSMSexper.split(';')
                X1 = []
                Y1 = []
                for l in range(len(f_MSMSexper)):
                    X1.append(float((f_MSMSexper[l].split(' ')[0])))
                    Y1.append((float((f_MSMSexper[l].split(' ')[1])) * 100))


                feature.fig2.clear()
                ax2 = feature.fig2.add_subplot(111)
                # 调整图像大小
                ax2.cla()
                msms_max = max(X1)
                ax2.bar(X1, Y1, width=msms_max / 150, color="black")
                ax2.set_ylim((0, 120))
                ax2.set_xlim((50, msms_max + 10))
                ax2.set_xlabel('m/z', fontsize=9)
                plt.xticks(size=9)
                plt.yticks(size=9)

                text = 'CFI&CNL: ' + str(f_cfi_num)
                ax2.text((min(X1)+msms_max)/2, 130, text, horizontalalignment='center', verticalalignment='top', color='red',
                         fontsize=8)
                feature.fig2.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.92, wspace=None, hspace=None)
                feature.canvas2.draw()

    def Image_In_ProcessII(self, feature):
        feature.fig1.clear()
        feature.canvas1.draw()
        feature.fig2.clear()
        feature.label_2.setText('')
        feature.canvas2.draw()
class MyMainForm(QMainWindow, Ui_widget):

    def __init__(self):
        super().__init__()
        self.setupUi(self)  # 初始化窗体设置
        qr = self.frameGeometry()  #显示至页面中心
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        self.setWindowFlags(Qt.WindowMinimizeButtonHint | Qt.WindowCloseButtonHint)  # 显示关闭按钮和最小化按钮
        self.setWindowFlags(Qt.MSWindowsFixedSizeDialogHint)
        self.tableWidget_2.setColumnCount(1)
        self.tableWidget_2.setStyleSheet(
            "QHeaderView::section{background:#6FACFF;color: black;font:14pt \'Times New Roman\';color: black;}")
        self.tableWidget_2.setLineWidth(0)
        self.tableWidget_2.setHorizontalHeaderLabels([
            'Category'])
        self.tableWidget_2.setColumnWidth(0, 230)
        self.tableWidget.setStyleSheet(
            "QHeaderView::section{background:#6FACFF;color: black;font:14pt \'Times New Roman\';color: black;}")
        self.tableWidget.setLineWidth(0)

        self.pushButton.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_3.setCursor(QCursor(Qt.PointingHandCursor))  # 鼠标移动至按钮时，更改形状为手状
        self.pushButton_4.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_5.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_6.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_7.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_9.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_2.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_12.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_13.setCursor(QCursor(Qt.PointingHandCursor))
        self.pushButton_11.setCursor(QCursor(Qt.PointingHandCursor))



        self.pushButton.clicked.connect(self.getmzml_1)
        self.pushButton_3.clicked.connect(self.getmzml_2)
        self.pushButton_4.clicked.connect(self.getExcel)
        self.pushButton_5.clicked.connect(self.standard_suspectscreening)
        self.pushButton_6.clicked.connect(self.online_literature_suspectscreening)
        self.pushButton_7.clicked.connect(self.tps_nontargetscreening)
        self.pushButton_9.clicked.connect(self.cfi_nontargetscreening)

        self.pushButton_2.clicked.connect(self.getresult)

        self.pushButton_12.clicked.connect(self.ShowLast)
        self.pushButton_13.clicked.connect(self.ShowNext)
        self.pushButton_11.clicked.connect(self.Export)


        self.mzml_path_1 = False
        self.peaklistpath_1=False
        self.index_no = 0
        self.current_path = os.getcwd()
        self.standard_database = pd.read_excel('supports/standard_LCMS.xlsx')
        self.online_database = pd.read_excel('supports/online_LCMS.xlsx')
        self.literature_database = pd.read_excel('supports/literature_LCMS.xlsx')
        self.tps_database = pd.read_excel('supports/tps_LCMS.xlsx')

        self.cfi_database = pd.read_excel('supports/cfi_LCMS.xlsx')

        self.class_list=list(self.cfi_database['Category'])


        self.tableWidget_2.setColumnCount(1)
        self.tableWidget_2.setHorizontalHeaderLabels(['Category'])

        self.tableWidget_2.setRowCount(len(self.class_list))
        for i in range(len(self.class_list)):
            new_item = QTableWidgetItem(str(self.class_list[i]))
            self.tableWidget_2.setItem(i, 0, new_item)

        self.Process = DrawFunction()  # process对象包含了所有的信号处理函数及其画图
        self.ImageLayout()
    def Export(self):
        if  0<=self.index_no<=self.resultdf.shape[0]-1:
            self.resultdf.loc[self.index_no, 'Manual check'] = self.lineEdit_16.text()
            os.chdir(self.current_path)
            os.chdir('results')
            file_name = 'Manualcheck_'+str(self.resultpath.split('/')[-1])
            self.resultdf.to_excel(file_name, index=False)

    def ImageLayout(self):
        self.fig1 = plt.figure()
        self.canvas1 = FigureCanvas(self.fig1)
        layout1 = QVBoxLayout()  # 垂直布局
        layout1.addWidget(self.canvas1)
        self.graphicsView.setLayout(layout1)  # 设置好布局之后调用函数

        self.fig2 = plt.figure()
        self.canvas2 = FigureCanvas(self.fig2)
        layout2 = QVBoxLayout()  # 垂直布局
        layout2.addWidget(self.canvas2)
        self.graphicsView_2.setLayout(layout2)  # 设置好布局之后调用函数

    def getmzml_1(self):
        try:
            self.mzml_path_1 = QFileDialog.getExistingDirectory(None, "选择路径", os.getcwd())

        except Exception as e:
            print(e)
    def getmzml_2(self):
        try:
            self.mzml_path_2 = QFileDialog.getExistingDirectory(None, "选择路径", os.getcwd())
        except Exception as e:
            print(e)

    def getExcel(self):
        try:
            # 选择的文件绝对路径
            self.peaklistpath_1 = QFileDialog.getOpenFileName(None, "选择文件", os.getcwd(), 'excel(*.xlsx)')[0]
        except Exception as e:
            print(e)

    def standard_suspectscreening(self):

        ms1_error = float(self.lineEdit.text())
        rt_error = float(self.lineEdit_2.text())
        isotope_sc = 1-float(self.lineEdit_4.text())
        SIM_sc = float(self.lineEdit_3.text())
        if self.mzml_path_1 ==False :
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select the path of the mzML folder',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.mzml_path_1 != False:
            print(' ')
            print('......正在进行基于自建标样数据库匹配的可疑筛查......')
            print(' ')
            self.qthreading1 = ThreadOne()
            self.qthreading1.ms1_error = ms1_error
            self.qthreading1.rt_error = rt_error
            self.qthreading1.isotope_sc = isotope_sc
            self.qthreading1.SIM_sc = SIM_sc
            self.qthreading1.mzml_path_1 = self.mzml_path_1
            self.qthreading1.current_path = self.current_path
            self.qthreading1.standard_database = self.standard_database
            self.qthreading1.start()
            self.qthreading1.wait()
            print(' ')
            print('......已完成基于自建标样数据库匹配的可疑筛查......')
            print(' ')


    def online_literature_suspectscreening(self):

        if self.mzml_path_1 ==False :
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select the path of the mzML folder',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.mzml_path_1 != False:
            print(' ')
            print('......正在进行基于公共数据库匹配的可疑筛查......')
            print(' ')
            ms1_error = float(self.lineEdit_5.text())
            rt_error = float(self.lineEdit_6.text())
            isotope_sc = 1-float(self.lineEdit_7.text())
            SIM_online = float(self.lineEdit_8.text())
            frgnum_literature = float(self.lineEdit_9.text())
            self.qthreading2 = ThreadTwo()
            self.qthreading2.ms1_error = ms1_error
            self.qthreading2.rt_error = rt_error
            self.qthreading2.isotope_sc = isotope_sc
            self.qthreading2.SIM_online = SIM_online
            self.qthreading2.frgnum_literature = frgnum_literature
            self.qthreading2.mzml_path_1 = self.mzml_path_1
            self.qthreading2.current_path = self.current_path
            self.qthreading2.online_database = self.online_database
            self.qthreading2.literature_database = self.literature_database
            self.qthreading2.start()
            self.qthreading2.wait()
            print(' ')
            print('......已完成基于公共数据库匹配的可疑筛查......')
            print(' ')

    def tps_nontargetscreening(self):

        ms1_error = float(self.lineEdit_10.text())
        rt_error = float(self.lineEdit_12.text())
        isotope_sc = float(self.lineEdit_11.text())
        frgnum_cfmid = float(self.lineEdit_13.text())
        if self.peaklistpath_1 ==False :
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select a peak table file',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.mzml_path_2 ==False :
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select the path of the mzML folder',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.peaklistpath_1 != False and self.mzml_path_2 !=False:
            print(' ')
            print('......正在进行基于转化产物匹配的非靶向筛查......')
            print(' ')
            self.qthreading3 = ThreadThree()
            self.qthreading3.ms1_error = ms1_error
            self.qthreading3.rt_error = rt_error
            self.qthreading3.isotope_sc = isotope_sc
            self.qthreading3.frgnum_cfmid = frgnum_cfmid

            self.qthreading3.mzml_path_2 = self.mzml_path_2
            self.qthreading3.current_path = self.current_path
            self.qthreading3.tps_database = self.tps_database
            self.qthreading3.peaklistpath_1 = self.peaklistpath_1
            self.qthreading3.start()
            self.qthreading3.wait()
            print(' ')
            print('......已完成基于转化产物匹配的非靶向筛查......')
            print(' ')




    def cfi_nontargetscreening(self):

        ms2_error = float(self.lineEdit_14.text())
        inten_filter = float(self.lineEdit_15.text())
        if self.peaklistpath_1 == False:
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select a peak table file',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.mzml_path_2 == False:
            reply = QMessageBox.question(self, 'Warning',
                                         'Please select the path of the mzML folder',
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.Yes)
            if reply == QMessageBox.Yes:
                pass
        if self.peaklistpath_1 != False and self.mzml_path_2 != False:
            print(' ')
            print('......正在进行基于特征碎片匹配的非靶向筛查......')
            print(' ')
            self.qthreading4 = ThreadFour()
            self.qthreading4.ms2_error = ms2_error
            self.qthreading4.inten_filter = inten_filter

            self.qthreading4.mzml_path_2 = self.mzml_path_2
            self.qthreading4.current_path = self.current_path
            self.qthreading4.cfi_database = self.cfi_database
            self.qthreading4.peaklistpath_1 = self.peaklistpath_1
            self.qthreading4.start()
            self.qthreading4.wait()
            print(' ')
            print('......已完成基于特征碎片匹配的非靶向筛查......')
            print(' ')

    def getresult(self):
        self.index_no = 0
        self.resultpath = QFileDialog.getOpenFileName(None, "选择文件", os.getcwd(), 'excel(*.xlsx)')[0]
        self.resultdf=pd.read_excel(self.resultpath)
        print(' ')
        print('......正在人工审核定性结果......')
        print(' ')
        print(self.resultpath)
        try:
            self.lineEdit_16.setText(str(self.resultdf.loc[0, 'Manual check']))
        except:
            self.resultdf['Manual check']=''
            self.lineEdit_16.setText('')
        try:
            self.Name= self.resultdf['Name']
            self.Formula = self.resultdf['Formula']
            self.mztheor = self.resultdf['m/ztheor']
            self.tRtheor = self.resultdf['tR/mintheor']
            self.tRMS = self.resultdf['tR/min-MS']
            self.tRMSMS = self.resultdf['tR/min-MSMS']
            self.SIM = self.resultdf['SIM']
            self.CL= self.resultdf['Confidence level']
            self.EIC = self.resultdf['EIC']
            self.MSMStheor = self.resultdf['MS/MStheor']
            self.MSMSexper = self.resultdf['MS/MSexper']

            self.tableWidget.setColumnCount(7)
            uniform_width = 70  # 你可以根据需要调整这个值
            for column in range(self.tableWidget.columnCount()):
                self.tableWidget.setColumnWidth(column, uniform_width)
            self.tableWidget.setHorizontalHeaderLabels(['Name','Formula','tR_theor','tR_exp', 'm/z','SIM','CL'])
            # 设置水平表头标签左对齐
            for column in range(self.tableWidget.columnCount()):
                header_item = self.tableWidget.horizontalHeaderItem(column)
                if header_item is not None:
                    header_item.setTextAlignment(Qt.AlignLeft)
            self.tableWidget.setColumnWidth(0, 100)
            self.tableWidget.setColumnWidth(1, 100)
            self.tableWidget.setRowCount(len(self.Name))
            for i in range(len(self.Name)):
                new_item0 = QTableWidgetItem(str(self.Name[i]))
                self.tableWidget.setItem(i, 0, new_item0)

                new_item1 = QTableWidgetItem(str(self.Formula[i]))
                self.tableWidget.setItem(i, 1, new_item1)

                new_item2 = QTableWidgetItem(str(self.tRtheor[i]))
                self.tableWidget.setItem(i, 2, new_item2)

                new_item3 = QTableWidgetItem(str(self.tRMS[i]))
                self.tableWidget.setItem(i, 3, new_item3)

                new_item4 = QTableWidgetItem(str(self.mztheor[i]))
                self.tableWidget.setItem(i, 4, new_item4)

                new_item5 = QTableWidgetItem(str(self.SIM[i]))
                self.tableWidget.setItem(i, 5, new_item5)

                new_item6 = QTableWidgetItem(str(self.CL[i]))
                self.tableWidget.setItem(i, 6, new_item6)
            self.Process.Image_In_Process(self)
        except:
            try:
                self.Name = self.resultdf['Name']
                self.Formula = self.resultdf['Formula']
                self.mztheor = self.resultdf['m/ztheor']
                self.tRtheor = self.resultdf['tR/mintheor']
                self.tRMS = self.resultdf['tR/min-MS']
                self.tRMSMS = self.resultdf['tR/min-MSMS']
                self.SIM_num = self.resultdf['SIM-num']
                self.CL = self.resultdf['Confidence level']
                self.EIC = self.resultdf['EIC']
                self.Same_Fig = self.resultdf['SIM-Frg']
                self.MSMSexper = self.resultdf['MS/MSexper']

                self.tableWidget.setColumnCount(7)
                uniform_width = 70
                for column in range(self.tableWidget.columnCount()):
                    self.tableWidget.setColumnWidth(column, uniform_width)
                self.tableWidget.setHorizontalHeaderLabels(
                    ['Name', 'Formula', 'tR_theor', 'tR_exp', 'm/z', 'SIM_num', 'CL'])
                for column in range(self.tableWidget.columnCount()):
                    header_item = self.tableWidget.horizontalHeaderItem(column)
                    if header_item is not None:
                        header_item.setTextAlignment(Qt.AlignLeft)
                self.tableWidget.setColumnWidth(0, 100)
                self.tableWidget.setColumnWidth(1, 100)
                self.tableWidget.setRowCount(len(self.Name))
                for i in range(len(self.Name)):
                    new_item0 = QTableWidgetItem(str(self.Name[i]))
                    self.tableWidget.setItem(i, 0, new_item0)

                    new_item1 = QTableWidgetItem(str(self.Formula[i]))
                    self.tableWidget.setItem(i, 1, new_item1)

                    new_item2 = QTableWidgetItem(str(self.tRtheor[i]))
                    self.tableWidget.setItem(i, 2, new_item2)

                    new_item3 = QTableWidgetItem(str(self.tRMS[i]))
                    self.tableWidget.setItem(i, 3, new_item3)

                    new_item4 = QTableWidgetItem(str(self.mztheor[i]))
                    self.tableWidget.setItem(i, 4, new_item4)

                    new_item5 = QTableWidgetItem(str(self.SIM_num[i]))
                    self.tableWidget.setItem(i, 5, new_item5)

                    new_item6 = QTableWidgetItem(str(self.CL[i]))
                    self.tableWidget.setItem(i, 6, new_item6)
                self.Process.Image_In_Process(self)
            except:
                self.Name = self.resultdf['Scan']
                self.mzexper = self.resultdf['m/zexper']
                self.tR = self.resultdf['tR/min']
                self.EIC = self.resultdf['EIC']
                self.MSMSexper = self.resultdf['MS/MSexper']
                self.cfi_num = self.resultdf['CFI_NUM']
                self.tableWidget.setColumnCount(5)
                uniform_width = 100
                for column in range(self.tableWidget.columnCount()):
                    self.tableWidget.setColumnWidth(column, uniform_width)
                self.tableWidget.setHorizontalHeaderLabels(
                    ['m/zexper', 'tR/min', 'EIC', 'MS/MSexper', 'CFI_NUM'])
                for column in range(self.tableWidget.columnCount()):
                    header_item = self.tableWidget.horizontalHeaderItem(column)
                    if header_item is not None:
                        header_item.setTextAlignment(Qt.AlignLeft)
                self.tableWidget.setRowCount(len(self.mzexper))
                for i in range(len(self.mzexper)):
                    new_item0 = QTableWidgetItem(str(self.mzexper[i]))
                    self.tableWidget.setItem(i, 0, new_item0)

                    new_item1 = QTableWidgetItem(str(self.tR[i]))
                    self.tableWidget.setItem(i, 1, new_item1)

                    new_item2 = QTableWidgetItem(str(self.EIC[i]))
                    self.tableWidget.setItem(i, 2, new_item2)

                    new_item3 = QTableWidgetItem(str(self.MSMSexper[i]))
                    self.tableWidget.setItem(i, 3, new_item3)

                    new_item4 = QTableWidgetItem(str(self.cfi_num[i]))
                    self.tableWidget.setItem(i, 4, new_item4)
                self.Process.Image_In_Process(self)

    def ShowLast(self):
        if self.index_no > 0:
            try:
                self.resultdf.loc[self.index_no, 'Manual check'] = self.lineEdit_16.text()
            except:
                pass
            self.index_no = self.index_no-1
            if self.index_no > -1:
                try:
                    self.lineEdit_16.setText(str(self.resultdf.loc[self.index_no, 'Manual check']))
                    self.Process.Image_In_Process(self)
                except:
                    pass


    def ShowNext(self):
        if self.index_no < self.resultdf.shape[0] - 1:
            try:
                self.resultdf.loc[self.index_no, 'Manual check'] = self.lineEdit_16.text()
            except:
                pass
            self.index_no = self.index_no+1
            if self.index_no <= self.resultdf.shape[0]-1:
                try:
                    self.lineEdit_16.setText(str(self.resultdf.loc[self.index_no, 'Manual check']))
                    self.Process.Image_In_Process(self)
                except:
                    pass



