#include "mainwindow.h"

#include "waveviewer.hpp"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    setWindowTitle("WaveServer by Marcel Urban");
    setCentralWidget(new WaveViewer);
}

MainWindow::~MainWindow()
{
}

